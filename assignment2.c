/* Gets the neighbors in a cartesian communicator
* Orginally written by Mary Thomas
* - Updated Mar, 2015
* Link: https://edoras.sdsu.edu/~mthomas/sp17.605/lectures/MPI-Cart-Comms-and-Topos.pdf
* Modifications to fix bugs, include an async send and receive and to revise print output
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>
#include <assert.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <netinet/in.h>
#include <net/if.h>
#include <arpa/inet.h>
#include <pthread.h>
#include <mqueue.h>


#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1
// #define RANDOM_MAX_WATER_DEPTH 7000
// #define RANDOM_MIN_WATER_DEPTH 5000
#define TOLERENCERANGE 500
#define ALTIMETERRANGE 500
#define MAX_ARRAY_COUNT 5
#define LIMIT_THRESHHOLD 6000
#define SENDING_INFO_TAG 10
#define TERMINATE_TAG 5
#define IP_TAG 20
#define VM_NAME_TAG 30
#define SEND_REQ 35
#define RECV_REQ 50
#define THRESHOLD_TAG 60
// defines, like some settings
#define WAVE_GENERATION_INTERVAL 5 // in seconds - each 5 seconds
#define WAVE_TRANSMISSION_INTERVAL 4 // in seconds - each 4 seconds (to get some alerts)
#define WAVE_TABLE_SIZE 200 // count of records for table, works as 'fifo'
#define WAVE_MESSAGE_QUE_TO_THREAD (const char*)"/wave_message_que_to" // default message que name for write to tread
#define WAVE_MESSAGE_QUE_FROM_THREAD (const char*)"/wave_message_que_from" // default message que name for read from tread
#define WAVE_APPLICATION_WORK_TIME 10 // in seconds - time before exit

int master_io(MPI_Comm world_comm, MPI_Comm comm);
int slave_io(MPI_Comm master_comm, MPI_Comm comm,float Threshold);
float random_water_depth(int my_rank, float Threshold);
void* ProcessFunc(void *pArg); // Common function prototype
float simple_moving_average(int count,float arr[], float new_water_depth);
void WriteToFile(char *pFilename, char *txt, bool newline);
void WriteNumToFile(char *pFilename, int num, bool newline);
void WriteFloatToFile(char *pFilename, float num, bool small, bool newline);


int ndims = 2;
int *dims;
int nrows, ncols;
float simple_moving_average_val;
int thread_termination = 1;
MPI_Comm comm2D;
float Threshold;

// This struct is used to send information to the base station from tsnuami sensor
struct valuestruct {
    int timeStamp;
    int rank[5];
    int my_cord[2];
    int neighbour1_cord[2];
    int neighbour2_cord[2];
    int neighbour3_cord[2];
    int neighbour4_cord[2];
    float height[5];
    char ipAddress[15];
    char processor_name [50];
    char ipAddressNeighbour1[15];
    char ipAddressNeighbour2[15];
    char ipAddressNeighbour3[15];
    char ipAddressNeighbour4[15];
    char processor_name_neighb1[50];
    char processor_name_neighb2[50];
    char processor_name_neighb3[50];
    char processor_name_neighb4[50];
} ;

// used types of data and structs

// declare type of signle wave record
typedef struct
{
  // keep date time stamp
  int year;
  int month;
  int day;
  int hour;
  int minute;
  int second;
  // keep coordinates
  int x;
  int y;
  // keep wave height
  int height;
} wave_record_t;

// declare type of records table
typedef struct
{
  // waves array
  wave_record_t* records;
  // pointer for reading
  int pos_read;
  // pointer for writing
  int pos_write;
  // count of items
  int count;
  // size of table
  int size;
  //  mutex to protect wave table from multithread access
  pthread_mutex_t mutex;
} wave_table_t;

// devlare instance of used waves table - global for all
wave_table_t wave_table;


int main(int argc, char *argv[]) {

	int size, my_rank, ierr, fd;
    MPI_Comm new_comm;
    struct ifreq ifr;


	/* start up initial MPI environment */
	//MPI_Init(&argc, &argv);
  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Take in command line arguments
	if (argc == 4) {
		nrows = atoi (argv[1]);
		ncols = atoi (argv[2]);
    Threshold = atoi (argv[3]);
	}
  else{
    printf("Please enter the arguments to form cartesian grid  in the form nrows ncols Threshold\n");
    MPI_Finalize();
    return 0;
  }


  // Splitting the communicator so that the first node represnets the base station and all other nodes represent the tsunami sensors
  MPI_Comm_split( MPI_COMM_WORLD,my_rank == 0, 0, &new_comm);


  if (my_rank == 0) {
      master_io(MPI_COMM_WORLD, new_comm);
  }
  else{
      slave_io(MPI_COMM_WORLD, new_comm,Threshold);
  }
  MPI_Finalize();
  return 0;
}

/* Base station thread function to receive alerts from sensors*/
void* receiveMsg(void *pArg) 
{
    struct valuestruct alert;
    int sensorCoords[2];
    bool coordFound = false;
    bool match = false;
    clock_t timer;
    char str_time_ReportingNode[100];
    char str_date_ReportingNode[100];
    char str_time_BaseStation[100];
    char str_date_BaseStation[100];
    float waterHeight;
    int index;
    time_t reportedNodeTime;
    time_t baseStationTime;
    struct tm *tm;
    int numAdjacentMatches = 0;
    float altimeterHeight;
    int altimeterTimestamp;
    MPI_Datatype ReceivedType;
    MPI_Datatype type[18] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR};
    int blocklen[18] = { 1, 5, 2, 2, 2, 2, 2, 5, 15, 50, 15, 15, 15, 15, 50, 50, 50, 50};
    MPI_Aint disp[18];
    MPI_Status status;
    MPI_Request request;
    char* logFile = "Log.txt";
    int numMatched = 0;
    int numMismatched = 0;
    float totalCommTime = 0;

    MPI_Get_address(&alert.timeStamp, &disp[0]);
    MPI_Get_address(&alert.rank, &disp[1]);
    MPI_Get_address(&alert.my_cord, &disp[2]);
    MPI_Get_address(&alert.neighbour1_cord, &disp[3]);
    MPI_Get_address(&alert.neighbour2_cord, &disp[4]);
    MPI_Get_address(&alert.neighbour3_cord, &disp[5]);
    MPI_Get_address(&alert.neighbour4_cord, &disp[6]);
    MPI_Get_address(&alert.height, &disp[7]);
    MPI_Get_address(&alert.ipAddress, &disp[8]);
    MPI_Get_address(&alert.processor_name, &disp[9]);
    MPI_Get_address(&alert.ipAddressNeighbour1, &disp[10]);
    MPI_Get_address(&alert.ipAddressNeighbour2, &disp[11]);
    MPI_Get_address(&alert.ipAddressNeighbour3, &disp[12]);
    MPI_Get_address(&alert.ipAddressNeighbour4, &disp[13]);
    MPI_Get_address(&alert.processor_name_neighb1, &disp[14]);
    MPI_Get_address(&alert.processor_name_neighb2, &disp[15]);
    MPI_Get_address(&alert.processor_name_neighb3, &disp[16]);
    MPI_Get_address(&alert.processor_name_neighb4, &disp[17]);

    //Make relative
    disp[1]=disp[1]-disp[0];
    disp[2]=disp[2]-disp[0];
    disp[3]=disp[3]-disp[0];
    disp[4]=disp[4]-disp[0];
    disp[5]=disp[5]-disp[0];
    disp[6]=disp[6]-disp[0];
    disp[7]=disp[7]-disp[0];
    disp[8]=disp[8]-disp[0];
    disp[9]=disp[9]-disp[0];
    disp[10]=disp[10]-disp[0];
    disp[11]=disp[11]-disp[0];
    disp[12]=disp[12]-disp[0];
    disp[13]=disp[13]-disp[0];
    disp[14]=disp[14]-disp[0];
    disp[15]=disp[15]-disp[0];
    disp[16]=disp[16]-disp[0];
    disp[17]=disp[17]-disp[0];
    disp[0]=0;

    // Create MPI struct
    MPI_Type_create_struct(18, blocklen, disp, type, &ReceivedType);
    MPI_Type_commit(&ReceivedType);

	int i = 0, nsensors, size, numIterations = 5;
	int iter = 0;
	char buf[256], buf2[256];
	MPI_Comm_size(MPI_COMM_WORLD, &size );
	
	int* p = (int*)pArg;
	// nsensors = *p;

    for (iter = 0; iter < numIterations; iter++){
        MPI_Recv(&alert, 1, ReceivedType, MPI_ANY_SOURCE, SENDING_INFO_TAG, MPI_COMM_WORLD, &status);

        timer = clock();

        sensorCoords[0] = alert.my_cord[0];
        sensorCoords[1] = alert.my_cord[1];

        waterHeight = alert.height[0];

        reportedNodeTime = (time_t) alert.timeStamp;
        tm = localtime(&reportedNodeTime);

        strftime(str_time_ReportingNode, sizeof(str_time_ReportingNode), "%T", tm);
        strftime(str_date_ReportingNode, sizeof(str_date_ReportingNode), "%a %F", tm);

        // baseStationTime = time(NULL);
        time(&baseStationTime);
        tm = localtime(&baseStationTime);

        strftime(str_time_BaseStation, sizeof(str_time_BaseStation), "%T", tm);
        strftime(str_date_BaseStation, sizeof(str_date_BaseStation), "%a %F", tm);

       // Look for sensorCoords in altimeter data structure.
       for(i = 0; i < wave_table.size; i++){
           if (wave_table.records[i].x == sensorCoords[0]){
               if (wave_table.records[i].y == sensorCoords[1]){
                   coordFound = true;
                   index = i;
                   break;
               }
           }
       }

        if (coordFound){
            altimeterHeight = wave_table.records[index].height/1000;

            if ((waterHeight > (altimeterHeight - ALTIMETERRANGE)) && (waterHeight < (altimeterHeight + ALTIMETERRANGE))){
                match = true;
            }
        }

        // Data to be printed
        WriteToFile(logFile, "-----------------------------------------------", true);
        WriteToFile(logFile, "Iteration: ", false);
        WriteNumToFile(logFile, (iter + 1), true);

        WriteToFile(logFile, "Logged time: ", false);
        WriteToFile(logFile, str_date_BaseStation, false);
        WriteToFile(logFile, " ", false);
        WriteToFile(logFile, str_time_BaseStation, true);

        WriteToFile(logFile, "Alert reported time: ", false);
        WriteToFile(logFile, str_date_ReportingNode, false);
        WriteToFile(logFile, " ", false);
        WriteToFile(logFile, str_time_ReportingNode, true);

        WriteToFile(logFile, "Alert type: ", false);

        if (!coordFound || !match){
            WriteToFile(logFile, "Mismatch", true);
            WriteToFile(logFile, "", true);
            numMismatched += 1;
        }
        else {
            WriteToFile(logFile, "Match", true);
            WriteToFile(logFile, "", true);
            numMatched += 1;
        }

        WriteToFile(logFile, "Reporting Node     Coord     Height (m)     IPv4", true);
        WriteNumToFile(logFile, alert.rank[0], false);
        WriteToFile(logFile, "                  (", false);
        WriteNumToFile(logFile, alert.my_cord[0], false);
        WriteToFile(logFile, ",", false);
        WriteNumToFile(logFile, alert.my_cord[1], false);
        WriteToFile(logFile, ")     ", false);

        WriteFloatToFile(logFile, alert.height[0], true, false);
        WriteToFile(logFile, "        ", false);
        WriteToFile(logFile, alert.ipAddress, false);
        WriteToFile(logFile, " (", false);
        WriteToFile(logFile, alert.processor_name, false);
        WriteToFile(logFile, ")\n", true);

        WriteToFile(logFile, "Adjacent Nodes     Coord     Height (m)     IPv4", true);
        if ((alert.rank[1] != -2) && (alert.height[1] > 0)){
            numAdjacentMatches += 1;

            WriteNumToFile(logFile, alert.rank[1], false);
            WriteToFile(logFile, "                  (", false);
            WriteNumToFile(logFile, alert.neighbour1_cord[0], false);
            WriteToFile(logFile, ",", false);
            WriteNumToFile(logFile, alert.neighbour1_cord[1], false);
            WriteToFile(logFile, ")     ", false);

            WriteFloatToFile(logFile, alert.height[1], true, false);
            WriteToFile(logFile, "        ", false);
            WriteToFile(logFile, alert.ipAddressNeighbour1, false);
            WriteToFile(logFile, " (", false);
            WriteToFile(logFile, alert.processor_name_neighb1, false);
            WriteToFile(logFile, ")", true);
        }

        if ((alert.rank[2] != -2) && (alert.height[2] > 0)){
            numAdjacentMatches += 1;

            WriteNumToFile(logFile, alert.rank[2], false);
            WriteToFile(logFile, "                  (", false);
            WriteNumToFile(logFile, alert.neighbour2_cord[0], false);
            WriteToFile(logFile, ",", false);
            WriteNumToFile(logFile, alert.neighbour2_cord[1], false);
            WriteToFile(logFile, ")     ", false);

            WriteFloatToFile(logFile, alert.height[2], true, false);
            WriteToFile(logFile, "        ", false);
            WriteToFile(logFile, alert.ipAddressNeighbour2, false);
            WriteToFile(logFile, " (", false);
            WriteToFile(logFile, alert.processor_name_neighb2, false);
            WriteToFile(logFile, ")", true);
        }

        if ((alert.rank[3] != -2) && (alert.height[3] > 0)){
            numAdjacentMatches += 1;

            WriteNumToFile(logFile, alert.rank[3], false);
            WriteToFile(logFile, "                  (", false);
            WriteNumToFile(logFile, alert.neighbour3_cord[0], false);
            WriteToFile(logFile, ",", false);
            WriteNumToFile(logFile, alert.neighbour3_cord[1], false);
            WriteToFile(logFile, ")     ", false);

            WriteFloatToFile(logFile, alert.height[3], true, false);
            WriteToFile(logFile, "        ", false);
            WriteToFile(logFile, alert.ipAddressNeighbour3, false);
            WriteToFile(logFile, " (", false);
            WriteToFile(logFile, alert.processor_name_neighb3, false);
            WriteToFile(logFile, ")", true);
        }

        if ((alert.rank[4] != -2) && (alert.height[4] > 0)){
            numAdjacentMatches += 1;

            WriteNumToFile(logFile, alert.rank[4], false);
            WriteToFile(logFile, "                  (", false);
            WriteNumToFile(logFile, alert.neighbour4_cord[0], false);
            WriteToFile(logFile, ",", false);
            WriteNumToFile(logFile, alert.neighbour4_cord[1], false);
            WriteToFile(logFile, ")     ", false);

            WriteFloatToFile(logFile, alert.height[4], true, false);
            WriteToFile(logFile, "        ", false);
            WriteToFile(logFile, alert.ipAddressNeighbour4, false);
            WriteToFile(logFile, " (", false);
            WriteToFile(logFile, alert.processor_name_neighb4, false);
            WriteToFile(logFile, ")", true);
        }

        WriteToFile(logFile, "\nHas satellite altimeter recorded reading for node ", false);
        WriteNumToFile(logFile, alert.rank[0], false);
        WriteToFile(logFile, ": ", false);

        if (coordFound){
            WriteToFile(logFile, "Yes", true);
            WriteToFile(logFile, "Satellite altimeter reporting time: ", false);
            WriteNumToFile(logFile, wave_table.records[index].year, false);
            WriteToFile(logFile, "-", false);
            WriteNumToFile(logFile, wave_table.records[index].month, false);
            WriteToFile(logFile, "-", false);
            WriteNumToFile(logFile, wave_table.records[index].day, false);
            WriteToFile(logFile, " ", false);
            WriteNumToFile(logFile, wave_table.records[index].hour, false);
            WriteToFile(logFile, ":", false);
            WriteNumToFile(logFile, wave_table.records[index].minute, false);
            WriteToFile(logFile, ":", false);
            WriteNumToFile(logFile, wave_table.records[index].second, true);

            WriteToFile(logFile, "Satellite altimeter reporting height (m): ", false);
            WriteFloatToFile(logFile, altimeterHeight, true, true);
            WriteToFile(logFile, "Satellite altimeter reporting Coord: (", false);
            WriteNumToFile(logFile, alert.my_cord[0], false);
            WriteToFile(logFile, ", ", false);
            WriteNumToFile(logFile, alert.my_cord[1], false);
            WriteToFile(logFile, ")\n", true);
        }
        else{
            WriteToFile(logFile, "No\n", true);
        }

        timer = clock() - timer;
        double time_taken = ((double) timer) / CLOCKS_PER_SEC;
        totalCommTime += time_taken;

        WriteToFile(logFile, "Communication time (seconds): ", false);
        WriteFloatToFile(logFile, time_taken, false, true);
        WriteToFile(logFile, "Total messages sent between reporting node and base station: 1", true);
        WriteToFile(logFile, "Number of adjacent matches to reporting node: ", false);
        WriteNumToFile(logFile, numAdjacentMatches, true);
        numAdjacentMatches = 0;

        WriteToFile(logFile, "Max. tolerance range between nodes readings (m): ", false);
        WriteNumToFile(logFile, TOLERENCERANGE, true);
        WriteToFile(logFile, "Max. tolerance range between satellite altimeter and reporting node readings (m): ", false);
        WriteNumToFile(logFile, ALTIMETERRANGE, true);
        WriteToFile(logFile, "-----------------------------------------------", true);

        sleep(5);
    }

    WriteToFile(logFile, "-----------------------------------------------", true);
    WriteToFile(logFile, "                   SUMMARY\n", true);
    WriteToFile(logFile, "Total messages sent (excluding termination message): " , false);
    WriteNumToFile(logFile, numIterations, true);
    WriteToFile(logFile, "Number of matched alerts: ", false);
    WriteNumToFile(logFile, numMatched, true);
    WriteToFile(logFile, "Number of mismatched alerts: ", false);
    WriteNumToFile(logFile, numMismatched, true);
    WriteToFile(logFile, "Total communication time (seconds): ", false);
    WriteFloatToFile(logFile, totalCommTime, false, true);
    WriteToFile(logFile, "-----------------------------------------------", true);

    int termination_num = -1;
    for (i = 1; i < size; i++){
        MPI_Isend(&termination_num, 1, MPI_INT, i, TERMINATE_TAG, MPI_COMM_WORLD, &request);
    }
	return 0;
}


/* Base station simulation */
int master_io(MPI_Comm world_comm, MPI_Comm comm)
{
    char* logFile = "Log.txt";
    WriteToFile(logFile, "LOG FILE", true);

	int size, nsensors;
	MPI_Comm_size(world_comm, &size );
	nsensors = size - 1;

    pthread_t threadID_RecvMsg;
    pthread_create(&threadID_RecvMsg, 0, receiveMsg, &nsensors);
    altimeter_simulation();
    pthread_join(threadID_RecvMsg, NULL);


    return 0;
}



void* ProcessFunc(void *pArg){
    while (100){
        if (thread_termination == -2){
            break;
        }
        else{
            int message;
            MPI_Status probe_status;
            MPI_Iprobe(MPI_ANY_SOURCE,SEND_REQ,comm2D,&message,&probe_status);
            // printf("Entering Thread of \n");

            if (message){
                int rank = probe_status.MPI_SOURCE;
                int new_message;
                MPI_Recv(&new_message,1,MPI_INT,rank,SEND_REQ,comm2D,MPI_STATUS_IGNORE);
                // int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                MPI_Send(&simple_moving_average_val,1,MPI_FLOAT,rank,RECV_REQ,comm2D);
            }
        }
    }
    return 0;

}

int slave_io(MPI_Comm master_comm, MPI_Comm comm, float Threshold)
{
  // Initialising variables
	int size, my_rank, reorder, my_cart_rank, ierr,fd;
  int nbr_i_lo, nbr_i_hi;
	int nbr_j_lo, nbr_j_hi;
	int coord[ndims];
	int wrap_around[ndims];
  float water_column_array[5];
  int count_of_array = 0;

  // creating the struct for sending information for report from tsnumai sensors to base station
  struct valuestruct values;
  MPI_Datatype Valuetype;
  MPI_Datatype type[18] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR};
  int blocklen[18] = { 1, 5, 2, 2, 2, 2, 2, 5, 15, 50, 15, 15, 15, 15, 50, 50 , 50, 50};
  MPI_Aint disp[18];

  MPI_Get_address(&values.timeStamp, &disp[0]);
  MPI_Get_address(&values.rank, &disp[1]);
  MPI_Get_address(&values.my_cord, &disp[2]);
  MPI_Get_address(&values.neighbour1_cord, &disp[3]);
  MPI_Get_address(&values.neighbour2_cord, &disp[4]);
  MPI_Get_address(&values.neighbour3_cord, &disp[5]);
  MPI_Get_address(&values.neighbour4_cord, &disp[6]);
  MPI_Get_address(&values.height, &disp[7]);
  MPI_Get_address(&values.ipAddress, &disp[8]);
  MPI_Get_address(&values.processor_name, &disp[9]);
  MPI_Get_address(&values.ipAddressNeighbour1, &disp[10]);
  MPI_Get_address(&values.ipAddressNeighbour2, &disp[11]);
  MPI_Get_address(&values.ipAddressNeighbour3, &disp[12]);
  MPI_Get_address(&values.ipAddressNeighbour4, &disp[13]);
  MPI_Get_address(&values.processor_name_neighb1, &disp[14]);
  MPI_Get_address(&values.processor_name_neighb2, &disp[15]);
  MPI_Get_address(&values.processor_name_neighb3, &disp[16]);
  MPI_Get_address(&values.processor_name_neighb4, &disp[17]);

  //Make relative
  disp[1]=disp[1]-disp[0];
  disp[2]=disp[2]-disp[0];
  disp[3]=disp[3]-disp[0];
  disp[4]=disp[4]-disp[0];
  disp[5]=disp[5]-disp[0];
  disp[6]=disp[6]-disp[0];
  disp[7]=disp[7]-disp[0];
  disp[8]=disp[8]-disp[0];
  disp[9]=disp[9]-disp[0];
  disp[10]=disp[10]-disp[0];
  disp[11]=disp[11]-disp[0];
  disp[12]=disp[12]-disp[0];
  disp[13]=disp[13]-disp[0];
  disp[14]=disp[14]-disp[0];
  disp[15]=disp[15]-disp[0];
  disp[16]=disp[16]-disp[0];
  disp[17]=disp[17]-disp[0];
  disp[0]=0;

  // Create MPI struct
  MPI_Type_create_struct(18, blocklen, disp, type, &Valuetype);
  MPI_Type_commit(&Valuetype);

  // default coordinate values of neighbours
  values.neighbour1_cord[0] = -2;
  values.neighbour1_cord[1] = -2;

  values.neighbour2_cord[0] = -2;
  values.neighbour2_cord[1] = -2;

  values.neighbour3_cord[0] = -2;
  values.neighbour3_cord[1] = -2;

  values.neighbour4_cord[0] = -2;
  values.neighbour4_cord[1] = -2;

  // default ranks
  values.rank[1] = -2;
  values.rank[2] = -2;
  values.rank[3] = -2;
  values.rank[4] = -2;

  //default height
  values.height[1] = -2;
  values.height[2] = -2;
  values.height[3] = -2;
  values.height[4] = -2;


  struct ifreq ifr;

  /*AF_INET - to define network interface IPv4*/
  /*Creating soket for it.*/
  fd = socket(AF_INET, SOCK_DGRAM, 0);

  /*AF_INET - to define IPv4 Address type.*/
  ifr.ifr_addr.sa_family = AF_INET;

  /*eth0 - define the ifr_name - port name
  where network attached.*/
  memcpy(ifr.ifr_name, "eth0", IFNAMSIZ - 1);

  /*Accessing network interface information by
  passing address using ioctl.*/
  ioctl(fd, SIOCGIFADDR, &ifr);
  /*closing fd*/
  close(fd);

  /*Extract IP Address*/
  strcpy(values.ipAddress, inet_ntoa(((struct sockaddr_in*)&ifr.ifr_addr)->sin_addr));


  // getting the size of the slave communicator
  MPI_Comm_size(comm, &size);

  // getting the rank of the slave communicator
  MPI_Comm_rank(comm, &my_rank);

	// Creating dimensions
	dims = (int*)malloc((ndims) * sizeof(int));

  // raises a error if the rows and columns arent equal to the size.
  if((nrows*ncols) != size) {
  printf("ERROR: nrows*ncols)=%d * %d = %d != %d\n", nrows, ncols, nrows*ncols,size);
  MPI_Finalize();
  return 0;
	}
  else{
    dims[0] = nrows; /* number of rows */
    dims[1] = ncols; /* number of columns */
  }


    /*************************************************************/
	/* create cartesian topology for processes */
	/*************************************************************/
	MPI_Dims_create(size, ndims, dims);

    /* create cartesian mapping */
	wrap_around[0] = 0;
	wrap_around[1] = 0; /* periodic shift is .false. */
	reorder = 0;
	ierr =0;
	ierr = MPI_Cart_create(comm, ndims, dims, wrap_around, reorder, &comm2D);
	if(ierr != 0) {
        printf("ERROR[%d] creating CART\n",ierr);
    }


    /* find my coordinates in the cartesian communicator group */
	MPI_Cart_coords(comm2D, my_rank, ndims, coord); // coordinated is returned into the coord array
	/* use my cartesian coordinates to find my rank in cartesian group*/
	MPI_Cart_rank(comm2D, coord, &my_cart_rank);
	/* get my neighbors; axis is coordinate dimension of shift */
	/* axis=0 ==> shift along the rows: P[my_row-1]: P[me] : P[my_row+1] */
	/* axis=1 ==> shift along the columns P[my_col-1]: P[me] : P[my_col+1] */

	MPI_Cart_shift( comm2D, SHIFT_ROW, DISP, &nbr_i_lo, &nbr_i_hi );
	MPI_Cart_shift( comm2D, SHIFT_COL, DISP, &nbr_j_lo, &nbr_j_hi );

    // printf("Cerated the cart \n");


    // We get the name of our processor using this
	char processor_name[50];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
  strncpy(values.processor_name, processor_name, name_len);
  values.processor_name[name_len] = '\0';


  // Initialise variables
  MPI_Status receive_status[4];

  MPI_Request ip_send_request[4];
  MPI_Status ip_send_status[4];

  MPI_Request processor_name_send_request[4];
  MPI_Status processor_name_send_status[4];


  MPI_Request ip_receive_request[4];
  MPI_Status ip_receive_status[4];


  MPI_Request processor_name_receive_request[4];
  MPI_Status processor_name_receive_status[4];


  // In this block current processor sends its ip address and processor name to its neighbours
  MPI_Isend(&values.ipAddress,15,MPI_CHAR,nbr_i_lo,IP_TAG,comm2D,&ip_send_request[0]);
  MPI_Isend(&values.ipAddress,15,MPI_CHAR,nbr_i_hi,IP_TAG,comm2D,&ip_send_request[1]);
  MPI_Isend(&values.ipAddress,15,MPI_CHAR,nbr_j_lo,IP_TAG,comm2D,&ip_send_request[2]);
  MPI_Isend(&values.ipAddress,15,MPI_CHAR,nbr_j_hi,IP_TAG,comm2D,&ip_send_request[3]);

  MPI_Isend(&values.processor_name,50,MPI_CHAR,nbr_i_lo,VM_NAME_TAG,comm2D,&processor_name_send_request[0]);
  MPI_Isend(&values.processor_name,50,MPI_CHAR,nbr_i_hi,VM_NAME_TAG,comm2D,&processor_name_send_request[1]);
  MPI_Isend(&values.processor_name,50,MPI_CHAR,nbr_j_lo,VM_NAME_TAG,comm2D,&processor_name_send_request[2]);
  MPI_Isend(&values.processor_name,50,MPI_CHAR,nbr_j_hi,VM_NAME_TAG,comm2D,&processor_name_send_request[3]);

  // In this block other neighbour processor recieves  ip address and processor name of its neighbour
  MPI_Irecv(&values.ipAddressNeighbour1, 15, MPI_CHAR, nbr_i_lo, IP_TAG, comm2D, &ip_receive_request[0]);
  MPI_Irecv(&values.ipAddressNeighbour2, 15, MPI_CHAR, nbr_i_hi, IP_TAG, comm2D, &ip_receive_request[1]);
  MPI_Irecv(&values.ipAddressNeighbour3, 15, MPI_CHAR, nbr_j_lo, IP_TAG, comm2D, &ip_receive_request[2]);
  MPI_Irecv(&values.ipAddressNeighbour4, 15, MPI_CHAR, nbr_j_hi, IP_TAG, comm2D, &ip_receive_request[3]);

  MPI_Irecv(&values.processor_name_neighb1, 50, MPI_CHAR, nbr_i_lo, VM_NAME_TAG, comm2D, &processor_name_receive_request[0]);
  MPI_Irecv(&values.processor_name_neighb2, 50, MPI_CHAR, nbr_i_hi, VM_NAME_TAG, comm2D, &processor_name_receive_request[1]);
  MPI_Irecv(&values.processor_name_neighb3, 50, MPI_CHAR, nbr_j_lo, VM_NAME_TAG, comm2D, &processor_name_receive_request[2]);
  MPI_Irecv(&values.processor_name_neighb4, 50, MPI_CHAR, nbr_j_hi, VM_NAME_TAG, comm2D, &processor_name_receive_request[3]);

  MPI_Waitall(4, ip_send_request, ip_send_status);
  MPI_Waitall(4, ip_receive_request, ip_receive_status);
  MPI_Waitall(4, processor_name_send_request, processor_name_send_status);
  MPI_Waitall(4, processor_name_receive_request, processor_name_receive_status);

    


  pthread_t tid;
  pthread_create(&tid, 0, ProcessFunc, NULL); // Create the thread

  while (100){

      MPI_Request termination_status;
      int termination_message_from_base;
      // It terminates the while loop when it gets message from the base station.
      MPI_Irecv(&termination_message_from_base,1,MPI_INT,0,TERMINATE_TAG,master_comm,&termination_status);
      if (termination_message_from_base == -1){
          thread_termination = -2;
          MPI_Comm_free(&comm2D);
          break;
      }

      // getting water depth detected by the node
      float water_depth = random_water_depth(my_rank,Threshold);


      // checking if array is full then we have to readjust the array to calculate simple moving average
      if (count_of_array< MAX_ARRAY_COUNT){
          simple_moving_average_val = simple_moving_average(count_of_array,water_column_array,water_depth);
          count_of_array++;

      }
      // checking if array is not full then add the new water depth detected by node and calculate the simple moving average
      else{
          simple_moving_average_val = simple_moving_average(count_of_array,water_column_array,water_depth);
          // printf("Working when array is full\n");
      }

      // checks if its limit threshold for tsunami occuring
      if (simple_moving_average_val > Threshold){
          int message = my_rank;
          int neighbours[4] = {nbr_i_lo,nbr_i_hi,nbr_j_lo,nbr_j_hi};

          for(int i =0 ; i <4; i++){
              if (neighbours[i] >= 0){
                  MPI_Send(&message,1,MPI_INT,neighbours[i],SEND_REQ,comm2D);
              }
          }



          // In this block processors receives simple movinge average value from its neighbours
          float recvValL = -1, recvValR = -1, recvValT = -1, recvValB = -1;
          if (nbr_i_lo >= 0){
              MPI_Recv(&recvValL,1,MPI_FLOAT,nbr_i_lo,RECV_REQ,comm2D,&receive_status[0]);
          }
          if (nbr_i_hi >= 0){
              MPI_Recv(&recvValR,1,MPI_FLOAT,nbr_i_hi,RECV_REQ,comm2D,&receive_status[1]);
          }
          if (nbr_j_lo >= 0){
              MPI_Recv(&recvValT,1,MPI_FLOAT,nbr_j_lo,RECV_REQ,comm2D,&receive_status[2]);
          }
          if (nbr_j_hi >= 0){
              MPI_Recv(&recvValB,1,MPI_FLOAT,nbr_j_hi,RECV_REQ,comm2D,&receive_status[3]);

          }

          int count_neighbours = 0;
          if ( nbr_i_lo >=0 && abs(simple_moving_average_val - recvValT) <= TOLERENCERANGE){
                count_neighbours++;
          }

          if (nbr_i_hi >=0  && abs(simple_moving_average_val - recvValB) <= TOLERENCERANGE){
              count_neighbours++;
          }

          if (nbr_j_lo && abs(simple_moving_average_val - recvValL) <= TOLERENCERANGE){
              count_neighbours++;
          }

          if (nbr_j_hi && abs(simple_moving_average_val - recvValR)<= TOLERENCERANGE){
              count_neighbours++;
              
          }

          // if we have atleast two or more neighbours whos  simple_moving_average_val when compared with the simple_moving_average_val of the tsnumai sensor who
          // requested neighbour simple_moving_average_val is less or equal to Threshold then report to base station
          if (count_neighbours >= 2){
            values.timeStamp = time(NULL);
            values.rank[0] = my_rank;
            values.my_cord[0] = coord[0];
            values.my_cord[1] = coord[1];
            values.height[0] =simple_moving_average_val;


            // check if top neighbour exists in the grid for the current sensor  
            if (nbr_i_lo >=0 ){
                MPI_Cart_coords(comm2D, nbr_i_lo, ndims, values.neighbour1_cord); // coordinated is returned into the coord array
                values.rank[1] = nbr_i_lo;
                values.height[1] = recvValT;
              

            }
            // check if bottom neighbour exists in the grid for the current sensor  
            if (nbr_i_hi >=0 ){
                MPI_Cart_coords(comm2D, nbr_i_hi, ndims, values.neighbour2_cord); // coordinated is returned into the coord array
                values.rank[2] = nbr_i_hi;
                values.height[2] = recvValB;
            
            }
            // check if left neighbour exists in the grid for the current sensor 
            if (nbr_j_lo >=0 ){
                MPI_Cart_coords(comm2D, nbr_j_lo, ndims, values.neighbour3_cord); // coordinated is returned into the coord array
                values.rank[3] = nbr_j_lo;
                values.height[3]= recvValL;
                
            

            }
            // check if right neighbour exists in the grid for the current sensor 
            if (nbr_j_hi >=0 ){
                MPI_Cart_coords(comm2D, nbr_j_hi, ndims, values.neighbour4_cord); // coordinated is returned into the coord array
                values.rank[4] = nbr_j_hi;
                values.height[4] = recvValR;

            }

            // send the data to the base station so it can write the report
            MPI_Send(&values, 1, Valuetype, 0, SENDING_INFO_TAG, master_comm);
        }
    }
    sleep(5);

    }

  pthread_join(tid, NULL); // Wait for the thread to complete.
  return 0;

}

// This function generates random float values which represents the  water column height of tsnuami sensor according to Threshold
float random_water_depth(int my_rank, float Threshold){
    srand((unsigned int)time(NULL) * (my_rank));
    float scale = rand() / (float) RAND_MAX; /* [0, 1.0] */
    return (Threshold-1000) + scale * ((Threshold+1000) - (Threshold-1000));      /* [min, max] */
}


// This function calculates the simple mobing average of tsnuami sensors
float simple_moving_average(int count,float arr[],float new_water_depth){

    int i;
    float sum = 0.0;
    float avg;
    if (count == MAX_ARRAY_COUNT){
        for(i=0; i < MAX_ARRAY_COUNT-1;i++){
            sum +=  arr[i+1];
            arr[i] = arr[i+1];
        }
        arr[MAX_ARRAY_COUNT-1] = new_water_depth;
        sum += new_water_depth;
        avg = sum/MAX_ARRAY_COUNT;
    }
    else{
        arr[count] = new_water_depth;
        for(i=0; i < count+1;i++){
            sum +=  arr[i];
        }
        avg = sum/(count+1);
    }
    return avg;
}

/* Function to write into the file */
void WriteToFile(char *pFilename, char *txt, bool newline)
{
    
	FILE *pFile = fopen(pFilename, "a");
	fprintf(pFile, "%s", txt);

    // If newline is true, output a newline at the end
	if (newline){
		fprintf(pFile, "\n");
	}

	fclose(pFile);
}

/* Function to write numbers into the file. */
void WriteNumToFile(char *pFilename, int num, bool newline)
{
	FILE *pFile = fopen(pFilename, "a");
	fprintf(pFile, "%d", num);

	// If newline is true, output a newline at the end
	if (newline){
		fprintf(pFile, "\n");
	}

	fclose(pFile);
}

/* Function to write floats into the file. */
void WriteFloatToFile(char *pFilename, float num, bool small, bool newline)
{
	FILE *pFile = fopen(pFilename, "a");

    if (small){
        fprintf(pFile, "%.2f", num);
    }
	else{
        fprintf(pFile, "%.6f", num);
    }

	// If newline is true, output a newline at the end
	if (newline){
		fprintf(pFile, "\n");
	}

	fclose(pFile);
}





//just print alert message to console using single printf
void alert(const char* msg)
{
  printf("%s\n", msg);
}


// defines used in code for messages ques
#define MSG_BUF_MAX_SIZE 256
#define MSG_DEFAULT_PRIORITY 1

// opens messages que for writing (from one proc to another)
// returns message handler
int message_que_open_for_write(const char* que_name)
{
  return mq_open(que_name, O_WRONLY);
}

// opens messages que for reading (from other to one)
// returns message handler
int message_que_open_for_read(const char* que_name)
{
  // fill que attributes
  struct mq_attr attr;
  attr.mq_flags = 0;
  attr.mq_maxmsg = 10;
  attr.mq_msgsize = MSG_BUF_MAX_SIZE;
  attr.mq_curmsgs = 0;

  return mq_open(que_name, O_CREAT | O_RDONLY  | O_NONBLOCK, 0644, &attr);
}

// get message from que and place it to passed pointer
// returns 1 if new message received
// returns 0 if no messages received
int get_message_from_que(int mh, char* msg)
{
  // priority of message
  unsigned int p;
  // read bytes
  return mq_receive(mh, msg, MSG_BUF_MAX_SIZE, &p);
}

// put message to que from passed pointer
void put_message_to_que(int mh, char* msg)
{
  // send bytes
  mq_send(mh, msg, MSG_BUF_MAX_SIZE, MSG_DEFAULT_PRIORITY);
}

// closes messages que
void message_que_close(int mh)
{
  // close message handler
  mq_close (mh);
}

// unlink messages que
void message_que_unlink(const char* name)
{
  // unlink message handler
  mq_unlink(name);
}

// fill wave record with neccesary data
void fill_wave_record(wave_record_t* record, int x, int y, int h)
{
  // sanity check - passed pointer should be not null
  if(!record)
    return;

  // set time to record
  time_t t = time(NULL);
  struct tm * lt = localtime(&t);

  record->year = lt->tm_year + 1900;
  record->month = lt->tm_mon + 1;
  record->day = lt->tm_mday;
  record->hour = lt->tm_hour;
  record->minute = lt->tm_min;
  record->second = lt->tm_sec;

  // set passed parameters
  record->height = h;
  record->x = x;
  record->y =y;
}

// represent table record as signle string
void convert_record_to_string(wave_record_t* record, char* str)
{
  // sanity check - passed pointer should be not null
  if ((!record) || (!str))
    return;
  
  // print values for each field into the passed char pointer
  sprintf(str, "%04d | %02d | %02d | %02d | %02d | %02d | %02d, %02d | %05.03f",
    record->year, record->month, record->day,
    record->hour, record->minute, record->second,
    record->x, record->y, (record->height/1000.0));
  
}

// makes printf of record values to screen
void print_record(wave_record_t* record)
{
  // sanity check - passed pointer should be not null
  if (!record)
    return;

  // variable for printing head first time
  static int head_output = 1;
  if (head_output)
  {
    // printf("YYYY | MM | DD | HH | MM | SS | XX, YY | HEIGHT\n");
    // printf("----------------------------------------------\n");
    head_output = 0;
  }
  
  // create char* array and convert record to printable string
  char str[MSG_BUF_MAX_SIZE] = {0};
  convert_record_to_string(record, str);
  
  // print out record
  printf("%s\n", str);
}

// copy fileds from one record to another
void copy_wave_record(wave_record_t* from, wave_record_t* to)
{
  // sanity check - passed pointer should be not null
  if ((!from) || (!to))
    return;

  // copy values for each field
  memcpy(to, from, sizeof(wave_record_t));
}

// initilize waves table
void wave_table_init(wave_table_t* table, int size)
{
  // sanity check - passed pointer should be not null
  if(!table)
    return;

  // lock mutex to lock access to table in other functions
  pthread_mutex_lock(&table->mutex);

  // allocate memory for storing wave records
  table->records = (wave_record_t*)malloc(sizeof(wave_record_t)*size);

  // reset all table parameters to zero
  table->pos_read = 0;
  table->pos_write = 0;
  table->count = 0;

  // setup table size
  table->size = size;

  // unlock mutex to allow access to table in other functions
  pthread_mutex_unlock(&table->mutex);
}

// function that adds wave record to table (place to fifo)
void wave_table_add(wave_table_t* table, wave_record_t* record)
{
  // sanity check - passed pointer should be not null
  if ((!record) || (!table))
    return;

  // lock mutex to lock access to table in other functions
  pthread_mutex_lock(&table->mutex);

  print_record(record);

  // write record to write position
  copy_wave_record(record, &table->records[table->pos_write]);

  // increment table records count
  table->count++;

  // check records count not greate than table size
  if (table->count > table->size)
    table->count = table->size;

  // shift write position to next (this code does not overlaps over table->size
  // and works in range [0 ... table->size - 1])
  table->pos_write++;
  table->pos_write = table->pos_write % table->size;

  // in case if 'pos_read' will be overwrited by 'pos write'
  // shift pos_read to next after pos_write
  // this works when new data placed before old was read
  if (table->pos_read == table->pos_write)
  {
    table->pos_read++;
    table->pos_read = table->pos_read % table->size;
  }

  // unlock mutex to allow access to table in other functions
  pthread_mutex_unlock(&table->mutex);
}

// function that get single wave from table 
// if get was success function will return 1, else 0
int wave_table_get(wave_table_t* table, wave_record_t* record)
{
  // ret code declaration, by default 0 returned
  int ret_code = 0;

  // sanity check - passed pointer should be not null
  if ((!record) || (!table))
    return ret_code;

  // lock mutex to lock access to table in other functions
  pthread_mutex_lock(&table->mutex);

  // check we have records to put to pointer
  if (table->count > 0)
  {
    //put first record that was not previously read
    copy_wave_record(&table->records[table->pos_read], record);

    // decrement table records count
    table->count--;

    // shift read position to next (this code do not overlapss over table->size
    // and works in range [0 ... table->size - 1])
    table->pos_read++;
    table->pos_read = table->pos_read % table->size;

    // setup ret code
    ret_code = 1;
  }

  // unlock mutex to allow access to table in other functions
  pthread_mutex_unlock(&table->mutex);

  // return 1 or 0
  return ret_code;
}

// used to check passed string table record with last geted
// return 1 if compared was ok or 0 if compared was not ok
int wave_table_check_with_str(wave_table_t* table, char* str)
{
  // ret code declaration, by default 0 returned
  int ret_code = 0;

  // sanity check - passed pointer should be not null
  if ((!str) || (!table))
    return ret_code;

  // lock mutex to lock access to table in other functions
  pthread_mutex_lock(&table->mutex);

  // calculate last read position 
  // shift read position to previous one (this code do not overlapss over table->size
  // and works in range [0 ... table->size - 1])
  int last_read_pos = table->pos_read;
  if (last_read_pos == 0)
    last_read_pos = table->size - 1;
  else
    last_read_pos--;
  
  // get record on this position and make string line with it is data
  char record_str[MSG_BUF_MAX_SIZE] = {0};
  convert_record_to_string(&table->records[last_read_pos], record_str);
  
  // compare strings - equal or not?
  if (strcmp(record_str, str) == 0)
    ret_code = 1; // if yes - we will return 1
  else
    ret_code = 0; // otherwice - 0
  
  // unlock mutex to allow access to table in other functions
  pthread_mutex_unlock(&table->mutex);

  // return 1 or 0
  return ret_code;
}

// Generates waves for specified grid layout (X x Y) with random height
// Fill generated data to passed table
//
void periodical_generation(wave_table_t* table)
{
  int x;
  int y;

  // loop on x
  for(x = 0; x < nrows; x++)
  {
    // loop on y
    for (y = 0; y < ncols; y++)
    {
      // generate random height
      int height = Threshold+1 + rand() / (RAND_MAX / (Threshold+1000 - Threshold+1 + 1) + 1);
      height = height*1000;

        wave_record_t record;
        fill_wave_record(&record, x, y, height);

        // and add that record to table
        wave_table_add(table, &record);
      }
    }
  }

// Thread foo that implements generation of waves and placing them to table
// Reads messages and depending on message provide data or exit from loop
void* wave_table_thread(void* param)
{
  // open message que, store message handler
  int mr = message_que_open_for_read(WAVE_MESSAGE_QUE_TO_THREAD);
  if (mr < 0)
  {
    // this is error condition
    printf ("erro creating que for reader in thread\n");
    exit(1);
  }

  // loop until que will be opened for communicatiobn with main loop
  int mw;
  do
  {
    mw = message_que_open_for_write(WAVE_MESSAGE_QUE_FROM_THREAD);
  } while (mw == -1);
  
  // in loop flag - used to break loop
  int in_loop = 1;

  // time of start of the loop
  time_t t_generation = time(NULL);
  time_t t_transmision = time(NULL);

  // go to infinite loop
  while (in_loop)
  {
    // used for receiving messages
    char msg[MSG_BUF_MAX_SIZE] = {0};

    // read message and if message received - process it
    if (get_message_from_que(mr, msg) > 0)
    {
      // check message - if it is 'exit' - reset in_loop to zero
      if (strcmp(msg, "exit") == 0)
      {
        in_loop = 0;
        break;
      }
    }

    // check time spended in loop
    int timeSpent;
    time_t t2 = time(NULL);

    //if spent more than interval of generation
    timeSpent = (int)difftime(t2, t_generation);
    if (timeSpent >= WAVE_GENERATION_INTERVAL)
    {
      // generate new records to table
      periodical_generation(&wave_table);
      // and update time for future interval calculation
      t_generation = time(NULL);
    }
    
    //if spent more than interval of transmission
    timeSpent = (int)difftime(t2, t_transmision);
    if (timeSpent >= WAVE_TRANSMISSION_INTERVAL)
    {
      // transmit record from table
      wave_record_t record;
      
      // get record from table, if table not empty it will return value = 1
      if (wave_table_get(&wave_table, &record) > 0)
      {
        // convert record data to char* (string) array 
        char msg[MSG_BUF_MAX_SIZE] = {0};
        convert_record_to_string(&record, msg);
        
        // and send as message into the que
        put_message_to_que(mw, msg);
      }
        
      // and update time for future interval calculation
      t_transmision = time(NULL);
    }
    
    // periodically sleep in thread to get time for other can do theirs job
    usleep(1000);
  }
  // close message que
  message_que_close(mr);

  // exit from thread
  pthread_exit(0);
}


// function that should be called to start table generation
void altimeter_simulation()
{
  // turn on ques in linux system - you need to be as root to do that
  system("mkdir -p /dev/mqueue");
  system("mount -t mqueue none /dev/mqueue");
 
  // initialize random timer
  srand(time(NULL));

  // intialize waves table (fifo)
  wave_table_init(&wave_table, WAVE_TABLE_SIZE);
  
  // unlink previously created ques if exist
  message_que_unlink(WAVE_MESSAGE_QUE_TO_THREAD);
  message_que_unlink(WAVE_MESSAGE_QUE_FROM_THREAD);
  
  // create message que for getting lines from wave table 
  int mr = message_que_open_for_read(WAVE_MESSAGE_QUE_FROM_THREAD);
  if (mr < 0)
  {
    // this is error condition
    printf ("erro creating que for reader in main loop\n");
    exit(0);
  }
  
  // start thread for wave table records genration
  pthread_t tid;
  pthread_create(&tid, NULL, wave_table_thread, NULL);

  // loop until que will be open for communicatiobn with thread
  int mw;
  do
  {
    mw = message_que_open_for_write(WAVE_MESSAGE_QUE_TO_THREAD);
  } while (mw == -1);

  // then go to the loop where we will wait N minutes then send 'exit' to thread and quit the append
  time_t t_start = time(NULL);
  time_t t_receive = time(NULL);
  
  int in_loop = 1;
  while (in_loop)
  {
    // get current time and calculate time spen from moment t1
    int timeSpent;
    time_t t2 = time(NULL);
    
    // if time spent reached total time allowed to app be in work state
    timeSpent = (int)difftime(t2, t_start);
    if (timeSpent >= WAVE_APPLICATION_WORK_TIME)
    {
      // reset in loop flag
      in_loop = 0;

      // send message to thread to exit
      char msg[MSG_BUF_MAX_SIZE] = {0};
      sprintf(msg, "%s", "exit");
      put_message_to_que(mw, msg);

      // quit from loop
      break;
    }
    
    // if time spent greater than transmission interval - read the data 
    timeSpent = (int)difftime(t2, t_receive);
    if (timeSpent >= WAVE_TRANSMISSION_INTERVAL)
    {
      char msg[MSG_BUF_MAX_SIZE] = {0};
      
      // read message (it should present) and if message received - process it
      if (get_message_from_que(mr, msg) > 0)
      {
        // if message present - compare with original in table
        int comparation_result = wave_table_check_with_str(&wave_table, msg);
        
        // and if comparation result != 1
        if (comparation_result != 1)
        {
          // make alert with notice
          alert("fault detection - record wrong");
        }
      }
      else
      {
        // no message present - make alert
        alert ("no data received from node");
      }

      // and update time for future interval calculation
      t_receive = time(NULL);
    }
    
    // periodically sleep in thread to get time for other can do theirs job
    usleep(1000);
  }

  // wait till wave thread will done
  pthread_join(tid,NULL);

  // close message que
  message_que_close(mw);
  message_que_unlink(WAVE_MESSAGE_QUE_TO_THREAD);
  message_que_unlink(WAVE_MESSAGE_QUE_FROM_THREAD);

}



