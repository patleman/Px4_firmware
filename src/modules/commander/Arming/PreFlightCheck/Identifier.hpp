
// This file has the function to check for truthness of hardware being used
/* Files in play :
1)HardwareChange.txt
2)

*/
#include<stdio.h>
#include<stdarg.h>
using namespace std;
#include<string.h>
#include<inttypes.h>
#define DRONE_ID "ABCD123PA456";
#include <uORB/topics/sensor_gps.h>


//function for maintaining hardware integrity
int RPAS_identifier(){
    // taking gps id
	int vehi_gps=orb_subscribe(ORB_ID(sensor_gps));
        struct sensor_gps_s raw;
        orb_copy(ORB_ID(sensor_gps),vehi_gps,&raw);
        int deviceID=raw.lat;
	printf("Device ID of GPS is: %d ",deviceID);

	return 1;

}
