/****************************************************************************
 *
 *   Copyright (c) 2013-2015 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file main.cpp
 *
 * Example implementation of a rover steering controller.
 *
 * @author Lorenz Meier <lorenz@px4.io>
 */

#include <px4_platform_common/px4_config.h>
#include <px4_platform_common/tasks.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>
#include <poll.h>
#include <uORB/topics/indi_prc.h>
#include <uORB/Subscription.hpp>
#include <uORB/SubscriptionInterval.hpp>


#include <uORB/topics/parameter_update.h>
#include <parameters/param.h>
#include <lib/ecl/geo/geo.h>
#include <perf/perf_counter.h>
#include <systemlib/err.h>
#include <matrix/math.hpp>
//#include "json_support.h"

/* process-specific header files */
#include "registration_check.hpp"

using namespace time_literals;

/* Prototypes */

/**
 * Daemon management function.
 *
 * This function allows to start / stop the background task (daemon).
 * The purpose of it is to be able to start the controller on the
 * command line, query its status and stop it, without giving up
 * the command line to one particular process or the need for bg/fg
 * ^Z support by the shell.
 */
extern "C" __EXPORT int registration_check_main(int argc, char *argv[]);

/**
 * Mainloop of daemon.
 */
int registration_check_thread_main(int argc, char *argv[]);

/**
 * Print the correct usage.
 */
static void usage(const char *reason);


/* Variables */
static bool thread_should_exit = false;		/**< Daemon exit flag */
static bool thread_running = false;		/**< Daemon status flag */
static int deamon_task;				/**< Handle of deamon task / thread */



/* Main Thread */
int registration_check_thread_main(int argc, char *argv[])
{



	/*
	 * PX4 uses a publish/subscribe design pattern to enable
	 * multi-threaded communication.
	 *
	 * The most elegant aspect of this is that controllers and
	 * other processes can either 'react' to new data, or run
	 * at their own pace.
	 *
	 * PX4 developer guide:
	 * https://pixhawk.ethz.ch/px4/dev/shared_object_communication
	 *
	 * Wikipedia description:
	 * http://en.wikipedia.org/wiki/Publish–subscribe_pattern
	 *
	 */
    
	FILE *fptr_registration;
	fptr_registration=fopen("/fs/microsd/UUID.txt","r");
	// publish the identtiy check status to the uorb topic  

    bool ics_updated;
	
	int ics_sub = orb_subscribe(ORB_ID(initial_check_status));

		struct initial_check_status_s ics_fetch;
	memset(&ics_fetch, 0, sizeof(ics_fetch));
    int aux;
	while(1){
		orb_check(ics_sub, &ics_updated);
		if(ics_updated){
		
	        orb_copy(ORB_ID(initial_check_status), ics_sub, &ics_fetch);
			aux=ics_fetch.identifier_status;
        
             break;
		}
		usleep(500000);
	}
	
    struct initial_check_status_s ics;
	memset(&ics, 0, sizeof(ics));
		
	orb_advert_t ics_pub = orb_advertise(ORB_ID(initial_check_status), &ics);  


	if(fptr_registration!=NULL)
	{	//file is present
		fclose(fptr_registration);
		char DroneID_container_1[30];
		char TAG_DRONE_ID_1[10]="DroneID";
		char DRONE_ID_FILE_1[100]="/fs/microsd/UUID.txt";
		int check_validity_DroneID_1=file_read(DRONE_ID_FILE_1,TAG_DRONE_ID_1,DroneID_container_1,1);
		// if non tampered -->then continue
		FILE *fptr;
		fptr=fopen("/fs/microsd/debug_uuid.txt","w");
		fprintf(fptr,"status  :::: %d\n",check_validity_DroneID_1);
		fprintf(fptr,"status1  :::: %d\n",ics_fetch.identifier_status);
		fclose(fptr);
		if (check_validity_DroneID_1==0){
			//UUID file is not valid, remove the invalid file
			remove("/fs/microsd/UUID.txt");
			ics.registration_status=3;
			ics.identifier_status=ics_fetch.identifier_status;
	
		}else{
			//file is valid, check for strings
			if(strcmp(DroneID_container_1,"ABBDDJEDNDJK")==0){
			//Drone is registered
				printf("\n\nDrone is registered\n\n");
				ics.registration_status=1;
				ics.identifier_status=aux;
				
			}else{
				// UUID.txt file is for different drone
				// not UUID.txt of this drone
				ics.registration_status=2;
				ics.identifier_status=ics_fetch.identifier_status;
				remove("/fs/microsd/UUID.txt");
			}

		}



	}
	else{
	// registration has not been done, also recommend to connect to management client and
	// internet
	// UUID.txt file not present
           ics.registration_status=0;
		   ics.identifier_status=ics_fetch.identifier_status;
		
	}
	ics.timestamp=hrt_absolute_time();
	orb_publish(ORB_ID(initial_check_status), ics_pub, &ics);



	return 0;
}

/* Startup Functions */

static void
usage(const char *reason)
{
	if (reason) {
		fprintf(stderr, "%s\n", reason);
	}

	fprintf(stderr, "usage: rover_steering_control {start|stop|status}\n\n");
}

/**
 * The daemon app only briefly exists to start
 * the background job. The stack size assigned in the
 * Makefile does only apply to this management task.
 *
 * The actual stack size should be set in the call
 * to px4_task_spawn_cmd().
 */
int registration_check_main(int argc, char *argv[])
{
	if (argc < 2) {
		usage("missing command");
		return 1;
	}

	if (!strcmp(argv[1], "start")) {

		if (thread_running) {
			warnx("running");
			/* this is not an error */
			return 0;
		}
        
		thread_should_exit = false;
		deamon_task = px4_task_spawn_cmd("registration_check",
						 SCHED_DEFAULT,
						 SCHED_PRIORITY_MAX - 20,
						 35048,
						 registration_check_thread_main,
						 (argv) ? (char *const *)&argv[2] : (char *const *)nullptr);
		thread_running = true;
		return 0;
	}

	if (!strcmp(argv[1], "stop")) {
		thread_should_exit = true;
		return 0;
	}

	if (!strcmp(argv[1], "status")) {
		if (thread_running) {
			warnx("running");

		} else {
			warnx("not started");
		}

		return 0;
	}

	usage("unrecognized command");
	return 1;
}



