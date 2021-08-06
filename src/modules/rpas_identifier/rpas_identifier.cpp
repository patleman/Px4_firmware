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
#include "rpas_identifier.hpp"

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
extern "C" __EXPORT int rpas_identifier_main(int argc, char *argv[]);

/**
 * Mainloop of daemon.
 */
int rpas_identifier_thread_main(int argc, char *argv[]);

/**
 * Print the correct usage.
 */
static void usage(const char *reason);


/* Variables */
static bool thread_should_exit = false;		/**< Daemon exit flag */
static bool thread_running = false;		/**< Daemon status flag */
static int deamon_task;				/**< Handle of deamon task / thread */



/* Main Thread */
int rpas_identifier_thread_main(int argc, char *argv[])
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
    int identity_check;
	FILE *fptr;
	fptr=fopen("/fs/microsd/hardware_status.txt","w");
	identity_check= RPAS_identifier(fptr);//file_read(fname,tag,res,0);
	fclose(fptr);
    

	// publish the identtiy check status to the uorb topic  
	struct initial_check_status_s ics;
	memset(&ics, 0, sizeof(ics));
		
	orb_advert_t ics_pub = orb_advertise(ORB_ID(initial_check_status), &ics);  
	switch (identity_check)
	{
		case 0:ics.identifier_status=0;
			break;
		case 1:ics.identifier_status=1;
			break;
		case 2:ics.identifier_status=2;
			break;
			
	}
	
	  ics.timestamp=hrt_absolute_time();
	orb_publish(ORB_ID(initial_check_status), ics_pub, &ics);


	FILE *droneid;
	droneid=fopen("/fs/microsd/DroneID.txt","r");
	if(droneid==NULL){
		// the file is not there in the system, it has to be created
		// currently without signing, later has to be with signing support(Key_pair_C)

		DroneIDcreation( );//DroneIDcreation

	}else{
		fclose(droneid);
		//file is there in the system check if it was witten by RPAS or not.
		// by checking for the signature
		char DroneID_container[30];
		char TAG_DRONE_ID[10]="DroneID";
		char DRONE_ID_FILE[100]="/fs/microsd/DroneID.txt";
		int check_validity_DroneID=file_read(DRONE_ID_FILE,TAG_DRONE_ID,DroneID_container,0);
		// if non tampered -->then continue
		if (check_validity_DroneID==0){
			remove("/fs/microsd/DroneID.txt");
			DroneIDcreation( );
		}
	// if tampered --> then delete it and create a genuine one. with a valid signature
	}


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
int rpas_identifier_main(int argc, char *argv[])
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
		deamon_task = px4_task_spawn_cmd("rpas_identifier",
						 SCHED_DEFAULT,
						 SCHED_PRIORITY_MAX - 20,
						 35048,
						 rpas_identifier_thread_main,
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



