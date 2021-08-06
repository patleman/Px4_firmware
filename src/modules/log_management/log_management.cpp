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
#include "log_management.hpp"

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
extern "C" __EXPORT int log_management_main(int argc, char *argv[]);

/**
 * Mainloop of daemon.
 */
int log_management_thread_main(int argc, char *argv[]);

/**
 * Print the correct usage.
 */
static void usage(const char *reason);


/* Variables */
static bool thread_should_exit = false;		/**< Daemon exit flag */
static bool thread_running = false;		/**< Daemon status flag */
static int deamon_task;				/**< Handle of deamon task / thread */



/* Main Thread */
int log_management_thread_main(int argc, char *argv[])
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

bool ics_updated;
	
	int ics_sub = orb_subscribe(ORB_ID(initial_check_status));

		struct initial_check_status_s ics_fetch;
	memset(&ics_fetch, 0, sizeof(ics_fetch));

	while(1){
		orb_check(ics_sub, &ics_updated);
		if(ics_updated){
		
	        orb_copy(ORB_ID(initial_check_status), ics_sub, &ics_fetch);
		//	aux=ics_fetch.identifier_status;
        
             break;
		}
		usleep(500000);
	}
	
    struct initial_check_status_s ics;
	memset(&ics, 0, sizeof(ics));
		
	orb_advert_t ics_pub = orb_advertise(ORB_ID(initial_check_status), &ics);  


//function 4 : recentPA.txt for flight log management.
//check for recentPA.txt:::: this is for bundling purpose, do we need to start bundling or not
//
int recentPA_presence=0;

char PA_ID[100];

FILE *fptr_recent;



FILE *fptr_pa;

char pa_address[50]="/fs/microsd/log/permission_artifact_breach.xml";

fptr_pa=fopen(pa_address,"r");

if(fptr_pa==NULL){
	
	ics.pa_present=0;
}else{
	fclose(fptr_pa);


	ics.pa_present=1;
}

fptr_recent=fopen("/fs/microsd/recentPA.txt","r");

	if(fptr_recent!=NULL){
        
		//if file is present then there might be the possiblity upcoming PA to be the old one.
		fclose(fptr_recent);
		int action=check_recentPA();

    if(action!=3){
			char log_name[100];
			int needTosign=check_for_sign(log_name);

			if(needTosign){
			
				signing_support_log(log_name);// later 

			}

			int needTobundle=check_for_bundle();
			if(needTobundle){
				Bundling_begins();
				char yes[2]="1";
			    update_recentPA(0,yes);// setting fetch required to 1
			    recentPA_presence++;
			    ics.recent_pa_present=1;
			}
	    }
		//checks to be done :
		//1)if the current time lies inside the time period of start time and end time
		// this will initiate the flight logs bundling and then the drone will wait for fetched.txt file
		//updated through a parameter in recentPA.txt

		if(action==0 || action==1){
			if(!needTobundle){
				// start bundling and update recentPA.txt with fetchrequired=0;
			    //Bundling_begins();//how to bundle part/modifying and updating recentPA.txt
				char yes[2]="1";
				update_recentPA(0,yes);// setting fetch required to 1
				recentPA_presence++;
				ics.recent_pa_present=1;
			}
		}else if(action==3){
			// file is not valid 
			ics.recent_pa_present=0;
			remove("/fs/microsd/recentPA.txt");

		}else{//no need to start bundling 2

			//make sure upcoming PA has same pa_id. as in recentPA.txt
			ics.recent_pa_present=1;
			recentPA_presence++;

		}

		//check inside recentPA.txt if any previous log file is needed to be signed


		// coming here means recentPA is present, and there might be a chance pa is also present at this moment
		// case1 pa is present 
		//case 2 pa is not present

	}else{
		//file not present
		// upcoming permission artefact is new, a fresh recentPA.txt file would be made.
	   
	   
	    ics.recent_pa_present=0;
	   
	   


	}



/// to check for permission artifact (is it present here or not) and its validation
/// then verifying for geo coordinates and time period
///
//Now first read recentPA.txt to know if  fetch of fetched.txt is required or not

int check_status_fetch=0;
int status_fetch=0;
if(recentPA_presence)
{
	// to know if we are supposed to look for fetched.txt
	status_fetch=read_for_fetch();

	if(status_fetch!=1){
		//fetch is not required (no need to check for fetched.txt file_)
	}else{
		
		// fetch is required (need to check for fetched.txt file)
		check_status_fetch=check_fetch();
		if(check_status_fetch==0 || check_status_fetch==1 ){
			ics.fetch_required=1;
		}else{
			ics.fetch_required=0;
			char no[3]="0";
			update_recentPA(0,no);
		}
		//if check_status==0 : no file present, return false
		//               ==1  : file is present but not valid return false(attempt to hack)
		//               ==2  : file is present and valid, delete recentPA.txt and set recentPA_presence=0
	}
}
    ics.registration_status=ics_fetch.registration_status;
	ics.identifier_status=ics_fetch.identifier_status;
	
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
int log_management_main(int argc, char *argv[])
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
		deamon_task = px4_task_spawn_cmd("log_management",
						 SCHED_DEFAULT,
						 SCHED_PRIORITY_MAX - 20,
						 20048,
						 log_management_thread_main,
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



