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
#include "pa_valid.hpp"

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
extern "C" __EXPORT int pa_valid_main(int argc, char *argv[]);

/**
 * Mainloop of daemon.
 */
int pa_valid_thread_main(int argc, char *argv[]);

/**
 * Print the correct usage.
 */
static void usage(const char *reason);


/* Variables */
static bool thread_should_exit = false;		/**< Daemon exit flag */
static bool thread_running = false;		/**< Daemon status flag */
static int deamon_task;				/**< Handle of deamon task / thread */



/* Main Thread */
int pa_valid_thread_main(int argc, char *argv[])
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

    char canon_content[2000];
   if(ics_fetch.pa_present==1){

		//check if file is present or not
		int valid_status=pa_validation(canon_content);

		if(ics_fetch.recent_pa_present==1 && valid_status==1)
		{
			fetching_publish_padata();

		}else if(ics_fetch.recent_pa_present==0 && valid_status==1)
		{
			//case of new pa
			data_fetch_delete(canon_content);

		}else{
			//this is the case when pa is not valid
		}

   }else
   {
	   //pa not present
	   if(ics_fetch.recent_pa_present==1){
		   //this is done to check if incase bundling is needed to be done
		   fetching_publish_padata();    
	   }else{
		   // this is the case when there is no pa nor recentPA
	   }
   }
	//check 


    if(ics_fetch.pa_present){
		// check DroneID
		int pa_same_status;
		
		if(ics_fetch.recent_pa_present )
		{
				//recentPA present
			char fname[40]="/fs/microsd/log/recentPA.txt";
			char tag_paid[25]="permissionArtifactId";
			char paid_con[80];
			fetch_tag(NULL,tag_paid,paid_con,fname);
			pa_same_status= DroneIDverification(paid_con);


		}else{
			pa_same_status= DroneIDverification(NULL);

		}

		if(pa_same_status==1){
			ics.same_pa=1;
			ics.drone_id=1;
		}else if(pa_same_status==2){
			ics.same_pa=0;
		}
		else{
			// bad drone id
			ics.drone_id=0;
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
int pa_valid_main(int argc, char *argv[])
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
		deamon_task = px4_task_spawn_cmd("pa_valid",
						 SCHED_DEFAULT,
						 SCHED_PRIORITY_MAX - 20,
						 35048,
						 pa_valid_thread_main,
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



