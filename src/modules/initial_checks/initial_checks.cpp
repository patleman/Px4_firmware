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
#include "initial_checks.hpp"

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
extern "C" __EXPORT int initial_checks_main(int argc, char *argv[]);

/**
 * Mainloop of daemon.
 */
int initial_checks_thread_main(int argc, char *argv[]);

/**
 * Print the correct usage.
 */
static void usage(const char *reason);


/* Variables */
static bool thread_should_exit = false;		/**< Daemon exit flag */
static bool thread_running = false;		/**< Daemon status flag */
static int deamon_task;				/**< Handle of deamon task / thread */



/* Main Thread */
int initial_checks_thread_main(int argc, char *argv[])
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

	struct vehicle_global_position_s vgp;
	memset(&vgp, 0, sizeof(vgp));


	int vgp_sub = orb_subscribe(ORB_ID(vehicle_global_position));


        orb_set_interval(vgp_sub,1000);

	hrt_abstime then=hrt_absolute_time();

	struct pollfd fds[1] {};

	fds[0].fd = vgp_sub;

	fds[0].events = POLLIN;

	usleep(10000000);// 10 seconds break to get valid data from gps
	while (!thread_should_exit) {

		/*
		 * Wait for a sensor or param update, check for exit condition every 500 ms.
		 * This means that the execution will block here without consuming any resources,
		 * but will continue to execute the very moment a new attitude measurement or
		 * a param update is published. So no latency in contrast to the polling
		 * design pattern (do not confuse the poll() system call with polling).
		 *
		 * This design pattern makes the controller also agnostic of the attitude
		 * update speed - it runs as fast as the attitude updates with minimal latency.
		 */
		int ret = poll(fds, 1, 1000);

		if (ret < 0) {
			/*
			 * Poll error, this will not really happen in practice,
			 * but its good design practice to make output an error message.
			 */
			warnx("poll error");

		} else if (ret == 0) {
			/* no return value = nothing changed for 500 ms, ignore */
		} else {

			if (fds[0].revents & POLLIN) {

				//this section runs  when a new estimate from vehicle_global_position arrives

                                if (hrt_elapsed_time(&then)>1_s)
				{
					orb_copy(ORB_ID(vehicle_global_position),vgp_sub,&vgp);


					//Function 1
					// for this hardware parts would be needed.
					// First we check for RPAS identifier, if any device has been changed or not
					// there is a true value for drone id and a value calculated every time powering on the system
					// if they both matched,then RPAS has not been tampered
					// if they didnt , RPAS has been tampered , not allow to fly, A file is created and send to
					//management client (signed by key/pair_C), then from management client to Management server.
					//file informaton: 1)timestamp 2)Component's device id that didnt match.

					static int RPAS_identifier_check=0;
						return 0;
					if (RPAS_identifier_check==0){
						FILE *fptr;
						fptr=fopen("/fs/microsd/gps_id.txt","w");
						int identity_check= RPAS_identifier(fptr);
						fprintf(fptr,"ide  %d",identity_check);
						fclose(fptr);
					

						/*if(identity_check==0){

							mavlink_log_critical(mavlink_log_pub, "Unautorized Hardware parts attached");

							return false;
						}
						if(identity_check==2){
							// this is the case when HArdwareInuse.tx fie is removed from the firmware
							mavlink_log_critical(mavlink_log_pub, "Hardware refference file not found. Please update the firmware");
							return false;
						}
						if(identity_check==1){
							// all hardwares are authorized
							RPAS_identifier_check=1;

						}*/
					}



					then=hrt_absolute_time();
				}


				// check if landing has taken place or not





			}
		}
	}




	thread_running = false;


	fflush(stdout);

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
int initial_checks_main(int argc, char *argv[])
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
        usleep(10000000);
	int indi_sub_fd = orb_subscribe(ORB_ID(indi_prc));
        while(1){
			struct indi_prc_s indi;
			/* copy sensors raw data into local buffer */
			orb_copy(ORB_ID(indi_prc), indi_sub_fd, &indi);
			if(indi.valid==1){
				break;
			}
			usleep(2000000);
	}
		thread_should_exit = false;
		deamon_task = px4_task_spawn_cmd("initial_checks",
						 SCHED_DEFAULT,
						 SCHED_PRIORITY_MAX - 20,
						 8048,
						 initial_checks_thread_main,
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



