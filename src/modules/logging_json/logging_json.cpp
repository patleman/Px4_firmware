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
#include "logging_json.h"

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
extern "C" __EXPORT int logging_json_main(int argc, char *argv[]);

/**
 * Mainloop of daemon.
 */
int logging_json_thread_main(int argc, char *argv[]);

/**
 * Print the correct usage.
 */
static void usage(const char *reason);


/* Variables */
static bool thread_should_exit = false;		/**< Daemon exit flag */
static bool thread_running = false;		/**< Daemon status flag */
static int deamon_task;				/**< Handle of deamon task / thread */



/* Main Thread */
int logging_json_thread_main(int argc, char *argv[])
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




	/*
	 * Declare and safely initialize all structs to zero.
	 *
	 * These structs contain the system state and things
	 * like attitude, position, the current waypoint, etc.
	 */


	Geo_tag *content_holder=(Geo_tag*) malloc(1500*sizeof(Geo_tag));



	//this is for knowing current time
	struct vehicle_gps_position_s vgps;
	memset(&vgps,0,sizeof(vgps));

	//vehicle global postion: got after estimation
	// this will be compared with the coordinates taken from recentPA.txt
	// (2x(lat,long), alt)
	struct vehicle_global_position_s vgp;
	memset(&vgp, 0, sizeof(vgp));

	//home position : this is to take reference ground level,
	//later helpful for knowing geofence area and allowed time
	// later helpful for altitude comparison
	struct pa_data_s pad;
	memset(&pad, 0, sizeof(pad));

	//takeoff status: to note down the instance at takeoff
	//to know when take off is taking place
	struct takeoff_status_s TO_status;
	memset(&TO_status, 0, sizeof(TO_status));


	//vehicle land detected : to note down the instance at landing
	// to know when landing is taking place
	struct vehicle_land_detected_s VLD;
	memset(&VLD, 0, sizeof(VLD));

	// this is to initiate RTL when geo breach/time breach is experienced
	//
	//struct vehicle_command_s v_cmd;
	//memset(&v_cmd, 0, sizeof(v_cmd));





	/*
	 * Later used to Advertise RTL vehicle command
	 * */
	//orb_advert_t v_cmd_pub = orb_advertise(ORB_ID(vehicle_command), &v_cmd);

        static int pas_time=0;
	/* subscribe to topics. */
	int vgp_sub = orb_subscribe(ORB_ID(vehicle_global_position));

	int vgps_sub = orb_subscribe(ORB_ID(vehicle_gps_position));

	int pad_sub = orb_subscribe(ORB_ID(pa_data));

	int TO_status_sub = orb_subscribe(ORB_ID(takeoff_status));

	int VLD_sub = orb_subscribe(ORB_ID(vehicle_land_detected));

	orb_set_interval(vgp_sub,1000);
	orb_set_interval(vgps_sub,1000);
	/* Setup of loop */

	struct pollfd fds[1] {};

	fds[0].fd = vgp_sub;

	fds[0].events = POLLIN;

//////////////////////////////////////////////////////
reset_check:
	thread_should_exit=0;
	int rtl_pass=0;
	int geo_tag_count=0;
	memset(content_holder,0,1500);

	// Now control will go into the  second loop only when two conditions are met, which are
	// 1)pa_data is updated
	// 2)takeoff_state is 5


	while(1){
		orb_copy(ORB_ID(pa_data),pad_sub,&pad);
		orb_copy(ORB_ID(takeoff_status),TO_status_sub,&TO_status);
		//printf(",hey patlee\n");

		if((pad.updated==1)&&(TO_status.takeoff_state==5)){
			break;
		}
		usleep(1000000);

	}
        printf("\n\nwhile loop break\n\n");

	if(TO_status.takeoff_state==5){

		orb_copy(ORB_ID(vehicle_global_position),vgp_sub,&vgp);// for lat,lon,altitude
		orb_copy(ORB_ID(vehicle_gps_position),vgps_sub,&vgps);// for timestamp
		// notedown the takeoff instance in log file
		printf("\n\nwhile loop break  66\n\n");
		content_holder[geo_tag_count].Entrytype=1;
		content_holder[geo_tag_count].Timestamp=(vgps.time_utc_usec)/1000000;
		content_holder[geo_tag_count].Lattitude=vgp.lat;
		content_holder[geo_tag_count].Longitude=vgp.lon;
		content_holder[geo_tag_count].Altitude=vgp.alt;
		geo_tag_count++;
		memset(&vgps,0,sizeof(vgps));
		memset(&vgp,0,sizeof(vgp));


	}
	hrt_abstime then=hrt_absolute_time();
       printf("\n\ntakeoff detected block passed %d\n\n",TO_status.takeoff_state);
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

                                if (hrt_elapsed_time(&then)>1_s){
					orb_copy(ORB_ID(vehicle_global_position),vgp_sub,&vgp);
					orb_copy(ORB_ID(vehicle_gps_position),vgps_sub,&vgps);
					orb_copy(ORB_ID(pa_data),pad_sub,&pad);

					int geobreach_status=check_Geobreach(pad,vgp);
					int timebreach_status=check_Timebreach(pad,vgps);


					//check for geo_breach
					if(geobreach_status==1){
						pas_time++;
						content_holder[geo_tag_count].Entrytype=0;
						content_holder[geo_tag_count].Timestamp=(vgps.time_utc_usec)/1000000;
						content_holder[geo_tag_count].Lattitude=vgp.lat;
						content_holder[geo_tag_count].Longitude=vgp.lon;
						content_holder[geo_tag_count].Altitude=vgp.alt;
						geo_tag_count++;

					}
					//check for time breach
					if(timebreach_status==1){
						content_holder[geo_tag_count].Entrytype=2;
						content_holder[geo_tag_count].Timestamp=(vgps.time_utc_usec)/1000000;
						content_holder[geo_tag_count].Lattitude=vgp.lat;
						content_holder[geo_tag_count].Longitude=vgp.lon;
						content_holder[geo_tag_count].Altitude=vgp.alt;
						geo_tag_count++;

					}
					//if(timebreach_status || geobreach_status){

					memset(&vgps,0,sizeof(vgps));
					memset(&vgp,0,sizeof(vgp));
					memset(&pad,0,sizeof(pad));

					//}
					if((timebreach_status || geobreach_status) && (!rtl_pass) ){
						//send RTL command

						vehicle_command_s vcmd{};
						//vcmd.command = vehicle_command_s::VEHICLE_CMD_NAV_RETURN_TO_LAUNCH;
						vcmd.command =vehicle_command_s::VEHICLE_CMD_DO_SET_MODE;
						vcmd.param1 = 1;
						vcmd.param2 = 4;
						vcmd.param3 = 5;


						uORB::SubscriptionData<vehicle_status_s> vehicle_status_sub{ORB_ID(vehicle_status)};
						vcmd.source_system = vehicle_status_sub.get().system_id;
						vcmd.target_system = vehicle_status_sub.get().system_id;
						vcmd.source_component = vehicle_status_sub.get().component_id;
						vcmd.target_component = vehicle_status_sub.get().component_id;

						uORB::Publication<vehicle_command_s> vcmd_pub{ORB_ID(vehicle_command)};
						vcmd.timestamp = hrt_absolute_time();
						vcmd_pub.publish(vcmd);
						rtl_pass=1;
					}
					then=hrt_absolute_time();
				}


				// check if landing has taken place or not
				orb_copy(ORB_ID(vehicle_land_detected), VLD_sub, &VLD);
				char bnm[4]="hi";
				if(VLD.landed){

					orb_copy(ORB_ID(vehicle_global_position),vgp_sub,&vgp);// for lat,lon,altitude
					orb_copy(ORB_ID(vehicle_gps_position),vgps_sub,&vgps);// for timestamp
					// landing has taken place
					//notedown the instance in log
					content_holder[geo_tag_count].Entrytype=3;
					content_holder[geo_tag_count].Timestamp=(vgps.time_utc_usec)/1000000;
					content_holder[geo_tag_count].Lattitude=vgp.lat;
					content_holder[geo_tag_count].Longitude=vgp.lon;
					content_holder[geo_tag_count].Altitude=vgp.alt;
					geo_tag_count++;

					memset(&vgps,0,sizeof(vgps));
					memset(&vgp,0,sizeof(vgp));
                                        printf("\nhey its landing\n");

					// exit from the loop
					update_recentPA(1,bnm);
					 printf("\nhey recent PA updaed after landing\n");
					thread_should_exit=1;
					// update the flight frequency


				}
				if(VLD.freefall){
					// drone is in freefall state
					// notedown the instance in log
					orb_copy(ORB_ID(vehicle_global_position),vgp_sub,&vgp);// for lat,lon,altitude
					orb_copy(ORB_ID(vehicle_gps_position),vgps_sub,&vgps);// for timestamp

					//notedown the instance in log
					content_holder[geo_tag_count].Entrytype=4;
					content_holder[geo_tag_count].Timestamp=vgps.time_utc_usec;
					content_holder[geo_tag_count].Lattitude=vgp.lat;
					content_holder[geo_tag_count].Longitude=vgp.lon;
					content_holder[geo_tag_count].Altitude=vgp.alt;
					geo_tag_count++;

					memset(&vgps,0,sizeof(vgps));
					memset(&vgp,0,sizeof(vgp));



					//exit from loop
					thread_should_exit=1;
					//updating the frequency
					update_recentPA(1,bnm);
				}
				memset(&VLD,0,sizeof(VLD));




			}
		}
	}
        printf("\n\ncount geobreach :::::::%d\n\n",pas_time);

        main_json_file_writing(content_holder,geo_tag_count);
	free(content_holder);

	//Bundling_begins();

	goto reset_check;




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
int logging_json_main(int argc, char *argv[])
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
		deamon_task = px4_task_spawn_cmd("logging_json",
						 SCHED_DEFAULT,
						 SCHED_PRIORITY_MAX - 20,
						 3548,
						 logging_json_thread_main,
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



