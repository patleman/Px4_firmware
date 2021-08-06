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
#include "key_rotation.hpp"

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
extern "C" __EXPORT int key_rotation_main(int argc, char *argv[]);

/**
 * Mainloop of daemon.
 */
int key_rotation_thread_main(int argc, char *argv[]);

/**
 * Print the correct usage.
 */
static void usage(const char *reason);


/* Variables */
static bool thread_should_exit = false;		/**< Daemon exit flag */
static bool thread_running = false;		/**< Daemon status flag */
static int deamon_task;				/**< Handle of deamon task / thread */



/* Main Thread */
int key_rotation_thread_main(int argc, char *argv[])
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
     
// Key rotation start, if needed

	FILE *fptr_key;
	fptr_key=fopen("/fs/microsd/KeyRotation.txt","r");

	if(fptr_key!=NULL)
	{	// if present, then check for its validity


		fclose(fptr_key);
		char DroneID_container_2[30];
		char FILE_ID_container[30];
		char TAG_DRONE_ID_2[10]="DroneID";
		char TAG_KEY_ROT_ID[20]="FILE_ID";
		char DRONE_ID_FILE_2[100]="/fs/microsd/KeyRotation.txt";
		int check_validity_DroneID_2=file_read(DRONE_ID_FILE_2,TAG_DRONE_ID_2,DroneID_container_2,1);
		// if non tampered -->then continue
	//	check_validity_DroneID_2=file_read(DRONE_ID_FILE_2,TAG_KEY_ROT_ID,FILE_ID_container,1);
        fetch_tag(NULL,TAG_KEY_ROT_ID,FILE_ID_container,DRONE_ID_FILE_2);
		
		if (check_validity_DroneID_2==0){
			//KeyRotation.txt file is not valid, remove the invalid file
			remove("/fs/microsd/KeyRotation.txt");


		}else{
			//file is valid, check for strings
			if(strcmp(DroneID_container_2,"ABBDDJEDNDJK")==0){
				//Checking for reusage of the file
				int reusage_check=CHECK_REUSAGE(FILE_ID_container);
				if (reusage_check==0){
					//File_ID is new, first time used
					//Update KeyLog.txt with the new File_ID
					call_KeyLog_Regen(FILE_ID_container);
					// start key rotation
					Key_rotation_start(FILE_ID_container);
					remove("/fs/microsd/KeyRotation.txt");
					
				}else{// The file for key rotation has already been used
				      // key rotation will not take place
					remove("/fs/microsd/KeyRotation.txt");

				}

			}else{
				// KeyRotation.txt file is for different drone
				//
				remove("/fs/microsd/KeyRotation.txt");
			}

		}

	}else
	{
			//no need of key rotation as there is no keyRotation.txt file present.
			//check for KeyChangePerm.txt
			FILE *fptr_keychange;
			char fname_keychange[40]="/fs/microsd/KeyChangePerm.txt";

			fptr_keychange=fopen(fname_keychange,"r");
			if(fptr_keychange!=NULL){
			fclose(fptr_keychange);
			//file present, check fot its validity, drone_id, KEY_ID
			//if all okay--> then make the change by calling the function
			// if not okay-->reasons 1)File is for different drone
			//			 2)File is not valid
			//                       3)file KEY_ID doesnt match
			//two files in play PublicKeyNew.txt and KeyChangePerm.txt
			char DroneID_container_3[30];
			char KEY_ID_container[30];
			char TAG_DRONE_ID_2[10]="DroneID";
			char TAG_KEY_ID[20]="KEY_ID";
			char MODULUS_FILE[720];
			char MODULUS_tag[20]="Modulus";
			char DRONE_ID_FILE_2[60]="/fs/microsd/KeyChangePerm.txt";
			int check_validity_DroneID_3=file_read(DRONE_ID_FILE_2,TAG_DRONE_ID_2,DroneID_container_3,1);
			// if non tampered -->then continue
			//check_validity_DroneID_3=file_read(DRONE_ID_FILE_2,TAG_KEY_ID,KEY_ID_container,1);

			//check_validity_DroneID_3=file_read(DRONE_ID_FILE_2,MODULUS_tag,MODULUS_FILE,1);
			fetch_tag(NULL,TAG_KEY_ID,KEY_ID_container,DRONE_ID_FILE_2);
			fetch_tag(NULL,MODULUS_tag,MODULUS_FILE,DRONE_ID_FILE_2);
			
			if(check_validity_DroneID_3==1){
				//file is valid
				//Now check for DroneID,KEY_ID and MODULUS
				//taking above from PublicKeyNew.txt
				char PKN_Modulus[720];
				char PKN_DroneID[30];
				char PKN_KEY_ID[30];
				char PKN_file[60]="/fs/microsd/PublicKeyNew.txt";
				check_validity_DroneID_3=file_read(PKN_file,MODULUS_tag,PKN_Modulus,0);
			//	check_validity_DroneID_3=file_read(PKN_file,TAG_DRONE_ID_2,PKN_DroneID,0);
			//	check_validity_DroneID_3=file_read(PKN_file,TAG_KEY_ID,PKN_KEY_ID,0);
				fetch_tag(NULL,TAG_KEY_ID,PKN_KEY_ID,PKN_file);
			    fetch_tag(NULL,TAG_DRONE_ID_2,PKN_DroneID,PKN_file);
			
				if(check_validity_DroneID_3==0){
					//not valid PublicKeyNew.txt, remove it
					remove("/fs/microsd/PublicKeyNew.txt");
				}else{
					int valid_sum=0;
					if(strcmp(PKN_Modulus,MODULUS_FILE)==0){
						valid_sum++;
					}
					if(strcmp(PKN_DroneID,DroneID_container_3)==0){
						valid_sum++;
					}
					if(strcmp(PKN_KEY_ID,KEY_ID_container)==0){
						valid_sum++;
					}
					if(valid_sum==3){
						// file is valid and belongs to this drone
						// begin the changes that are needed
						KEY_CHANGE_INITIATION();
						

					}else{	//
						//file is valid but does not belong to this drone.
						remove("/fs/microsd/KeyChangePerm.txt");
					}
				}

			}else{
				//file is not valid, remove the invalid file
				remove("/fs/microsd/KeyChangePerm.txt");
			}

		}else{//file is not present, no need to change keys

		}

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
int key_rotation_main(int argc, char *argv[])
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
		deamon_task = px4_task_spawn_cmd("key_rotation",
						 SCHED_DEFAULT,
						 SCHED_PRIORITY_MAX - 20,
						 35048,
						 key_rotation_thread_main,
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



