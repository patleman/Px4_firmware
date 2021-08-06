/****************************************************************************
 *
 *   Copyright (c) 2013-2017 PX4 Development Team. All rights reserved.
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
 * @file main.c
 *
 * Example implementation of a fixed wing attitude controller. This file is a complete
 * fixed wing controller for manual attitude control or auto waypoint control.
 * There is no need to touch any other system components to extend / modify the
 * complete control architecture.
 *
 * @author Lorenz Meier <lm@inf.ethz.ch>
 */

#include "params.h"

#include <poll.h>

#include <drivers/drv_hrt.h>
#include <lib/ecl/geo/geo.h>
#include <matrix/math.hpp>
#include <px4_platform_common/px4_config.h>
#include <px4_platform_common/tasks.h>
#include <systemlib/err.h>
#include <parameters/param.h>
#include <perf/perf_counter.h>
#include <uORB/Subscription.hpp>
#include <uORB/SubscriptionInterval.hpp>
#include <uORB/topics/actuator_controls.h>
#include <uORB/topics/manual_control_setpoint.h>
#include <uORB/topics/parameter_update.h>
#include <uORB/topics/position_setpoint_triplet.h>
#include <uORB/topics/vehicle_attitude.h>
#include <uORB/topics/initial_check_status.h>
#include <uORB/topics/vehicle_attitude_setpoint.h>
#include <uORB/topics/vehicle_global_position.h>
#include <uORB/topics/vehicle_global_position.h>
#include <uORB/topics/vehicle_rates_setpoint.h>
#include <uORB/topics/vehicle_status.h>

////////////////

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<inttypes.h>

#include <uORB/topics/takeoff_status.h>
#include <uORB/topics/home_position.h>
#include <uORB/topics/vehicle_command.h>
#include <uORB/topics/vehicle_land_detected.h>
#include <uORB/topics/vehicle_global_position.h>
#include <uORB/topics/vehicle_gps_position.h>
#include <uORB/topics/vehicle_status.h>
#include <uORB/topics/pa_data.h>
#include <uORB/topics/indi_prc.h>
#include <drivers/drv_hrt.h>
#include<motion_planning/main_utility.hpp>
#include<motion_planning/print_hello.hpp>

/////////////
#include <sys/random.h> 
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

#include <examples/logging_json/logging_json.hpp>

#include <uORB/topics/parameter_update.h>
#include <parameters/param.h>
#include <lib/ecl/geo/geo.h>
#include <perf/perf_counter.h>
#include <systemlib/err.h>
#include <matrix/math.hpp>


////////////////

using namespace time_literals;

/**
 * Daemon management function.
 *
 * This function allows to start / stop the background task (daemon).
 * The purpose of it is to be able to start the controller on the
 * command line, query its status and stop it, without giving up
 * the command line to one particular process or the need for bg/fg
 * ^Z support by the shell.
 */
extern "C" __EXPORT int ex_fixedwing_control_main(int argc, char *argv[]);

/**
 * Mainloop of daemon.
 */
int fixedwing_control_thread_main(int argc, char *argv[]);

/**
 * Print the correct usage.
 */
static void usage(const char *reason);

/* Variables */
static int deamon_task;				/**< Handle of deamon task / thread */


static bool thread_should_exit = false;		/**< Daemon exit flag */
static bool thread_running = false;		/**< Daemon status flag */

/* Main Thread */
int fixedwing_control_thread_main(int argc, char *argv[])
{
/*	Geo_tag content_holder[10];
	memset(&content_holder[0],0,sizeof(content_holder));

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





	
	 // Later used to Advertise RTL vehicle command
	 
	//orb_advert_t v_cmd_pub = orb_advertise(ORB_ID(vehicle_command), &v_cmd);

        static int pas_time=0;
	// subscribe to topics. 
	int vgp_sub = orb_subscribe(ORB_ID(vehicle_global_position));

	int vgps_sub = orb_subscribe(ORB_ID(vehicle_gps_position));

	int pad_sub = orb_subscribe(ORB_ID(pa_data));

	int TO_status_sub = orb_subscribe(ORB_ID(takeoff_status));

	int VLD_sub = orb_subscribe(ORB_ID(vehicle_land_detected));

	orb_set_interval(vgp_sub,1000);
	orb_set_interval(vgps_sub,1000);
	

	struct pollfd fds[1] {};

	fds[0].fd = vgp_sub;

	fds[0].events = POLLIN;
*/

/*
char file_id[20]="rt567";
Key_rotation_start(file_id);*/
  /*script for generating prime number
  FILE *file;
  file=fopen("/fs/microsd/key_new.txt","w");
   fprintf(file,"\n%s\n","prime number generation");
    mp_int  p1;
	mp_int q1;
    char    buf[1000];

    mp_init(&p1);
    mp_init(&q1);

    mp_prime_rand(&p1,1,1024,MP_PRIME_2MSB_ON,file);
    mp_to_decimal(&p1,buf,sizeof(buf));
    fprintf(file,"\nPrime 1 :%s\n",buf);
    memset(buf,0,sizeof(buf));
    
	mp_prime_rand(&q1,1,1024,MP_PRIME_2MSB_ON,file);
    mp_to_decimal(&q1,buf,sizeof(buf));
    fprintf(file,"\n\nPrime 2 :%s\n",buf);
    memset(buf,0,sizeof(buf)); 
	mp_clear(&p1);
	mp_clear(&q1);*/

//usleep(3000000);
  /* script for testing file method for random genration
    FILE *file;
    file=fopen("/fs/microsd/key_new.txt","w");
    
    FILE *f;
	f=fopen("/fs/microsd/log/pprime.dat","r");
	if(f){
		fprintf(file,"\n%s\n","opening....");
		fclose(f);
	}else{
		fprintf(file,"\n%s\n","not opening....");
	}
   mp_int  p1,q1;
   char    buf[1096];
   int     k=40;
   int  li=3;
   
   

 

   load_tab();
   int err;
   err=mp_init(&p1);
   err=mp_init(&q1);

   
  err= pprime(k, li, &p1, &q1);
   int r=rand();
   fprintf(file,"\nrandom number generated == %d\n", r);
  err= mp_to_decimal(&p1, buf, sizeof(buf));
   fprintf(file,"\nP == %s\n", buf);
  err= mp_to_decimal(&q1, buf, sizeof(buf));
   fprintf(file,"Q == %s\n", buf);
    mp_clear(&q1);
	mp_clear(&p1);

    mp_int s;
	mp_init(&s);
	unsigned char buff[50]="100010010101";
	err=mp_from_ubin(&s,buff, sizeof(buf));
	mp_read_radix(&s,buf,10);
	fprintf(file,"err code == %d\n", err);
	fprintf(file,"decimal of binary == %s\n", buf);

	mp_clear(&s);*/
 // lib_dumpbuffer("Random values", (FAR const uint8_t*)buffer, nread);
  //close(fd);
 
	//fclose(file);

//pa_validation();
/*
key key_rfm;

FILE *fptr1;

fptr1=fopen("/fs/microsd/debug.txt","w");
get_RFM_Key(&key_rfm);
fprintf(fptr1,"\n%s\n",key_rfm.modulus);
fprintf(fptr1,"\n%s\n",key_rfm.private_exponent);

char f[60]="/fs/microsd/log/recentPA.txt";
int valid_status=Validating_File(f,key_rfm,fptr1);

fprintf(fptr1,"\n \n status :%d\n",valid_status);
//fetching_publish_padata();

fclose(fptr1);	*/			
 /*FILE *fptr;
      fptr=fopen("/fs/microsd/result.txt","w");

      char file_name[]="/fs/microsd/log/permission_artifact_breach.xml";
      char *Digest_Value0=(char*) malloc(100*sizeof(char));
      char Sha_of_Reference[70];
      char Sha_of_SignedInfo[70];
      valid_chunk_refe(fptr,Sha_of_Reference,1);
      fprintf(fptr,"\n\n%s\n",Sha_of_Reference);

      valid_chunk_refe(fptr,Sha_of_SignedInfo,2);
      fprintf(fptr,"\n\n%s\n",Sha_of_SignedInfo);

      char Tag_Digest_value0[12]="DigestValue";
      getTagvalue(Tag_Digest_value0,Digest_Value0,file_name);
	  fprintf(fptr,"\n\n Digest value %s\n",Digest_Value0);

      char *Digest_Value_in_hex=(char*) malloc(100*sizeof(char));
      base64decoder(Digest_Value0, Digest_Value_in_hex);
      fprintf(fptr,"\n\n Digest value %s\n",Digest_Value_in_hex);
      free(Digest_Value0);



         char Signature_Value_in_hex[515];//;=(char*) malloc(500*sizeof(char));
         char Modulus[514];
            //  mp_int ; dgca public key modulus
				strcpy(Modulus,"d6912a773335d3a193c742071762794e26bcac49d78b6b65a784e1c18d18f2f88b9f13d7fc41a5b9e11e3d75905b0f8373dbac5658962940d104711467814b7319774014644df82e9764b4e852cb7547d56de3224aada000db231b1356ffe86391b3eca2ddf67d44f40ea6e84092cc67387db3a2042487b8dacfe5b588738973a59f5fa813a33e4bb8fbafa407794db990a307201b0a0fbd92ddd181a868a6cfa64625353714e9b1de6f81f3addbf5b202030dbd9e51385db9314591f01af01f07cc92f18ebe60545cea53eb54438e251fd88e7380b5a0612c29ccb7dfd88f7a7d7f0dd07a9b596923602798ad9f1bbb89fffd3cdb8fb94c48640d27fa681e47");
				// char Modulus[513]="ab9d5c8d1fe67207749d63b7dcedd233ce32bb70d175a1bc38c612ab33e2c58e51f83f2788e4d52d9bceb5a1513929de3f526650071a067e6c161b05c60a495fc3ba79ed26f4fa8b2fe2ca8dec44b39759f39206f06a85f9424005a29f05e4cf3a0239340c28c993c1a61cf1b2b6b57c7d8e576ae86827f812b327625baec9ecbf55f1651d35600b9f955f6c2f3bea3aa5852ecdd36a0af818c19acc1030979bed3c89993faa92e0aa0502413b3ca86bbf63477f12ac069aff7137cb72c57f886da79033bbb3b4df0f6cc7fcc18e343aa76036681a566311e267c03b65c98abc91e58f090020c67f776199c0eb76d7e6363687475d3da36ff050f85275607fdd";
            // then processing is of PA
          //  char file_name[]="/fs/microsd/log/permission_artifact_breach.xml";
	         char Tag_Signed_Value0[17]="SignatureValue";
				//  char Signature_Value0[550];// in base 64
				char *Signature_Value0=(char*) malloc(550*sizeof(char));

				getTagvalue(Tag_Signed_Value0,Signature_Value0, file_name);
 fprintf(fptr,"\n\n Digest value %s\n",Signature_Value0);
				base64decoder(Signature_Value0, Signature_Value_in_hex);
 fprintf(fptr,"\n\n Digest value %s\n",Signature_Value_in_hex);
				free(Signature_Value0);

				int cb;
				mp_int message,modulus,public_key,Decrypted;
				cb= mp_init_multi(&message,&modulus,&public_key,&Decrypted,NULL);
				cb= mp_read_radix(&message,Signature_Value_in_hex,16);
				// free(Signature_Value_in_hex);

				/// Public key of dgca: This has to be taken from a reserved file in directory./ firmware update

				cb= mp_read_radix(&modulus,Modulus,16);

				// mp_int ;
				cb= mp_read_radix(&public_key,"65537",10);

				//mp_int ;
				cb= mp_exptmod(&message, &public_key,&modulus,&Decrypted);   //                       this part for decrypting encrypted text
				printf("%d",cb);
            char message_string_hex[513];
            mp_to_hex(&Decrypted,message_string_hex,sizeof(message_string_hex));
            mp_clear_multi(&message,&modulus,&public_key,&Decrypted,NULL);

              fprintf(fptr,"\n\n Digest mesa %s\n",message_string_hex);
      char *Extracted_sha_from_Signature_value=(char*) malloc(100*sizeof(char));
      char *messa=(char*) malloc(750*sizeof(char));

     check_debug_mp(messa);

	 fprintf(fptr,"\n\n Digest mesa 2 %s\n",messa);
	 char padding_SHA256[447]="1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff003031300d060960864801650304020105000420";
      for(int o=0;o<444;o++){
            if(padding_SHA256[o]!=small_letter(messa[o])){
               //return 2;
               fprintf(fptr,"%s","\n Bad SIGN \n");
               break;
            }
      }

      fprintf(fptr,"%s","\n SIGN is good \n");
      int k_cout=0;
      for(int i2=445;messa[i2]!='\0';i2++){
            Extracted_sha_from_Signature_value[k_cout]=small_letter(messa[i2]);
            k_cout++;
         // printf("%c",message_string_hex[i]);
      }
      Extracted_sha_from_Signature_value[k_cout]='\0';//
      fprintf(fptr,"\n%s\n",Extracted_sha_from_Signature_value);
      if(strcmp(Extracted_sha_from_Signature_value,Sha_of_SignedInfo)==0){

         if(strcmp(Sha_of_Reference,Digest_Value_in_hex)==0){
                  fprintf(fptr,"%s","\n 100000000000PA is valid and non tampered \n");

            }else{
                  fprintf(fptr,"%s","\n PA is not valid\n");

            }

         }else{
            fprintf(fptr,"%s","\n PA is not valid \n");
         //  PX4_INFO("PA is valid???????????");

         }
      free(Extracted_sha_from_Signature_value);
      free(messa);
      free(Digest_Value_in_hex);
    

	 fclose(fptr);*/
/*
      char padding_SHA256[447]="1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff003031300d060960864801650304020105000420";
      for(int o=0;o<444;o++){
            if(padding_SHA256[o]!=small_letter(message[o])){
               //return 2;
               fprintf(fptr,"%s","\n Bad SIGN \n");
               break;
            }
      }

      fprintf(fptr,"%s","\n SIGN is good \n");
      int k_cout=0;
      for(int i2=445;message[i2]!='\0';i2++){
            Extracted_sha_from_Signature_value[k_cout]=small_letter(message[i2]);
            k_cout++;
         // printf("%c",message_string_hex[i]);
      }
      Extracted_sha_from_Signature_value[k_cout]='\0';//
      fprintf(fptr,"\n%s\n",Extracted_sha_from_Signature_value);
      if(strcmp(Extracted_sha_from_Signature_value,Sha_of_SignedInfo)==0){

         if(strcmp(Sha_of_Reference,Digest_Value_in_hex)==0){
                  fprintf(fptr,"%s","\n 100000000000PA is valid and non tampered \n");

            }else{
                  fprintf(fptr,"%s","\n PA is not valid\n");

            }

         }else{
            fprintf(fptr,"%s","\n PA is not valid \n");
         //  PX4_INFO("PA is valid???????????");

         }
      free(Extracted_sha_from_Signature_value);
      free(message);
      free(Digest_Value_in_hex);
      fclose(fptr);*/
	  
//usleep(3000000);

//fetching_publish_padata();
/*struct indi_prc_s indi;
memset(&indi, 0, sizeof(indi));
orb_advert_t indi_pub = orb_advertise(ORB_ID(indi_prc), &indi);
indi.valid=1;
orb_publish(ORB_ID(indi_prc), indi_pub, &indi);*/
/*key RFM_key;
get_RFM_Key(&RFM_key);
char fname[60]="/fs/microsd/log/HardwareInuse.txt";
char res[40];
char tag[30]="GPS_ID";
fetch_tag(NULL, tag,res, fname);
int identity_check= Validating_File(fname,&RFM_key);
FILE *fptr;
fptr=fopen("/fs/microsd/hardware_valid.txt","w");
fprintf(fptr,"\nthe validation status : %d\n",identity_check);
fprintf(fptr,"\ntag gps_id %s\n",res);
identity_check= RPAS_identifier(fptr);
fprintf(fptr,"\n\n\nthe validation status : %d\n",identity_check);
fprintf(fptr,"\ntag gps_id %s\n",res);
fclose(fptr);*/


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
	{//file is present
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
//////////////////////////////////////////////////////
/*reset_check:
	thread_should_exit=0;
	int rtl_pass=0;
	int geo_tag_count=0;
	//memset(content_holder,0,1500);

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

	orb_copy(ORB_ID(takeoff_status),TO_status_sub,&TO_status);
//phirse_check:
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


	}*//*else{
		usleep(5000000);
		goto phirse_check;
		//px4_usleep(10000);

	}*/
	/*
	hrt_abstime then=hrt_absolute_time();
       printf("\n\ntakeoff detected block passed %d\n\n",TO_status.takeoff_state);
	while (!thread_should_exit) {*/

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
	/*	int ret = poll(fds, 1, 1000);

		if (ret < 0) {
			
			 //Poll error, this will not really happen in practice,
			 // but its good design practice to make output an error message.
			 
			warnx("poll error");

		} else if (ret == 0) {
			 //no return value = nothing changed for 500 ms, ignore 
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
					//content_holder[geo_tag_count].Entrytype=3;
					//content_holder[geo_tag_count].Timestamp=(vgps.time_utc_usec)/1000000;
					//content_holder[geo_tag_count].Lattitude=vgp.lat;
					//content_holder[geo_tag_count].Longitude=vgp.lon;
					//content_holder[geo_tag_count].Altitude=vgp.alt;
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
					//content_holder[geo_tag_count].Entrytype=4;
					//content_holder[geo_tag_count].Timestamp=vgps.time_utc_usec;
					//content_holder[geo_tag_count].Lattitude=vgp.lat;
					//content_holder[geo_tag_count].Longitude=vgp.lon;
					//content_holder[geo_tag_count].Altitude=vgp.alt;
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

      //  main_json_file_writing(content_holder,geo_tag_count);
	//free(content_holder);

	//Bundling_begins();

	goto reset_check;




	thread_running = false;




	fflush(stdout);

	//return 0;*/
}

/* Startup Functions */

static void
usage(const char *reason)
{
	if (reason) {
		fprintf(stderr, "%s\n", reason);
	}

	fprintf(stderr, "usage: ex_fixedwing_control {start|stop|status}\n\n");
}

/**
 * The daemon app only briefly exists to start
 * the background job. The stack size assigned in the
 * Makefile does only apply to this management task.
 *
 * The actual stack size should be set in the call
 * to px4_task_spawn_cmd().
 */
int ex_fixedwing_control_main(int argc, char *argv[])
{
	if (argc < 2) {
		usage("missing command");
		return 1;
	}

	if (!strcmp(argv[1], "start")) {

		if (thread_running) {
			printf("ex_fixedwing_control already running\n");
			/* this is not an error */
			return 0;
		}

		thread_should_exit = false;
		deamon_task = px4_task_spawn_cmd("ex_fixedwing_control",
						 SCHED_DEFAULT,
						 SCHED_PRIORITY_MAX - 20,
						 25048,
						 fixedwing_control_thread_main,
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
			printf("\tex_fixedwing_control is running\n");

		} else {
			printf("\tex_fixedwing_control not started\n");
		}

		return 0;
	}

	usage("unrecognized command");
	return 0;
}
