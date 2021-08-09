/****************************************************************************
 *
 *   Copyright (c) 2019-2020 PX4 Development Team. All rights reserved.
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
 * @file PreFlightCheck.cpp
 */

#include "PreFlightCheck.hpp"
//#include "MP_INT.hpp"
#include <motion_planning/main_utility.hpp>
#include <drivers/drv_hrt.h>
#include <HealthFlags.h>

#include <px4_defines.h>
#include <lib/parameters/param.h>
#include <systemlib/mavlink_log.h>
#include <uORB/Subscription.hpp>
#include <sys/random.h>

#include <uORB/topics/vehicle_gps_position.h>

#include <uORB/topics/sensor_gps.h>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include<stdlib.h>
#include<time.h>


using namespace time_literals;


///////////////////////////////////////////
static constexpr unsigned max_mandatory_gyro_count = 1;
static constexpr unsigned max_optional_gyro_count = 4;
static constexpr unsigned max_mandatory_accel_count = 1;
static constexpr unsigned max_optional_accel_count = 4;
static constexpr unsigned max_mandatory_mag_count = 1;
static constexpr unsigned max_optional_mag_count = 4;
static constexpr unsigned max_mandatory_baro_count = 1;
static constexpr unsigned max_optional_baro_count = 4;




bool PreFlightCheck::preflightCheck(orb_advert_t *mavlink_log_pub, vehicle_status_s &status,
				    vehicle_status_flags_s &status_flags, bool report_failures, const bool prearm,
				    const hrt_abstime &time_since_boot)





{



static int RPAS_identifier_check=0;

if (RPAS_identifier_check==0){

int identity_check=check_ID();

if(identity_check==0){

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

}
}
//Function 1
// for this hardware parts would be needed.
// First we check for RPAS identifier, if any device has been changed or not
// there is a true value for drone id and a value calculated every time powering on the system
// if they both matched,then RPAS has not been tampered
// if they didnt , RPAS has been tampered , not allow to fly, A file is created and send to
//management client (signed by key/pair_C), then from management client to Management server.
//file informaton: 1)timestamp 2)Component's device id that didnt match.


/*
static int RPAS_identifier_check=0;

if (RPAS_identifier_check==0){
int identity_check= RPAS_identifier();

if(identity_check==0){

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

}
}*/
/*


//Function 2
// Drone_Id file is generated at each power cycle (if its not there in the directory)
// This is fetched by the MC at each power on cycle, therefore it is important to make sure
// the file exist
//DroneID.txt file generator function

static int drone_id_generation=0;

if(drone_id_generation==0)
{
//This has to be a function (taking in account all the unique things of the RPAS)

FILE *droneid;
droneid=fopen("/fs/microsd/log/DroneID.txt","r");
if(droneid==NULL){
  // the file is not there in the system, it has to be created
  // currently without signing, later has to be with signing support(Key_pair_C)

   call_DroneIDcreation( );
drone_id_generation=1;
}else{
	//file is there in the system check if it was witten by RPAS or not.
	// by checking for the signature
	char DroneID_container[30];
	char TAG_DRONE_ID[10]="DroneID";
	char DRONE_ID_FILE[100]="/fs/microsd/log/DroneID.txt";
	int check_validity_DroneID=call_file_read(DRONE_ID_FILE,TAG_DRONE_ID,DroneID_container,0);
	// if non tampered -->then continue
	if (check_validity_DroneID==0){
		remove("/fs/microsd/log/DroneID.txt");
		call_DroneIDcreation( );
	}else{
		printf("\n\nDRONEID.txt already made\n\n");
	}
	// if tampered --> then delete it and create a genuine one. with a valid signature
drone_id_generation=1;
}

}

// Function 3
// Till now RPAS identifier and Drone id is set
// Now its important to check that if the Drone has been registered or not
// Algo for that:
// Check for the UUID.txt file  : this file would be signed using MS key (Key/pair C)
// and send by the MC at each power on. Thats how the integrity of  our RPAS and MC would be
// maintained. In other way, we can say that our RPAS wont get ready to fly until our MC is not in use.
//
// It will contain following things: 1) Drone_ID 2)UUID 3) Signature(key/pair_C)
// This function will check for the Drone_Id in the file and will check its validity using
// the Key/pair_C
// if Drone_ID comes out to be same as the one mentioned in the above two function and file
// is non tampered(checked using Key/pair_C) then we can say that Registration has been done

static int UUID_check=0;
if (UUID_check==0)
{
	FILE *fptr_registration;
	fptr_registration=fopen("/fs/microsd/log/UUID.txt","r");
	printf("\ninside UUID check\n");
	if(fptr_registration!=NULL)
	{//file is present
		fclose(fptr_registration);
		char DroneID_container_1[30];
		char TAG_DRONE_ID_1[10]="DroneID";
		char DRONE_ID_FILE_1[100]="/fs/microsd/log/UUID.txt";
		int check_validity_DroneID_1=call_file_read(DRONE_ID_FILE_1,TAG_DRONE_ID_1,DroneID_container_1,1);
		// if non tampered -->then continue
		if (check_validity_DroneID_1==0){
			//UUID file is not valid, remove the invalid file
			remove("/fs/microsd/log/UUID.txt");
			mavlink_log_critical(mavlink_log_pub, "Drone registration not done, Connect Management Client and Internet to begin registration");

		}else{
			//file is valid, check for strings
			if(strcmp(DroneID_container_1,"ABBDDJEDNDJK")==0){
				UUID_check++;//Drone is registered
				printf("\n\nDrone is registered\n\n");
			}else{
				// UUID.txt file is for different drone
				// not UUID.txt of this drone
				remove("/fs/microsd/log/UUID.txt");
				mavlink_log_critical(mavlink_log_pub, "Drone registration not done, Connect Management Client and Internet to begin registration");
			}

		}



	}
	else{
	// registration has not been done, also recommend to connect to management client and
	// internet
	// UUID.txt file not present

		mavlink_log_critical(mavlink_log_pub, "Drone registration not done, Connect Management Client and Internet to begin registration");
		return false;
	}
}




// 3) Key rotation start, if needed

static int checky=0;
if(checky==0)
{
	FILE *fptr_key;
	fptr_key=fopen("/fs/microsd/log/KeyRotation.txt","r");

	if(fptr_key!=NULL)
	{// if present, then check for its validity


		fclose(fptr_key);
		char DroneID_container_2[30];
		char FILE_ID_container[30];
		char TAG_DRONE_ID_2[10]="DroneID";
		char TAG_KEY_ROT_ID[20]="FILE_ID";
		char DRONE_ID_FILE_2[100]="/fs/microsd/log/KeyRotation.txt";
		int check_validity_DroneID_2=call_file_read(DRONE_ID_FILE_2,TAG_DRONE_ID_2,DroneID_container_2,1);
		// if non tampered -->then continue
		check_validity_DroneID_2=call_file_read(DRONE_ID_FILE_2,TAG_KEY_ROT_ID,FILE_ID_container,1);

		if (check_validity_DroneID_2==0){
			//KeyRotation.txt file is not valid, remove the invalid file
			remove("/fs/microsd/log/KeyRotation.txt");


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
					remove("/fs/microsd/log/KeyRotation.txt");
					checky=1;
				}else{// The file for key rotation has already been used
				      // key rotation will not take place
					remove("/fs/microsd/log/KeyRotation.txt");

				}

			}else{
				// KeyRotation.txt file is for different drone
				//
				remove("/fs/microsd/log/KeyRotation.txt");
			}

		}

	}else{
	//no need of key rotation as there is no keyRotation.txt file present.
		//check for KeyChangePerm.txt
		FILE *fptr_keychange;
		char fname_keychange[40]="/fs/microsd/log/KeyChangePerm.txt";

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
			char DRONE_ID_FILE_2[60]="/fs/microsd/log/KeyChangePerm.txt";
			int check_validity_DroneID_3=call_file_read(DRONE_ID_FILE_2,TAG_DRONE_ID_2,DroneID_container_3,1);
		// if non tampered -->then continue
			check_validity_DroneID_3=call_file_read(DRONE_ID_FILE_2,TAG_KEY_ID,KEY_ID_container,1);

			check_validity_DroneID_3=call_file_read(DRONE_ID_FILE_2,MODULUS_tag,MODULUS_FILE,1);
			if(check_validity_DroneID_3==1){
				//file is valid
				//Now check for DroneID,KEY_ID and MODULUS
				//taking above from PublicKeyNew.txt
				char PKN_Modulus[720];
				char PKN_DroneID[30];
				char PKN_KEY_ID[30];
				char PKN_file[60]="/fs/microsd/log/PublicKeyNew.txt";
				check_validity_DroneID_3=call_file_read(PKN_file,MODULUS_tag,PKN_Modulus,0);
				check_validity_DroneID_3=call_file_read(PKN_file,TAG_DRONE_ID_2,PKN_DroneID,0);
				check_validity_DroneID_3=call_file_read(PKN_file,TAG_KEY_ID,PKN_KEY_ID,0);
				if(check_validity_DroneID_3==0){
					//not valid PublicKeyNew.txt, remove it
					remove("/fs/microsd/log/PublicKeyNew.txt");
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
						checky++;

					}else{	//
						//file is valid but does not belong to this drone.
						remove("/fs/microsd/log/KeyChangePerm.txt");
					}
				}

			}else{
				//file is not valid, remove the invalid file
				remove("/fs/microsd/log/KeyChangePerm.txt");
			}

		}else{//file is not present, no need to change keys

		}

	}
}




/// Some parameters are very crucial and therefore acess for the same is not given
/// to the user to write over the same Parameter update and set
static int para_set=0;
if(para_set==0){
	FILE *fptr_ParamChange;
	fptr_ParamChange=fopen("/fs/microsd/log/ParamChangePerm.txt","r");
	if(fptr_ParamChange!=NULL)
	{	fclose(fptr_ParamChange);
		//a modified ParamInuse.txt is needed to be formed
		//also check the presence of ParamInuse.txt
		FILE *fptr_ParamInuse;
		fptr_ParamInuse=fopen("/fs/microsd/log/ParamInuse.txt","r");
		if(fptr_ParamInuse!=NULL){
			fclose(fptr_ParamInuse);
			//both the files are present
			// validate both of them
			//ParamInuse.txt has to be validated with RFM public key
			//ParamChangePerm.txt has to be validated with FMP public key
			char DroneID_container_4[30];
			char DroneID_container_5[30];
			char TAG_DRONE_ID_4[10]="DroneID";
			char PARAM_IN_FILE_2[30]="/fs/microsd/log/ParamInuse.txt";
			char PARAM_CHANGE_FILE_2[100]="/fs/microsd/log/ParamChangePerm.txt";
			int check_validity_DroneID_3=call_file_read(PARAM_IN_FILE_2,TAG_DRONE_ID_4,DroneID_container_4,0);
			// if non tampered -->then continue
			if(check_validity_DroneID_3==1){
				if(strcmp(DroneID_container_4,"ABBDDJEDNDJK")==0){
					//droneID is valid
				}else{
					//droneID is not of this drone, remove the file
					remove("/fs/microsd/log/ParamInuse.txt");
					mavlink_log_critical(mavlink_log_pub, "Authentic Parameter file missing");

					return false;
				}
			}else{
				//remove the file
				remove("/fs/microsd/log/ParamInuse.txt");
				mavlink_log_critical(mavlink_log_pub, "Authentic Parameter file missing");

				return false;
			}
			int check_validity_DroneID_4=call_file_read(PARAM_CHANGE_FILE_2,TAG_DRONE_ID_4,DroneID_container_5,1);
			if(check_validity_DroneID_4==1){
				if(strcmp(DroneID_container_5,"ABBDDJEDNDJK")==0){
					//droneID is valid
					// At this point both the files have been validated ()
					//Begin the modification of ParamInuse.txt
					ParamInuseModify();
					//Once the new ParamInuse.txt file is formed,set the params
					//as specified by the file.
					ParamSetfile();
					para_set=1;
				}else{
					//droneID is not valid, remove the file
					remove("/fs/microsd/log/ParamChangePerm.txt");
					goto ParamInuse_way;
				}
			}else{
				//remove the file
				remove("/fs/microsd/log/ParamChangePerm.txt");
				goto ParamInuse_way;
			}

		}else{
			//Absence of ParamInuse.txt
			//this is a serious issue, firmware update is required or
			// a function needs to be made to generate ParamInuse.txt from initial
			// important parameters
			mavlink_log_critical(mavlink_log_pub, "Authentic Parameter file missing");
			return false;

		}
	}else
	{ParamInuse_way:
		//use the already made ParamInuse.txt in the RFM to set parameters
		//check the presence of ParamInuse.txt
		FILE *fptr_ParamInuse;
		fptr_ParamInuse=fopen("/fs/microsd/log/ParamInuse.txt","r");
		if(fptr_ParamInuse!=NULL){
			// validate ParamInuse.txt
			//and then set the param and values accordingly
			char DroneID_container_4[30];
			char TAG_DRONE_ID_4[10]="DroneID";
			char PARAM_IN_FILE[30]="/fs/microsd/log/ParamInuse.txt";
			int check_validity_DroneID_3=call_file_read(PARAM_IN_FILE,TAG_DRONE_ID_4,DroneID_container_4,0);
			// if non tampered -->then continue
			if(check_validity_DroneID_3==1){
				if(strcmp(DroneID_container_4,"ABBDDJEDNDJK")==0){
					//droneID is valid
					//Begin reading the file and start setting the parameters
					ParamSetfile();
					para_set=1;


				}else{
					//droneID is not valid, remove the file
					//ParamInuse.txt is for another drone
					remove("/fs/microsd/log/ParamInuse.txt");
					mavlink_log_critical(mavlink_log_pub, "Authentic Parameter file missing");


				}
			}else{
				//ParamInuse.txt file is not valid, you need to remove it
				remove("/fs/microsd/log/ParamInuse.txt");
				mavlink_log_critical(mavlink_log_pub, "Authentic Parameter file missing");

			}


		}else{   // this file has to be there inside the RFM, otherwise Drone wont fly.
			// this is an issue which requires to be logged into the SystemLog.txt
			//Drone will not get armed to fly before it does not receive the Param
			mavlink_log_critical(mavlink_log_pub, "Authentic Parameter file missing");

			return false;
		}

	}



}






//function 4 : recentPA.txt for flight log management.
//check for recentPA.txt:::: this is for bundling purpose, do we need to start bundling or not
//
static int check_recent=0;
static int recentPA_presence=0;
static char PA_ID[100];
if(check_recent==0){
	FILE *fptr_recent;

	fptr_recent=fopen("./log/recentPA.txt","r");
	if(fptr_recent!=NULL){

		//if file is present then there might be the possiblity upcoming PA to be the old one.

		fclose(fptr_recent);
		int action=check_recentPA(PA_ID);
		//checks to be done :
		//1)if the current time lies inside the time period of start time and end time
		// this will initiate the flight logs bundling and then the drone will wait for fetched.txt file
		//updated through a parameter in recentPA.txt

		if(action==0 || action==1){
			// start bundling and update recentPA.txt with fetchrequired=0;
		//	Bundling_begins();//how to bundle part/modifying and updating recentPA.txt
			char yes[20]="1";
			update_recentPA(0,yes);
			recentPA_presence++;

			check_recent=1;

		}else if(action==3){
			mavlink_log_critical(mavlink_log_pub, "Invalid pa handling file");
			remove("./log/recentPA.txt");
			return false;


		}else{//no need to start bundling 2
			check_recent=1;
			//make sure upcoming PA has same pa_id. as in recentPA.txt
			recentPA_presence++;

		}

	}else{
		//file not present
		// upcoming permission artefact is new, a fresh recentPA.txt file would be made.
		check_recent=1;


	}
}
*/

/// to check for permission artifact (is it present here or not) and its validation
/// then verifying for geo coordinates and time period
///
//Now first read recentPA.txt to know if  fetch of fetched.txt is required or not
/*
int check_status_fetch=0;
int status_fetch=0;
if(recentPA_presence){
	// to know if we are supposed to look for fetched.txt
	status_fetch=read_for_fetch();

	if(status_fetch!=1){
		//fetch is not required (no need to check for fetched.txt file_)
	}else{
		// fetch is required (need to check for fetched.txt file)
		check_status_fetch=check_fetch();
		//if check_status==0 : no file present, return false
		//               ==1  : file is present but not valid return false(attempt to hack)
		//               ==2  : file is present and valid, delete recentPA.txt and set recentPA_presence=0
	}
}
*/

// Tackling first thing : 1) Check for valid PA
//static int checky0=0;

static int checky1=0;
static int tampercheck=0;
//int time_verification=0;
//int drone_id_verification;
int verification=0;


if(!checky1)
{

	if((time_since_boot > 10_s))
	{   
		
		FILE *file;
		file = fopen("/fs/microsd/log/permission_artifact_breach.xml", "r");// ./log/ for posix /log/ for nuttx
		if (file){

			fclose(file);
			if(tampercheck==0)
			{	
				verification=Is_PA_VAlid();// if pa
				//char file_name[]="/fs/microsd/log/permission_artifact_breach.xml";
				//char Tag_Digest_value0[12]="DigestValue";
				//char *Digest_Value0=(char*) malloc(100*sizeof(char));
			/*	char *Reference_canonilized=(char*) malloc(5000*sizeof(char));
				char *xmlExeclusive=(char*) malloc(5000*sizeof(char));
  				char *output=(char*) malloc(5000*sizeof(char));
				char *Sha_of_Reference=(char*) malloc(300*sizeof(char));
				Reference_canon(file_name,Reference_canonilized);
                                cleanerXML(Reference_canonilized,output);
                                xmlExeclusiveCanon(output,xmlExeclusive);
				free(Reference_canonilized);
                                free(output);
                                Sha256_implement(xmlExeclusive,Sha_of_Reference);
				free(Sha_of_Reference);
				free(xmlExeclusive);
				char *SignedInfo_canonilized=(char*) malloc(5000*sizeof(char));
                                char *outputS=(char*) malloc(5000*sizeof(char));
                                char *Sha_of_SignedInfo=(char*) malloc(5000*sizeof(char));
				SignedInfo_canon(file_name,SignedInfo_canonilized);
        			cleanerXML(SignedInfo_canonilized,outputS);
        			//printf("\nSigned Info :%s\n",SignedInfo_canonilized);
				//free(outputS);
        			free(SignedInfo_canonilized);
				Sha256_implement(outputS,Sha_of_SignedInfo);
				free(Sha_of_SignedInfo);

                                char Tag_Digest_value0[12]="DigestValue";
				char *Digest_Value0=(char*) malloc(100*sizeof(char));
				char *Digest_Value_in_hex=(char*) malloc(300*sizeof(char));
				getTagvalue(Tag_Digest_value0,Digest_Value0,file_name);

				base64decoder(Digest_Value0, Digest_Value_in_hex);

				free(Digest_Value0);
				free(Digest_Value_in_hex);*/

				// storing Signed value
				/* char Tag_Signed_Value0[17]="SignatureValue";

				//  char Signature_Value0[550];// in base 64
				char *Signature_Value0=(char*) malloc(550*sizeof(char));
				char Signature_Value_in_hex[500];//;=(char*) malloc(500*sizeof(char));
				getTagvalue(Tag_Signed_Value0,Signature_Value0, file_name);

				base64decoder(Signature_Value0, Signature_Value_in_hex);

				free(Signature_Value0);
				//free( Signature_Value_in_hex);

  int cb;
   mp_int message,modulus,public_key,Decrypted;
   cb= mp_init_multi(&message,&modulus,&public_key,&Decrypted,NULL);
   cb= mp_read_radix(&message,Signature_Value_in_hex,16);
  // free(Signature_Value_in_hex);

   /// Public key of dgca: This has to be taken from a reserved file in directory./ firmware update

 //  mp_int ;
    char Modulus[513]="d6912a773335d3a193c742071762794e26bcac49d78b6b65a784e1c18d18f2f88b9f13d7fc41a5b9e11e3d75905b0f8373dbac5658962940d104711467814b7319774014644df82e9764b4e852cb7547d56de3224aada000db231b1356ffe86391b3eca2ddf67d44f40ea6e84092cc67387db3a2042487b8dacfe5b588738973a59f5fa813a33e4bb8fbafa407794db990a307201b0a0fbd92ddd181a868a6cfa64625353714e9b1de6f81f3addbf5b202030dbd9e51385db9314591f01af01f07cc92f18ebe60545cea53eb54438e251fd88e7380b5a0612c29ccb7dfd88f7a7d7f0dd07a9b596923602798ad9f1bbb89fffd3cdb8fb94c48640d27fa681e47";
  // char Modulus[513]="ab9d5c8d1fe67207749d63b7dcedd233ce32bb70d175a1bc38c612ab33e2c58e51f83f2788e4d52d9bceb5a1513929de3f526650071a067e6c161b05c60a495fc3ba79ed26f4fa8b2fe2ca8dec44b39759f39206f06a85f9424005a29f05e4cf3a0239340c28c993c1a61cf1b2b6b57c7d8e576ae86827f812b327625baec9ecbf55f1651d35600b9f955f6c2f3bea3aa5852ecdd36a0af818c19acc1030979bed3c89993faa92e0aa0502413b3ca86bbf63477f12ac069aff7137cb72c57f886da79033bbb3b4df0f6cc7fcc18e343aa76036681a566311e267c03b65c98abc91e58f090020c67f776199c0eb76d7e6363687475d3da36ff050f85275607fdd";
  cb= mp_read_radix(&modulus,Modulus,16);

  // mp_int ;
  cb= mp_read_radix(&public_key,"65537",10);

   //mp_int ;
  //cb= mp_exptmod(&message, &public_key,&modulus,&Decrypted);   //                       this part for decrypting encrypted text
  if(cb!=0){
          mavlink_log_critical(mavlink_log_pub, "Problem in bignum ");
  }*/
  // char message_string_hex[513];
  // mp_to_hex(&Decrypted,message_string_hex,sizeof(message_string_hex));
  // free(message_string_hex);
  // mp_clear_multi(&message,&modulus,&public_key,&Decrypted,NULL);



				//getTagvalue(Tag_Digest_value0,Digest_Value0,file_name);
				//free(Digest_Value0);
				//char Reference_canonilized[2000]="patlee";//=(char*) malloc(15000*sizeof(char));
				//printf("Refer %s",Reference_canonilized);

				//Reference_canon(file_name,Reference_canonilized);
				//free(Reference_canonilized);
				if (verification==1)
				{
					//printf("bale bale\n");
					mavlink_log_info(mavlink_log_pub," Permission Artefact is non tampered (valid)");
					tampercheck=1;
					checky1=1;
					//fetching_publish_padata();
					//data_fetch_delete();
					//checky0=1;
				}else if(verification==2)
				{
					mavlink_log_critical(mavlink_log_pub, "Permission Artefact  has bad Sign");
					//remove("./log/permission_artifact_breach.xml");
					return false;

				}
				else{
					//
					mavlink_log_critical(mavlink_log_pub, "Permission Artefact  has been tampered(non valid)");
					//remove("./log/permission_artifact_breach.xml");
					return false;
				}
			}
		// 2) Now extracting date_time and geo_coordinates and checking for time and geofence breach.


/*

			time_verification=date_time_extract_and_check();
			printf("\ntime verification ::::%d\n",time_verification);
			if (time_verification==0){
				printf("\nNot in time limit\n");
				mavlink_log_critical(mavlink_log_pub, "Preflight Fail: Not allowed to fly, Time breach (bad time)");
				return false;
			}else if(time_verification==2)
			{

				printf("\nNot in correct place \n");
				mavlink_log_critical(mavlink_log_pub, "Preflight Fail: Not allowed to fly, Geo breach (bad geo)");
				return false;
			}else{
				printf("\n In time limit and permitted area. \n");
			}

			//if(recentPA_presence==0){
			drone_id_verification=DroneIDverification();
			//}else{
			//	drone_id_verification=DroneIDverification(PA_ID);
			//}
			if(drone_id_verification==0){
				printf("\nNot correct DroneID \n");
				mavlink_log_critical(mavlink_log_pub, "Preflight Fail: Not valid DroneID");
				return false;
			}else if(drone_id_verification==2){
				mavlink_log_critical(mavlink_log_pub, "Please upload the previous permission artefact in use, its end limit hasnt reached yet");
				return false;
			}
			else{
			//	if(recentPA_presence==0){

					// data_fetch_delete will only run for new PAs(when recentPA.txt is not present)
					// Now delete PA.xml and take its data inside recentPA.txt

					data_fetch_delete();// at this point pa_data.msg will get published
				// }else{
					// this suggests that recentPA.txt has the information of time and area limit
					// from recentPA.txt it would be published to pa_data.msg
					//fetching_publish_padata();

				//}
			checky0=1;
			checky1=1;
			}
	*/



		}else{
			mavlink_log_critical(mavlink_log_pub, "Preflight Fail: no Permission artefact found. Please send Permission artefact from Management client.");
			return false;
		}

	}else{
		printf("\nwaiting for gps signals....\n");
		return false;
	}

}


/*static int verification=1;
if(verification==1){
data_fetch_delete();
verification=0;
}*/
//printf("%d",verification);
//int in_fence = 0;
//param_get(param_find("IN_FENCE"),&in_fence);
//if (!in_fence) return false;
/*
char fi[40]="46790";
static int chadd=0;
if(chadd<=1){
	Key_rotation_start(fi);
	chadd++;
}*/


////
	report_failures = (report_failures && status_flags.condition_system_hotplug_timeout
			   && !status_flags.condition_calibration_enabled);

	bool failed = false;

	failed = failed || !airframeCheck(mavlink_log_pub, status);
	failed = failed || !sdcardCheck(mavlink_log_pub, status_flags.sd_card_detected_once, report_failures);

	/* ---- MAG ---- */
	{
		int32_t sys_has_mag = 1;
		param_get(param_find("SYS_HAS_MAG"), &sys_has_mag);

		if (sys_has_mag == 1) {

			/* check all sensors individually, but fail only for mandatory ones */
			for (unsigned i = 0; i < max_optional_mag_count; i++) {
				const bool required = (i < max_mandatory_mag_count) && (sys_has_mag == 1);
				bool report_fail = report_failures;

				int32_t device_id = -1;

				if (!magnetometerCheck(mavlink_log_pub, status, i, !required, device_id, report_fail)) {
					if (required) {
						failed = true;
					}

					report_fail = false; // only report the first failure
				}
			}

			// TODO: highest priority mag

			/* mag consistency checks (need to be performed after the individual checks) */
			if (!magConsistencyCheck(mavlink_log_pub, status, report_failures)) {
				failed = true;
			}
		}
	}

	/* ---- ACCEL ---- */
	{
		/* check all sensors individually, but fail only for mandatory ones */
		for (unsigned i = 0; i < max_optional_accel_count; i++) {
			const bool required = (i < max_mandatory_accel_count);
			bool report_fail = report_failures;

			int32_t device_id = -1;

			if (!accelerometerCheck(mavlink_log_pub, status, i, !required, device_id, report_fail)) {
				if (required) {
					failed = true;
				}

				report_fail = false; // only report the first failure
			}
		}

		// TODO: highest priority (from params)
	}

	/* ---- GYRO ---- */
	{
		/* check all sensors individually, but fail only for mandatory ones */
		for (unsigned i = 0; i < max_optional_gyro_count; i++) {
			const bool required = (i < max_mandatory_gyro_count);
			bool report_fail = report_failures;

			int32_t device_id = -1;

			if (!gyroCheck(mavlink_log_pub, status, i, !required, device_id, report_fail)) {
				if (required) {
					failed = true;
				}

				report_fail = false; // only report the first failure
			}
		}

		// TODO: highest priority (from params)
	}

	/* ---- BARO ---- */
	{
		int32_t sys_has_baro = 1;
		param_get(param_find("SYS_HAS_BARO"), &sys_has_baro);

		bool baro_fail_reported = false;

		/* check all sensors, but fail only for mandatory ones */
		for (unsigned i = 0; i < max_optional_baro_count; i++) {
			const bool required = (i < max_mandatory_baro_count) && (sys_has_baro == 1);
			bool report_fail = (required && report_failures && !baro_fail_reported);

			int32_t device_id = -1;

			if (!baroCheck(mavlink_log_pub, status, i, !required, device_id, report_fail)) {
				if (required) {
					baro_fail_reported = true;
				}

				report_fail = false; // only report the first failure
			}
		}
	}

	/* ---- IMU CONSISTENCY ---- */
	// To be performed after the individual sensor checks have completed
	{
		if (!imuConsistencyCheck(mavlink_log_pub, status, report_failures)) {
			failed = true;
		}
	}

	/* ---- AIRSPEED ---- */
	/* Perform airspeed check only if circuit breaker is not engaged and it's not a rotary wing */
	if (!status_flags.circuit_breaker_engaged_airspd_check &&
	    (status.vehicle_type == vehicle_status_s::VEHICLE_TYPE_FIXED_WING || status.is_vtol)) {

		int32_t airspeed_mode = 0;
		param_get(param_find("FW_ARSP_MODE"), &airspeed_mode);
		const bool optional = (airspeed_mode == 1);

		int32_t max_airspeed_check_en = 0;
		param_get(param_find("COM_ARM_ARSP_EN"), &max_airspeed_check_en);

		float airspeed_trim = 10.0f;
		param_get(param_find("FW_AIRSPD_TRIM"), &airspeed_trim);

		const float arming_max_airspeed_allowed = airspeed_trim / 2.0f; // set to half of trim airspeed

		if (!airspeedCheck(mavlink_log_pub, status, optional, report_failures, prearm, (bool)max_airspeed_check_en,
				   arming_max_airspeed_allowed)
		    && !(bool)optional) {
			failed = true;
		}
	}

	/* ---- RC CALIBRATION ---- */
	if (status.rc_input_mode == vehicle_status_s::RC_IN_MODE_DEFAULT) {
		if (rcCalibrationCheck(mavlink_log_pub, report_failures, status.is_vtol) != OK) {
			if (report_failures) {
				mavlink_log_critical(mavlink_log_pub, "RC calibration check failed");
			}

			failed = true;

			set_health_flags(subsystem_info_s::SUBSYSTEM_TYPE_RCRECEIVER, status_flags.rc_signal_found_once, true, false, status);
			status_flags.rc_calibration_valid = false;

		} else {
			// The calibration is fine, but only set the overall health state to true if the signal is not currently lost
			status_flags.rc_calibration_valid = true;
			set_health_flags(subsystem_info_s::SUBSYSTEM_TYPE_RCRECEIVER, status_flags.rc_signal_found_once, true,
					 !status.rc_signal_lost, status);
		}
	}

	/* ---- SYSTEM POWER ---- */
	if (status_flags.condition_power_input_valid && !status_flags.circuit_breaker_engaged_power_check) {
		if (!powerCheck(mavlink_log_pub, status, report_failures, prearm)) {
			failed = true;
		}
	}

	/* ---- Navigation EKF ---- */
	// only check EKF2 data if EKF2 is selected as the estimator and GNSS checking is enabled
	int32_t estimator_type = -1;

	if (status.vehicle_type == vehicle_status_s::VEHICLE_TYPE_ROTARY_WING && !status.is_vtol) {
		param_get(param_find("SYS_MC_EST_GROUP"), &estimator_type);

	} else {
		// EKF2 is currently the only supported option for FW & VTOL
		estimator_type = 2;
	}

	if (estimator_type == 2) {

		const bool ekf_healthy = ekf2Check(mavlink_log_pub, status, false, report_failures) &&
					 ekf2CheckSensorBias(mavlink_log_pub, report_failures);

		// For the first 10 seconds the ekf2 can be unhealthy, and we just mark it
		// as not present.
		// After that or if report_failures is true, we'll set the flags as is.

		if (!ekf_healthy && time_since_boot < 10_s && !report_failures) {
			set_health_flags(subsystem_info_s::SUBSYSTEM_TYPE_AHRS, true, false, false, status);

		} else {
			set_health_flags(subsystem_info_s::SUBSYSTEM_TYPE_AHRS, true, true, ekf_healthy, status);
		}

		failed |= !ekf_healthy;
	}

	/* ---- Failure Detector ---- */
	if (!failureDetectorCheck(mavlink_log_pub, status, report_failures, prearm)) {
		failed = true;
	}

	failed = failed || !manualControlCheck(mavlink_log_pub, report_failures);
	failed = failed || !cpuResourceCheck(mavlink_log_pub, report_failures);

	/* Report status */
	return !failed;
}
