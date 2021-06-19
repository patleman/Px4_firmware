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
#include "MP_INT.hpp"
#include "PA_EXTRACT.hpp"
#include "Identifier.hpp"

#include <drivers/drv_hrt.h>
#include <HealthFlags.h>

#include <px4_defines.h>
#include <lib/parameters/param.h>
#include <systemlib/mavlink_log.h>
#include <uORB/Subscription.hpp>

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
//Function 1
// for this hardware parts would be needed.
// First we check for RPAS identifier, if any device has been changed or not
// there is a true value for drone id and a value calculated every time powering on the system
// if they both matched,then RPAS has not been tampered
// if they didnt , RPAS has been tampered , not allow to fly, A file is created and send to
//management client (signed by key/pair_C), then from management client to Management server.
//file informaton: 1)timestamp 2)Component's device id that didnt match.

int identity_check= RPAS_identifier();

if(identity_check==0){
	printf("Hardware parts have been changed");
	return false;
}



//Function 2
// Drone_Id file is generated at each power cycle (if its not there in the directory)
// This is fetched by the MC at each power on cycle, therefore it is important to make sure
// the file exist
//DroneID.txt file generator function

static int drone_id_generation=0;
if(drone_id_generation==0){
//This has to be a function (taking in account all the unique things of the RPAS)
char Drone_ID[40]="11dd44@@334499200";// truth value
FILE *droneid;
droneid=fopen("./log/Drone_ID.txt","r");
if(droneid==NULL){
  // the file is not there in the system, it has to be created
  // currently without signing, later has to be with signing support(Key_pair_C)
droneid=fopen("./log/Drone_ID.txt","w");
fprintf(droneid, "%s", Drone_ID);
fclose(droneid);
drone_id_generation=1;
}else{
	//file is there in the system check if it was witten by RPAS or not.
	// by checking for the signature
	// if non tampered -->then continue
	// if tampered --> then delete it and create a genuine one. with a valid signature
}

}

// Function 3
// Till now RPAS identifier and Drone id is set
// Now its important to check that if the Drone has been registered or not
// Algo for that:
// Check for the UUID.txt file  : this file would be signed using MS key (Key/pair C)
// and send by the MC at each power on. Thats how the integrity of the our RPAS and MC would be
// maintained. In other way, we can say that our RPAS wont get ready to fly until our MC is not in use.
//
// It will contain following things: 1) Drone_ID 2)UUID 3) Signature(key/pair_C)
// This function will check for the Drone_Id in the file and will check its validity using
// the Key/pair_C
// if Drone_ID comes out to be same as the one mentioned in the above two function and file
// is non tampered(checked using Key/pair_C) then we can say that Registration has been done
FILE *fptr_registration;
char c_aux;
fptr_registration=fopen("./log/UUID.txt","r");
if(fptr_registration!=NULL){
	while(1){
		c_aux=fgetc(fptr_registration);
		if (c_aux==EOF){
			break;
		}
		else{
			fputc(c_aux,fptr_registration);
		}
	}
	fclose(fptr_registration);


}else{
	// registration has not been done, also recommend to connect to management client and
	// internet
	mavlink_log_critical(mavlink_log_pub, "Drone registration not done");
}

/// first priority is to check for permission artifact (is it present here or not) and its validation
/// then verifying for geo coordinates and time period
/// and then need for key rotation (separate app (decision at hold))


// Tackling first thing : 1) Check for valid PA
FILE *file;
file = fopen("./log/permission_artifact_breach.xml", "r");// ./log/ for posix /log/ for nuttx
if (file){
fclose(file);
static int checky0=0;
int verification=0;
if(!checky0)
{
verification=Is_PA_VAlid();
checky0=1;
}
if (verification==1){
    printf("bale bale\n");
    mavlink_log_info(mavlink_log_pub," Permission Artefact is non tampered");
}
/// 2) Now extracting date_time and geo_coordinates and checking for time and geofence breach.
static int checky2=0;
int time_verification;
if(!checky2 && (time_since_boot>10))
{
time_verification=date_time_extract_and_check();
if (time_verification==0){
    printf("\nNot in time limit\n");
    mavlink_log_critical(mavlink_log_pub, "Preflight Fail: Not allowed to fly, Time breach");
    return false;
}else{
    printf("\nIn time limit\n");
}
        if(time_verification==1){
          checky2=1;}
}
}else{
	mavlink_log_critical(mavlink_log_pub, "Preflight Fail: no Permission artefact found. Please send Permission artefact from Management client.");
	return false;
}



//int in_fence = 0;
//param_get(param_find("IN_FENCE"),&in_fence);
//if (!in_fence) return false;

// 3) Key rotation start, if needed
int key_rot_done=0;
param_get(param_find("KEY_ROT"),&key_rot_done);
static int checky=0;
if(!checky && key_rot_done!=1)
{
Key_rotation_start();
checky=1;
}


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
