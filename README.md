Augmenting PX4 software with PKI(Public Key Infrastructure) architecture:
Details lies in commits.
implemented following apis in px4

RFM API Reference
1. Get_rfmInfo
Parameter 1: IST date-time (fetched from GLPS data or any other source)
Output 1: RFM public key
Output 2: Digital Sky Public Key being used in RFM
Output 3: Firmware Version
Output 4: RFM version
Output 5: RPAS category
Output 6: Operator ID obtained during registration
[Provide RPAS info]
a. The API will return the required information for any external application.
b. The RFM Provider has to sign the information using RFM private key.


2. Apply Permission Artefact:
Parameter 1: Permission_artefact
Parameter 2: IST date-time (fetched from GLPS data or any other source)
[ provide permission artefact to the RFM ]
a. This API should be called by Flight Controller or host system every
time the RPAS is powered on.
b. RFM will verify the signed permission artefact using the Digital Sky
Public Key.
c. RFM, being stateless, will not store the permission artefacts in the
memory, so Flight Module Provider is responsible for providing relevant
permission artefact to the RFM on every power cycle. The permission
state is maintained by the RFM for the current power session only.
d. RFM will use timestamp provided by Flight Module Provider, the
artefact publish timestamp (marked by DSP when the artefact is
created), and the ttl (time to live) value to validate that the artefact is
eligible at that time and that the system time-date is not falsified.

3. Get_Geofence_restriction:
Parameter 1: IST date-time
[ provides geofence information to the Flight Controller or Host system after a
valid permission artefact has been applied ]
a. RFM will verify the date-time stamp to reassert the time validity of the
permission artefact
b. This API is provided by RFM for the Flight Controller or Host system.
The Flight Controller is expected to use this data for enforcing safety
restrictions or warning the pilot.
c. It is responsibility of Flight Module Provider to ensure that the RPAS
follows the geofence restrictions. In case of geofence or time limit
violations, the Flight Module Provider should provide the event details
to the RFM.

4. Get_time_restriction:
Parameter 1: IST date-time
[Provides allowed time limit after a valid permission artefact has been applied]
a. RFM will verify the date-time stamp to reassert the time validity of the
permission artefact.
b. It is responsibility of Flight Module Provider to ensure that the RPAS is
landed before the permitted time period expires. In case of time limit
violations, the Flight Module Provider should provide the event details
to the RFM

5. Log_takeoff_location:
Parameter1: IST date-time
Parameter2: GLPS coordinates
[To provide takeoff event information to library for internal logging.]
a. The Flight Module Provider is responsible to implement the
functionality to ensure that this API is called immediately after takeoff
event.
b. The Flight Module Provider should not provide any means to bypass
the notification to RFM on a takeoff event.
c. The RFM will store and package all such events in one flight log
against a permission artefact.

6. Log_Land_location:
Parameter1: IST date-time
Parameter2: GLPS coordinates
[To provide land event information to RFM for internal logging.]
a. The Flight Module Provider is responsible to implement the
functionality to ensure that this API is called immediately after land
event.
b. The Flight Module Provider should not provide any means to bypass
the notification to RFM on a land event.
c. The RFM will store and package all such events in one flight log
against a permission artefact.


7. Log_geofence_breach:
Parameter1: IST date-time
Parameter2: Breach_start_timestamp
Parameter3: Breach_stop_timestamp
Parameter2: GLPS coordinates
[To provide geofence breach event information to RFM for internal logging.]
a. The Flight Module Provider is responsible to implement the
functionality to ensure that this API is called immediately after a
geofence breach incidence.
b. The Flight Module Provider is responsible for implementing
functionality to call this API at 1 Hz, from the time when Geofence is
breached to when the drone lands.
c. The RFM will store and package all such events in one flight log
against a permission artefact.


8.Log_timelimit_breach:
Parameter1: IST date-time
Parameter2: time_overrun_start_timestamp
Parameter3: time_overrun_land_timestamp
[To provide time limit breach event information to RFM for internal logging.]

REVISION 1a. The Flight Module Provider is responsible to implement the
functionality to ensure that this API is called immediately for a time limit
overrun event.
b. During the first call, when the drone is still in flight, the
time_overrun_timestamp parameter will not have any value. The event
would be logged as start of time overrun.
c. When the drone lands this API should be called again to indicate end
of flight. In this case the API caller needs to provide
time_overrun_start_timestamp again. This event will be registered as
end of time-overrun event.
d. The RFM will store and package all such events in flight log against a
permission artefact.



9. Get_individual_flight_logs:
Parameter1: IST date-time
[To get the individual flight logs from the RFM for storage]
a. The RFM must prepare a flight log after each land, geofence breach, or
flight time overrun.
b. These flight logs will be digitally signed by the RFM to make sure that
flight logs are not tampered with during the transport. It is Flight Module
Provider's responsibility to implement the functionality required to keep
the private key secure. These reports may be verified against that
public key that was shared during registration of the RPAS.
c. It is responsibility of the Flight Module Provider to implement the
functionality to collect these reports immediately after every takeoff,
land, Geofence and time-limit breach event from RFM. In case of crash
where such events were not registered, the RFM on next power cycle,
should close the log with failed landing incidence.
d. The Flight Module Provider is responsible to implement the
functionality to store all the flight logs until they are uploaded to the
DSP. This storage should provide the ‘write access’ only to the
applications authorized by Flight Module Provider. The read access to
the storage should be available in case of accidents. The Flight Module
Provider is required to provide the authorities with the specialized
equipment required to read the flight logs from a crashed/damaged
RPAS’s onboard storage.
e. Flight Module Provider is responsible to provide communication
interface to upload the flight log directly to the DSP or through external
applications, such as Ground Control Station.
f. The Flight Module Provider is responsible to ensure the storage space
for logs. If storage space is not available for logging then the flight
should not be allowed.


10. Bundle_flight_logs:
Parameter1: IST time-stamp
Parameter2: List_of_individual_incidence_reports_from_storage
Parameter3: Permission Artefact
Output1: Bundled_incident_report_with_digital signature.

FIRST EDITION[To bundle the signed flight logs from storage into a single bundle (a bundle
per permission artefact).]
a. After the time-period of a permission artefact is over, the user has to
submit the flight logs to Digital Sky APIs within 3 days. The Flight
Module Provider is required to implement the functionality to provide all
the flight logs associated with a particular permission artefact and pass
them on to this API to get in return a signed bundle.
b. Flight Module Provider is responsible to provide communication
interface to upload this bundle directly to the DSP or through external
applications, such as Ground Control Station.
