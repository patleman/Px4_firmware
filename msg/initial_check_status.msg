uint64 timestamp                 # time since system start (microseconds)



# RPAS_identifier status possible values = 0,1,2
# 0-> Unauthorized parts attached
# 1-> Authorized parts attached
# 2-> file not find, serious problem 

uint32 identifier_status 

# Registration Check presence of UUID file signed from FMP key.

uint32 registration_status # 0,1,2,3 



# fetch required or not , wait for MC to fetch the flight logs. this will show message "Connect MC to fetch logs"
bool fetch_required 

#Permission Artefact status
bool pa_present # if pa is present or not

bool bad_sign  #if signature is not of the form "1ffffff...."

# this is the drone_id check done in pa_valid and checked in preflight_check

bool drone_id # drone id doesnt match to what is write inside the firmware

uint32 validate  # if the PA has come from DGCA and is non tampered




# this is the condition that to be checked in pa_valid (only if recentPA_present set to 1 in log_management) and retrieved in pref
#light check

bool same_pa

# this is set in log_management and retrieved in pa_valid
bool recent_pa_present

