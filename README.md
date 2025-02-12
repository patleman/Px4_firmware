# **RFM API Reference**

## 1. **Get_rfmInfo**

### **Parameters**:
- **IST date-time**: Fetched from GLPS data or any other source.

### **Outputs**:
- **RFM public key**
- **Digital Sky Public Key** being used in RFM
- **Firmware Version**
- **RFM Version**
- **RPAS Category**
- **Operator ID** obtained during registration

### **Description**:
- Provides RPAS info.
- The API returns the required information for any external application.
- The **RFM Provider** must sign the information using the **RFM private key**.

---

## 2. **Apply Permission Artefact**

### **Parameters**:
- **Permission_artefact**
- **IST date-time**: Fetched from GLPS data or any other source.

### **Description**:
- This API should be called by **Flight Controller** or **Host System** every time the RPAS is powered on.
- **RFM** will verify the signed **permission artefact** using the **Digital Sky Public Key**.
- **RFM** is stateless and does not store the permission artefacts in memory, so the **Flight Module Provider** must provide the relevant permission artefact on every power cycle.
- The **permission state** is maintained by **RFM** only for the current power session.
- The **RFM** will validate the artefact using the **timestamp** provided by the **Flight Module Provider**, the **artefact publish timestamp** (marked by DSP when the artefact is created), and the **TTL (time to live)** value to ensure its validity.

---

## 3. **Get_Geofence_Restriction**

### **Parameters**:
- **IST date-time**

### **Description**:
- Provides **geofence** information to the **Flight Controller** or **Host System** after a valid **permission artefact** has been applied.
- **RFM** will verify the **date-time stamp** to reassert the time validity of the permission artefact.
- The **Flight Controller** uses this data to enforce **safety restrictions** or **warning the pilot**.
- The **Flight Module Provider** is responsible for ensuring that the RPAS follows **geofence restrictions**.
- In case of **geofence** or **time limit** violations, the **Flight Module Provider** must provide the event details to the **RFM**.

---

## 4. **Get_time_restriction**

### **Parameters**:
- **IST date-time**

### **Description**:
- Provides the allowed **time limit** after a valid **permission artefact** has been applied.
- **RFM** will verify the **date-time stamp** to reassert the time validity of the **permission artefact**.
- The **Flight Module Provider** is responsible for ensuring that the RPAS is landed before the **permitted time** period expires.
- In case of **time limit** violations, the **Flight Module Provider** should provide the event details to the **RFM**.

---

## 5. **Log_takeoff_location**

### **Parameters**:
- **IST date-time**
- **GLPS coordinates**

### **Description**:
- Provides **takeoff event information** to the library for internal logging.
- The **Flight Module Provider** is responsible for ensuring that this API is called immediately after the **takeoff event**.
- The **Flight Module Provider** should not provide any means to bypass the notification to **RFM** on a **takeoff event**.
- **RFM** will store and package all such events in one **flight log** against a **permission artefact**.

---

## 6. **Log_Land_location**

### **Parameters**:
- **IST date-time**
- **GLPS coordinates**

### **Description**:
- Provides **land event information** to **RFM** for internal logging.
- The **Flight Module Provider** is responsible for ensuring that this API is called immediately after the **land event**.
- The **Flight Module Provider** should not provide any means to bypass the notification to **RFM** on a **land event**.
- **RFM** will store and package all such events in one **flight log** against a **permission artefact**.

---

## 7. **Log_geofence_breach**

### **Parameters**:
- **IST date-time**
- **Breach_start_timestamp**
- **Breach_stop_timestamp**
- **GLPS coordinates**

### **Description**:
- Provides **geofence breach event information** to **RFM** for internal logging.
- The **Flight Module Provider** is responsible for ensuring that this API is called immediately after a **geofence breach incidence**.
- The **Flight Module Provider** should call this API at **1 Hz** from the time when the **geofence** is breached until the drone lands.
- **RFM** will store and package all such events in one **flight log** against a **permission artefact**.

---

## 8. **Log_timelimit_breach**

### **Parameters**:
- **IST date-time**
- **time_overrun_start_timestamp**
- **time_overrun_land_timestamp**

### **Description**:
- Provides **time limit breach event information** to **RFM** for internal logging.

### **Revision 1a**:
- The **Flight Module Provider** is responsible for ensuring that this API is called immediately for a **time limit overrun event**.
- During the first call, when the drone is still in flight, the **time_overrun_timestamp** parameter will not have any value. The event will be logged as the **start of time overrun**.
- When the drone lands, this API should be called again to indicate the **end of flight**. In this case, the API caller needs to provide **time_overrun_start_timestamp** again.
- **RFM** will store and package all such events in the **flight log** against a **permission artefact**.

---

## 9. **Get_individual_flight_logs**

### **Parameters**:
- **IST date-time**

### **Description**:
- The **RFM** prepares a **flight log** after each **land**, **geofence breach**, or **flight time overrun**.
- These **flight logs** are **digitally signed** by the **RFM** to ensure that the logs are not tampered with during transport. 
- It is the **Flight Module Provider's** responsibility to implement the functionality required to keep the **private key secure**.
- The **Flight Module Provider** is responsible for collecting the flight logs immediately after every **takeoff**, **land**, **geofence**, and **time-limit breach** event from **RFM**.
- In the case of a **crash** where events were not registered, **RFM**, on the next power cycle, should close the log with a **failed landing incidence**.
- The **Flight Module Provider** is responsible for ensuring that the **flight logs** are stored securely and available for uploading to **Digital Sky Platform (DSP)** or through external applications, such as **Ground Control Station**.

---

## 10. **Bundle_flight_logs**

### **Parameters**:
- **IST timestamp**
- **List_of_individual_incidence_reports_from_storage**
- **Permission Artefact**

### **Outputs**:
- **Bundled_incident_report_with_digital_signature**

### **First Edition**:
- Bundles the signed **flight logs** from storage into a single bundle (one bundle per **permission artefact**).
- After the time period of a **permission artefact** is over, the user must submit the **flight logs** to **Digital Sky APIs** within 3 days.
- The **Flight Module Provider** must provide all the **flight logs** associated with a particular **permission artefact** and pass them on to this API to get a **signed bundle** in return.
- The **Flight Module Provider** is responsible for providing the communication interface to upload the **bundle** directly to the **DSP** or through external applications like **Ground Control Station**.

---

### **Important Notes**:
- All the **flight logs** should be stored and packaged securely by the **RFM**.
- The **Flight Module Provider** has the responsibility to implement appropriate functions for logging, permissions, and interactions with the **RFM** to ensure safe operations.
- **Digital Sky Public Key** is used to verify signatures and permissions for various functions.

---

