/****************************************************************************
 *
 *   Copyright (c) 2018 PX4 Development Team. All rights reserved.
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

/* program for writing individual log file

Given: 1) Log entries Array of Geo_tag
       2) Permission artefact ID
       3) previous log hash
       4) Private Key

Output: Signed Json flight Log

'''''''''''''''''''''''''''''''''''''
{
    "FlightLog"=
        {
            "PermissionArtefact"=string;
            "previous_log_hash"=string;
            "LogEntries"=[
                {},
                {},
            ]
        },
    "Signature"="hash"
}
*/


//////////////////////////////////////////////
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
#include <drivers/drv_hrt.h>
#include<motion_planning/main_utility.hpp>
#include<motion_planning/print_hello.hpp>


using namespace std;



struct Geo_tag{

    int Entrytype ;
    long Timestamp;
    float Longitude;
    float Lattitude;
    int Altitude;

};


double max(double a,double b){
    if(a<b){
        return b;
    }
    else{
        return a;
    }

}


int fetch_previous_log_hash(char *result){
    // this function will fetch hash value from recentPA.txt
    FILE *fptr;
    fptr=fopen("./log/recentPA.txt","r");
    signed char ch;
    int i=0;
    char *content=(char*)malloc(sizeof(char)*5000);
    while((ch=fgetc(fptr))!=EOF){
        content[i]=ch;
        i++;
    }
    content[i]='\0';
     fclose(fptr);
    char pr_log_hash_tag[50]="previous_log_hash";
    char pr_log_hash_contain[100];

    fetch_tag(content,pr_log_hash_tag,pr_log_hash_contain);
    free(content);
    if(strcmp( pr_log_hash_contain,"None")==0){
        // this is the first log for a particular PA.
        return 0;

    }else{
        strcpy(result,pr_log_hash_contain);
        return 1;
    }


}



void signing_log_content(char *HEX,char *cipher_string_hex){
      printf("\nproblem at line 668 possible mp iniialization %s\n",HEX);

    // modulus in public key
   mp_int modulus;
   mp_init(&modulus);
   char Modulus[513] ="c8f770862ad1bc391310984fde3bb0879590e44e25dfbf9177431b27af55bb21bb6e93358d22d6ea60a38ce9e08cb30077e3d8e6908a2cc46bde10e767c0cc8977d66322be0fc196e5f399675a22a666563ce7de7af613c9dc23f79df8abd8e99391d957de97e9950bd2e1a8671404e8a9d19c32cae9dd5f6d29474cf45b613a31596bead1db45c7a929994fc99bdf68426ab3e9857b61fee826c01ee81c9fc030754252e26822ee26d5931b8835d221cc7889cdefec20ccc4312934d52b4de3cdfcb258dfefbf73ad365fb04a18edd3b8d3f68196fb910b791f6dbcd60b5c25b7f646176f9b009305c47d1b6c2cfe1ebddb58c08ad2da9dc932851ba9e258ef";
	//char Modulus[513] = "bc07d529450214ef63a8d61966987e8ca0594d9a7ec4f1881117b4f8ecbdc74b8769f6c98bfe931c9474116be8bd36527acfd95f6633d12cc8a960ab3d3e7a0b4b3e4990b594ee61af3b56315337501225525fb997b65c38118d614601dcb8bd631673a510498f2c3dab44d723d8b6daa697d0108e7fcb4d27525f386e7fcd9ce29c4ab12c4258aa77872259a25804791a1eaef54b65226ec84765442ac839db30467d86910e700d802807de1f4fef5235738d66359cb0a2707cb9cd90e90bb1f2d0d807aafbd048b1ddbb156d4984cfbbaa9a435b9230d213140dd5be64b5e594945d474665eaf5267fc598a5f75b99f83b029971b80c4149891d43abe62b95";
	//char Modulus[513]="cfb7e3ed5fb094ca81a2da9b07db403fc3882a33acf45f34c57d8fe677b3c787e3f034b4cfeba70f8a7beb95a53520780a59b1494ce00393aa812edf96abdafab551dd6728cd5d4ac4371c6049acecc6b63972737e60773129fa86dd8c5277ed64d13febb808776f630f7afb895789f609b28ba37caf0b14381f48ab123c093dc2621c34c570144d912de27fbf1d5669d427a07c9cf027013c93167d7faab8e8b4c95cfd0a9be0e4a92ae132713891a63af36315fce8d6f376cd2be410ef2155be9182caba4772275c93b47d5b8eed6cfe3fda8be884f79d49de09d02e42d93f350a25486a90653b0c72b6fb4d7bf9d40a2bc1559785b4bd0dd24a70aca7ef93";
   mp_read_radix(&modulus,Modulus,16);

    //private key exponent
   mp_int Private_key;
   mp_init(&Private_key);
   char private_hex[513]="184885e9405d4d801c04a252ec488c2125fa770bd659bdfd26cb0e09f28eca68de0c136fa219369ce5867dad78fba75984231cff6731bb0d14f7a55540dd3419dc48247c7b38ce2c9ca69dbfb64d7f8bd819cdeebd2ee4df3c6180372f681c72c4e917b91d657fcd09bbb696b1b5e28df68f246fa2c33583a55e1a867af45bc004583335a5d5a79d046b2d9c74bb5b7282a4ffb0167c3f0c982051fe6e253a8c20d063e1eccf0dfcdcaf4c34a3f0ca2cf8ebbb37d03c8176ed157b55f738c5f47fe4d0d42b05c6c0e9a58cc3766ec87926ee8e38843b25171886ec413c644b931acbb92895c4665274a17355d1d72401339d59ee0ec69a65d613e9b5c215a7c9";
  // char private_hex[513] = "9b78ea83264133684182400d5eaca6aec68330cc97176712f7f71f3758210f61df44f9beead78372753987922f2e0c75a480aa1edc95e9d65ad0da529ce044ef83b6ac03507125ae75c2dd61098ac9d54730d65fd21702278633dd8392549c18548f22ee100a92aca50d316da68131a897691dac22f77df57c96fa8ee1a7212db313396410a5c9c8a31f6f940724cff2b2db5eb078eedad92b6ff29a8636fcd370e99773e96168f34839693f84b7a083597bfbe0f674c79b2348b038ca730cada30bcf2dd9cde27dd555891d3cc10b7831b23e7cda163570635727f11d569492a201f55c56d9a92d46f71b6ecea30f28f8c040f834a2da43f72a1ec927df9441";
   //char private_hex[513]="6beefbcaae7c4cf46524403f6a77ad0cf5075e1677fa8b361aa0c2135983db5c6b3eb7c4747dd8d3247c7bcfc886b0966f9a679ad50d5a0e72fca964992037ab2a689d892b147b338c7dae8b01fd8f133a40e38dcbcf48600d96165a2cbdf57f2f71e3ab1277a3c8074b55f63a497870965d665dcf3e0d9db603db78b902e5317897c3538a24057c29f25f98284c7262bd965f4923da6a117f37c88e78ec6ce34da38eae946c06cc4ac166df70a7dd24b1be2254c8a25797bde3c58918b9f8ae0b5bf8019f8f60fb9e7e087d20bcb935b0b9c3a43c0b7f64ca9027c8bfa9e83525781dfd7b529a88a305ae8b4ec56e2cadff1ded53f2ec461ae2bd56e26231f1";
   mp_read_radix(&Private_key,private_hex,16);//
   //char yu[]="b59cab400e9c64525f566e85c97c2d65a2d0620f0479cb8ef8e6eb1d3c1cf093be5eb827a4d01d9a585b0043f148d9e3c213676d819818be05bf39a7657cd953bd4a7b2e361d03e1d64d99c755036eb1f97223832b1d596ee131f0adb84879cbd6763f1fbe61c214ff6514551900ef6034feb20cf0daf995ebaf9d7e5c3566cf012e1aa32e3733747be37955ec0d127d32fbddfdd8451fc5254dfe09a4b1b26bd9ec3487cd04121fd4b965456c024e21086da38c81670b8cb57fdad93f8f90eb6bfa109d57ffda156c105fe15c4c78d588de116ff06abc1ea03546ef36b1ac711cfe4feb64250264ac5ff2704a76a6b1fe04480af906e9ee71836190c345bf67";
    /*const char *content_to_hash =(const char*) malloc(sizeof(char)*4000);

    int rsa_encryption_count_aux=0;
    while(HEX_format_Digest[rsa_encryption_count_aux]!='\0'){
        content_to_hash[rsa_encryption_count_aux]=HEX_format_Digest[rsa_encryption_count_aux];
        rsa_encryption_count_aux++;
    }*/
   //making the signing process pkcs.1.15 compatible
   char padding_SHA256[810]="1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff003031300d060960864801650304020105000420";
   strcat(padding_SHA256,HEX);
   printf("\npadded sha digest :%s\n",padding_SHA256);
   mp_int hash_to_sign;
   mp_init(&hash_to_sign);

   mp_read_radix(&hash_to_sign,padding_SHA256,16);
   char *aux_hash_ex=(char*) malloc(1014*sizeof(char));
   memset(aux_hash_ex,0,1014);
  // char aux_hash_ex[1014];
   mp_to_hex(&hash_to_sign,aux_hash_ex,1014);
   printf("\n   hash to sign ===\n%s\n\n",aux_hash_ex);
   free(aux_hash_ex);

   //Signing begins

   mp_int cipher;
   mp_init(&cipher);

   mp_exptmod(&hash_to_sign,&Private_key,&modulus,&cipher);
  //
   //char cipher_string[1013];
   char *cipher_string=(char*) malloc(1013*sizeof(char));
   mp_to_hex(&cipher,cipher_string,1013);
   printf("\ncipher string in hex %s\n",cipher_string);//this is the encrypted message
    strcpy(cipher_string_hex,cipher_string);
    free(cipher_string);
    mp_clear(&modulus);
     mp_clear(&Private_key);
      mp_clear(&hash_to_sign);
    mp_clear(&cipher);


}


void get_PA_ID(char *res){
    FILE *fptr;
    fptr=fopen("./log/recentPA.txt","r");
    if(fptr==NULL){
        printf("\nfile is not opening\n");
    }
    //char content[5000];
    char *content=(char*)malloc(sizeof(char)*5000);
    signed char ch;
    int i=0;

    while((ch=fgetc(fptr))!=EOF){
        content[i]=ch;
        i++;
    }
    content[i]='\0';
    fclose(fptr);
    char tag[50]="permissionArtifactId";

    fetch_tag(content,tag,res);
    free(content);


}

void writingKey_Value(char *str_ptr,char *Key,int *count,char *value,int last){

    int i=*count;

    str_ptr[i]='"';
    i++;

    for(int y=0;y<int(strlen(Key));y++){
        str_ptr[i]=*(Key+y);
        i++;
    }
    str_ptr[i]='"';
    i++;
    str_ptr[i]=':';
    i++;
    str_ptr[i]=' ';
    i++;


    //this is for string

    str_ptr[i]='"';
    i++;
    for(int y=0;y<int(strlen((char*)value));y++){
        str_ptr[i]=*((char*)value+y);
        i++;
    }
    str_ptr[i]='"';
    i++;


    if(last==1){
        str_ptr[i]='\n';
        i++;
    }else{
        str_ptr[i]=',';
        i++;
         str_ptr[i]='\n';
        i++;
    }
    *count=i;

}


void writingKey_Value(char *str_ptr,char *Key,int *count,int value,int last){

    int i=*count;

    str_ptr[i]='"';
    i++;

    for(int y=0;y<int(strlen(Key));y++){
        str_ptr[i]=*(Key+y);
        i++;
    }
    str_ptr[i]='"';
    i++;
    str_ptr[i]=':';
    i++;
    str_ptr[i]=' ';
    i++;

    //this is for integer
    char str[30];
    sprintf(str, "%d", value);

    //putting char from str to str_ptr
    int i_count=0;
    while(str[i_count]!='\0'){
    str_ptr[i]=str[i_count];
    i_count++;
    i++;
    }



    if(last==1){
        str_ptr[i]='\n';
        i++;
    }else{
        str_ptr[i]=',';
        i++;
         str_ptr[i]='\n';
        i++;
    }
    *count=i;

}


void writingKey_Value(char *str_ptr,char *Key,int *count,long value,int last){

    int i=*count;

    str_ptr[i]='"';
    i++;

    for(int y=0;y<int(strlen(Key));y++){
        str_ptr[i]=*(Key+y);
        i++;
    }
    str_ptr[i]='"';
    i++;
    str_ptr[i]=':';
    i++;
    str_ptr[i]=' ';
    i++;

    //this is for integer
    char str[30];
    sprintf(str, "%ld", value);
    printf("\n\n\n string unsigned int  ::::::::::  %s\n\n\n",str);

    //putting char from str to str_ptr
    int i_count=0;
    while(str[i_count]!='\0'){
    str_ptr[i]=str[i_count];
    i_count++;
    i++;
    }



    if(last==1){
        str_ptr[i]='\n';
        i++;
    }else{
        str_ptr[i]=',';
        i++;
         str_ptr[i]='\n';
        i++;
    }
    *count=i;

}


void writingKey_Value(char *str_ptr,char *Key,int *count,double value,int last){

    int i=*count;

    str_ptr[i]='"';
    i++;

    for(int y=0;y<int(strlen(Key));y++){
        str_ptr[i]=*(Key+y);
        i++;
    }
    str_ptr[i]='"';
    i++;
    str_ptr[i]=':';
    i++;
    str_ptr[i]=' ';
    i++;


    // this is for float
    char str[30];
    sprintf(str, "%f", value);

    //putting char from str to str_ptr
    int i_count=0;
    while(str[i_count]!='\0'){
    str_ptr[i]=str[i_count];
    i_count++;
    i++;
    }




    if(last==1){
        str_ptr[i]='\n';
        i++;
    }else{
        str_ptr[i]=',';
        i++;
         str_ptr[i]='\n';
        i++;
    }
    *count=i;

}


void writingKey(char *str_ptr,char *Key,int *count){

    int i=*count;

    str_ptr[i]='"';
    i++;

    for(int y=0;y<int(strlen(Key));y++){
        str_ptr[i]=*(Key+y);
        i++;
    }
    str_ptr[i]='"';
    i++;
    str_ptr[i]=':';
    i++;
    str_ptr[i]=' ';
    i++;
    *count=i;

}



void putSpace(int count_Space,char *ptr_string ,int *i_ptr){
    int y=*i_ptr;
    for(int u=0;u<count_Space;u++){
        ptr_string[y]=' ';
        y++;
    }
    *i_ptr=y;
}

void objectCreator(char *ptr,Geo_tag geo,int last,int space_count,int *count){

    /// this function will make object entries inside LogEntries
    /*
    {
        "Entry_type": "TAKEOFF/ARM",
        "TimeStamp": 41,
        "Longitude": 380.000000,
        "Latitude": 443.5000000,
        "Altitude": 704.3400,
    }
    */


    char Geolat[]="Latitude";// float  Latitude in Degrees East
    char Geolon[]="Longitude";//float  Longitude in Degrees North
    char GeoTime[]="TimeStamp";//milliseconds, type : integer
    char  GeoAlt[]="Altitude";//Ellipsoidal Height in Meters  type: integer
    char GeoEntryType[]="Entry_type";//entry_type

    int iterator=*count;

    ptr[iterator]='{';
    iterator++;
    ptr[iterator]='\n';
    iterator++;
    space_count=space_count+2;
    // putting spaces
    putSpace(space_count,ptr,&iterator);
      /// Lattitude
    writingKey_Value(ptr,Geolat,&iterator,geo.Lattitude,0);
    putSpace(space_count,ptr,&iterator);

    // starting with entrytype
    char entry[20];
    if(geo.Entrytype==0){
         strcpy(entry,"GEOFENCE_BREACH");
    }else if(geo.Entrytype==1){
         strcpy(entry,"TAKEOFF/ARM");
    }else if(geo.Entrytype==2){
         strcpy(entry,"TIME_BREACH");
    }else if(geo.Entrytype==3){
         strcpy(entry,"LAND");
    }else{
        strcpy(entry,"CRASH");
    }
    writingKey_Value(ptr,GeoEntryType,&iterator,entry,0);
    putSpace(space_count,ptr,&iterator);

    /// Altitude
    writingKey_Value(ptr,GeoAlt,&iterator,geo.Altitude,0);

    putSpace(space_count,ptr,&iterator);
    /// Longitude
    writingKey_Value(ptr,Geolon,&iterator,geo.Longitude,0);
    putSpace(space_count,ptr,&iterator);

     /// TimeStamp
    writingKey_Value(ptr,GeoTime,&iterator,geo.Timestamp,1);
    space_count=space_count-2;
    putSpace(space_count,ptr,&iterator);


    //closing

    ptr[iterator]='}';
    iterator++;

    *count=iterator;


}



void  log_naming_support(char *paID_firstTerm,char *done_freq){

    FILE *fptr;

    fptr=fopen("./log/recentPA.txt","r");
    signed char ch;
    int i=0;
    char *content=(char*)malloc(sizeof(char)*5000);
   // char content[5000];

    while((ch=fgetc(fptr))!=EOF){
        content[i]=ch;
        i++;
    }
    content[i]='\0';
    fclose(fptr);
    char tag_paID[50]="permissionArtifactId";
    char tag_paID_contain[80];
    char tag_done_freq[50]="frequencies_done";
    char tag_done_freq_contain[100];

    fetch_tag(content,tag_paID,tag_paID_contain);
    fetch_tag(content,tag_done_freq,tag_done_freq_contain);
    free(content);
    i=0;
    while(tag_paID_contain[i]!='-'){
        paID_firstTerm[i]=tag_paID_contain[i];
        i++;
    }
    paID_firstTerm[i]='\0';

    strcpy(done_freq,tag_done_freq_contain);

}

void main_json_file_writing(Geo_tag *data_array, int length_dat){
    char *buff=(char*) malloc(sizeof(char)*10000);
    //char buff[10000];// this is where json would be written, later copied to a file

    int space_count=0;
    //object {}
    char object_start='{';
  //  char object_end='}';

    //array []
    char array_start='[';
    //char array_end=']';
    // ""
   // char quote_start='"';
   // char quote_end='"';
    //colon :
  //  char colon=':';
    //comma
   // char comma=',';


    // writting JSON flight log
    int i=0;

    buff[i]=object_start;
    space_count=space_count+2;
    i++;
    buff[i]='\n';
    i++;
    // after { and new line spaces need to be put

    putSpace(space_count,&buff[0],&i);

    // starting with writing Flight log
    int write_signature_help=i;
    char firstKey[]="FlightLog";

    writingKey(&buff[0],firstKey,&i);

    int str_start=i;

    buff[i]=object_start;
    i++;
    buff[i]='\n';
    i++;
    space_count=space_count+2;

    // after { and new line, spaces need to be put

    putSpace(space_count,&buff[0],&i);



    ////// LogEntries////
    char Log[]="LogEntries";

    writingKey(&buff[0],Log,&i);

    buff[i]=array_start;
    i++;
    buff[i]='\n';
    i++;
    space_count=space_count+2;

    // after { and new line, spaces need to be put

    putSpace(space_count,&buff[0],&i);
    // now a loop will run
    for (int counter=0;counter<length_dat;counter++){
        objectCreator(&buff[0],data_array[counter],0,space_count,&i);
        if(counter==length_dat-1)
        {
            buff[i]='\n';
            i++;
            space_count=space_count-2;
        }else{
            buff[i]=',';
            i++;
            buff[i]='\n';
            i++;
        }

        putSpace(space_count,&buff[0],&i);
    }
   /* space_count=space_count-2;
    putSpace(space_count,&buff[0],&i);*/
    buff[i]=']';
    i++;
    buff[i]=',';
    i++;
    buff[i]='\n';
    i++;
    putSpace(space_count,&buff[0],&i);


    printf("\nproblem at line 516\n");


    // writting permission artefact id from recentPA.txt


    char PA_id[70];//="ABDHUBDUBBEBDIBWBID";

    get_PA_ID(PA_id);
printf("\nproblem at 525\n");
    char PA[]="PermissionArtefact";
    writingKey_Value(&buff[0],PA,&i,PA_id,0);
    printf("\nproblem at 528\n");
    putSpace(space_count,&buff[0],&i);
    printf("\nproblem at 530\n");


    // fetching previous log hash from recentPA.txt
    // if None then its the first log for the particular PA
    // writing previous log hash
    char previous_log_hash[70];//="aaaanfbofburbfubvoinvknvofn";
    char previous_log_hash_2[70];//="aaaanfbofburbfubvoinvknvofn";
    printf("\nproblem at 536\n");
    int hash_status=fetch_previous_log_hash(previous_log_hash_2);

    if(hash_status==0){
        // no previous log hash is present.
        strcpy(previous_log_hash,"First_log");

    }else{
        strcpy(previous_log_hash,previous_log_hash_2);
    }

    char prev_hash[]="previous_log_hash";
    writingKey_Value(&buff[0],prev_hash,&i,previous_log_hash,1);



   // putSpace(space_count,&buff[0],&i);

  printf("\nproblem at line 554\n");
   space_count=space_count-2;
   putSpace(space_count,&buff[0],&i);
   int str_end=i;/////// the point where value of flightlog key ends
   buff[i]='}';
   i++;
   buff[i]='\n';
   i++;
   space_count=space_count-2;
   putSpace(space_count,&buff[0],&i);
   buff[i]='}';
   i++;
   buff[i]='\0';
    // At this point, json file is completed
    // next step is to canonicalize the value(object) under "FLightlog" key
   char *canonicalized_flight=(char*) malloc(sizeof(char)*14000);
  // char canonicalized_flight[14000];
   int canon_count=0;


   for(int h=str_start;h<=str_end;h++)
    {
    /* canonicalize hack : all the serialization under objects have been taken care of
    earlier only, now just keep on taking the elements one by one from one array to another
    except new lines and spaces
    */
        if (buff[h]!=' ' && buff[h]!='\n'){
            canonicalized_flight[canon_count]=buff[h];
            canon_count++;
        }

    }
    canonicalized_flight[canon_count]='\0';


    // at this point canonicalization has been done
    // our next step is to implement SHA256 algo

    signed char ch;
    unsigned char *st=new unsigned char[10000];
    unsigned int i_sha=0;
    //// assigning char to unsigned char, this has to be done to implement sha256
    while(1){
	   ch=canonicalized_flight[i_sha];
      // printf("%c",ch);
	   if(ch=='\0'){
		   break;
	   }
        *(st+i_sha)=ch;
		printf("%c",*(st+i_sha));
		i_sha++;

    }
  free(canonicalized_flight);
  printf("\nproblem at line 607\n");

   SHA256 sha;

	sha.update(st,i_sha);

	uint8_t * digest = sha.digest();
    int it=0;//just an iterator
	//printing digest
	printf("\n");

    while(it<32){
		printf("%02x ",*(digest+it));
		it++;
	}
    char HEX_format_Digest[65];

    // code for making hex string
    int hex_count=0;
    for(int yui=0;yui<32;yui++){
         //printf("%d",*(digest+yui)>>4||0x0f);
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui)>>4 & 0x0f];
         hex_count++;
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui) & 0x0f];
         hex_count++;
      }
   HEX_format_Digest[hex_count]='\0';

   printf("\nhere is the string:%s\n\n",HEX_format_Digest);

   // this is the hash of the current log file, this will be updated in recentPA.txt.
   //
   update_recentPA(2,HEX_format_Digest);

    //At this point we have got the string format for the digest


//	delete[] digest;

   delete []st;
   st=NULL;
  // printf("\n\n");
   // printf("%s",canonicalized_flight);
   //printf("\n\n");
    //printf("%s",buff);
  // printf("\n");

/*
   FILE *fptr;
   fptr=fopen("flight_log1.json","w");// this is the resultant json file
   fprintf(fptr, "%s", buff);
   fclose(fptr);
*/


    // Now RSA encryption is required
    /* following elements are required
    1)Private key
    2) above calculated digest HEX_format_digest
    */
    //////////////////////// RSA Encryption ///////////////
    char *cipher_string_hex=(char*) malloc(1013*sizeof(char));

    signing_log_content(HEX_format_Digest,cipher_string_hex);





    //char ty[]="734CB799A64A670D68F3AA84B2542BBEE1F7FD5AC022460CDD63E297DF55D8B04CBCB1103BDC720EB85D090CCA091D067A5F54BFADBECB09A2EBDDA9D92E00AF83D5E9FED29B38C295A315CBD52CB4ED1645BC3707D4ED165E5E4D906A1C149E5073AE6A8D5FC74CC11D63885CEFE236E4464C0D387A5861234814FF31A8A3A0B2D6022B8A8FAEC79DF61B638B0B85392F152D4743BFD2779C4472B73CF952027C83A8DB7EEA32D6A1B96ED8F92CF7576EBBB35F8CA341C69E60C352BD79628C17585030D8BF07C8E6D73FC0AB51D2226B1E7C42D5DA0BC4EAF5063C8F72601AC0090D92B74C17365463B3F07B0397692AF204F58A0ACB665CBA4C336121C4CA";
    /*
    mp_int Public_key;
    mp_read_radix(&Public_key,"65537",10);
    mp_int message;
    mp_exptmod(&cipher, &Public_key,&modulus,&message);                          this part for decrypting encrypted text

   char message_string_hex[513];
    mp_to_hex(&message,message_string_hex,sizeof(message_string_hex));
    printf("\n%s\n",message_string_hex);//this is the decrypted message
    */


   /////////////////////// base64 Encoded //////////////// hex to base64
   char *res=(char*) malloc(sizeof(char)*2500);
   //char res[2500];
   base64Encoder(cipher_string_hex,strlen(cipher_string_hex),res);
   free(cipher_string_hex);
   printf("\n\nbase 64 encoder %s\n",res);
   char signature[700]="\"signature\": \"";
   strcat(signature,res);
   free(res);
   char t[4];
   char *gio=&t[0];
   t[0]='"';
   t[1]=',';
   t[2]='\0';
   strcat(signature,gio);
   printf("\n%s\n",signature);
   printf("the length of string %d",(int)strlen(signature));


   /////////////////////// Writing Signature ///////////////// open file and rewrite
   /* string insertion
   1) length of string that has to be inserted  l1   res
   2) length of string  where we want to insert  l2  buff
   */
    int i_sign=0;
    char *buff1=(char*) malloc(sizeof(char)*15000);
    //char buff1[15000];
    space_count=0;
    buff1[i_sign]=object_start;
    space_count=space_count+4;
    i_sign++;
    buff1[i_sign]='\n';
    i_sign++;
    // after { and new line, spaces need to be put

    putSpace(space_count,&buff1[0],&i_sign);
    int you=0;
    while(signature[you]!='\0'){
      buff1[i_sign]=signature[you];
      i_sign++;
      you++;
   }
   buff1[i_sign]='\n';
   i_sign++;
   putSpace(4,&buff1[0],&i_sign);
   for(int iu=write_signature_help;iu<int(strlen(buff));iu++){
      buff1[i_sign]=buff[iu];
      i_sign++;

   }
   free(buff);
   buff1[i_sign]='\n';
   i_sign++;
   buff1[i_sign]='\0';
   printf("\n\n%s\n\n",buff1);

   //deciding file name for logs
   // current approach : take note of the done frequencies inside recentPA.txt
   //and last four letters of pa id first term.
   //filename ==== paidfirst_term-log-(done_frequency+1)

   char paID_firstTerm[20];
   char done_freq[10];

   log_naming_support(paID_firstTerm,done_freq);

   char file_directory[10]="./log/";

   char filename[100];

   strcpy(filename,file_directory);
   strcat(filename,paID_firstTerm);
   strcat(filename,"_log_");
   char aux_start_filename[40];
   strcpy(aux_start_filename,filename);

   strcat(filename,done_freq);
   strcat(filename,".json");

   printf("%s\n\n",filename);

   update_log_of_logs(filename,HEX_format_Digest,done_freq,PA_id,aux_start_filename);

   FILE *fptr1;
   fptr1=fopen(filename,"w");// this is the resultant json file
   fprintf(fptr1, "%s", buff1);
   fclose(fptr1);
   free(buff1);

}



int check_Geobreach(pa_data_s data,vehicle_global_position_s vgp){
    // return 1 if geofence breached
    // return 0 if geofence not breached
    double lattitude1 =data.lattitude1;
    double longitude1 =data.longitude1;
    double lattitude2 =data.lattitude2;
    double longitude2 =data.longitude2;
    int home_altitude =data.home_altitude;

    double mini_lat=min(lattitude1,lattitude2);
    double maxi_lat=max(lattitude1,lattitude2);

    double mini_lon=min(longitude1,longitude2);
    double maxi_lon=max(longitude1,longitude2);
    double allowed_height=0.3048*(data.allowable_height);// converting feets to meters


    double main_lat=vgp.lat;
    double main_lon=vgp.lon;
    double main_alt=vgp.alt;

    if(main_lat>=mini_lat && main_lat<=maxi_lat){
        if(main_lon>=mini_lon && main_lon<=maxi_lon){
            if((main_alt-(double)home_altitude)>allowed_height){
                return 1;

            }else{
                return 0;
            }

        }else{
            return 1;
        }


    }else{
        return 1;
    }





}


int check_Timebreach(pa_data_s data,vehicle_gps_position_s vgp){
    // return 1 timebreach
    // return 0 no timebreach
    Date_time start;
    Date_time end;

    start.date=data.start_date;
    start.Hours=data.start_hours;
    start.Minutes=data.start_minutes;
    start.Month=data.start_month;
    start.Seconds=data.start_seconds;
    start.Year=data.start_year;

    end.date=data.end_date;
    end.Hours=data.end_hours;
    end.Minutes=data.end_minutes;
    end.Month=data.end_month;
    end.Seconds=data.end_seconds;
    end.Year=data.end_year;



    Date_time current;

    time_t timestamp;
    timestamp=(vgp.time_utc_usec)/1000000;//microsecond to seconds
	//printf("\ntimestamp is seconds: %u\n",timestamp);
    struct tm  ts;

    ts = *localtime(&timestamp);
    current.Year=ts.tm_year+1900;
    current.Month=ts.tm_mon+1;
    current.date=ts.tm_mday;
	current.Hours=ts.tm_hour;
	current.Minutes=ts.tm_min;
	current.Seconds=ts.tm_sec;

    int status=In_Time(current,start,end);

    if(status==0){
        return 1;

    }
    return 0;




}
