#include<stdio.h>
#include"MP_INT.hpp"
#include<string.h>
#include<inttypes.h>
using namespace std;

struct pair_set{
    char tag[30];
    char value[1200];
};

struct key{
char modulus[1000];
char private_exponent[1000];
};



//true value of this number stays in HArdwareInuse.txt, that is compared to the one that is created at each start up
int HardwareNumber;



// functions(next ) will write the file as per the given (tag and values)_pair and will also sign it

char HEX_ele[17]="0123456789abcdef";

void fetch_tag(char *content, char *target, char *result);

int isSubstring(char *s1, char *s2)///s1 is the sub string ; s2 is the larger string
{
    int M = strlen(s1);
 //  printf("11 %s jkl\n",s1);
    int N = strlen(s2);
  // printf("22 %s jkl\n",s2);

    for (int i = 0; i <= N - M; i++) {
        int j;

        for (j = 0; j < M; j++)
            if (s2[i + j] != s1[j])
                break;

        if (j == M)
            return i;
    }
    return -1;
}

int find_int_hex(char ptr){
    char char_set16[]="0123456789ABCDEF";
    int i=0;
    while(1){
        if(char_set16[i]==ptr){
            return i;
        }
        i++;
    }

}
void base64Encoder(char input_str[], int len_str, char *result)
{

    //character set of base 64 encoding scheme
    char char_set[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    char char_set16[]="0123456789ABCDEF";

  //  char *res_str=(char *) malloc(1000 * sizeof(char));//chr *res_str =new char(1000*1)

    if(len_str%2!=0){
        printf("the hex string is not valid for having conversion to base64 digits");
        exit(1);
    }

    int i=len_str;

    /* basic agenda
    an initial check that number of char in string are even or not (gouping of bytes(8 bits= 2 hex char) is done)
    first check if 6 continous hex digits could be taken or not
    if yes:
    1)  Take 6 hex digits :24 bits
    2)  24/6=4 then write 4 base64 characters
    if no:
    1)check how many hex char are available(<6)
    2)possible cases 2,4
    if 2 hex char available:
        -add 4 bits(valued 0) to make the bit length divisible by 6
    if 4 hex char available
        - add 2 bits(valued 0) to make the bit length divisible by 6
    */
    int j=0;
    int res_count=0;
    while(1){
        if(i>=6){//if yes code starts here

            int aux=0;
           // printf("%d\n",res_count);
            for(int index=j ; index<j+6;index++){
                aux=(aux<<4);
              //  printf("%d ",find_int_hex(input_str[index]));
                aux=aux|find_int_hex(input_str[index]);
              //  printf(" %08x ",aux);
            }
            j=j+6;
           // printf("\n");
            int aux1=0;
            for(int o=3;o>=0;o--){
                aux1=aux & 0x3f;
                result[res_count+o]=char_set[aux1];
                aux=aux>>6;
                //res_count--;
            }
            res_count=res_count+4;
            i=i-6;
        }else{

                if(i==0){
               //     printf("Done...");
                 //   printf("\n");
                    break;
                }
                if(i==2){
                    int aux=0;
                    for(int index=j ; index<j+2;index++){
                            aux=(aux<<4);
                        //  printf("%d ",find_int_hex(input_str[index]));
                            aux=aux|find_int_hex(input_str[index]);
                            //  printf(" %08x ",aux);
                        }
                        aux=aux<<4;
                    int aux1=0;
                    for(int o=1;o>=0;o--){
                        aux1=aux & 0x3f;
                        result[res_count+o]=char_set[aux1];
                        aux=aux>>6;
                    //res_count--;
                    }
                    res_count=res_count+2;
                    result[res_count]='=';
                    res_count++;
                    result[res_count]='=';
                    res_count++;
                    result[res_count]='\0';
                    break;

                }else{

                    int aux=0;
                    for(int index=j ; index<j+4;index++){
                            aux=(aux<<4);
                        //  printf("%d ",find_int_hex(input_str[index]));
                            aux=aux|find_int_hex(input_str[index]);
                            //  printf(" %08x ",aux);
                        }
                        aux=aux<<2;
                    int aux1=0;
                    for(int o=2;o>=0;o--){
                        aux1=aux & 0x3f;
                        result[res_count+o]=char_set[aux1];
                        aux=aux>>6;
                    //res_count--;
                    }
                    res_count=res_count+3;
                    result[res_count]='=';
                    res_count++;
                    result[res_count]='\0';
                    break;


                }

        }


    }


}


// signs the content with signKey(a key) and puts the result inside result
void signing_support(key signKey,char *content,char*result){
    printf("\n\n%s\n\n",content);
    char ch;
    unsigned char *st=new unsigned char[3000];
    unsigned int i_sha=0;
    //// assigning char to unsigned char, this has to be done to implement sha256
    int i_sha_a=0;
    while(1){

	   ch=content[i_sha_a];
      // printf("%c",ch);
	   if(ch=='\0'){
		   break;
	   }
       if(ch!='\n'){
        *(st+i_sha)=ch;
	//	printf("%c",*(st+i_sha));
		i_sha++;}
        i_sha_a++;

    }

    printf("ooooLLL\n\n%s\n\n",st);

    SHA256 sha;

	sha.update(st,i_sha);

	uint8_t * digest = sha.digest();
    int it=0;//just an iterator

    char HEX_format_Digest[65];

    // code for making hex string
    int hex_count=0;
    for(int yui=0;yui<32;yui++){
     //    printf("%d ",*(digest+yui)>>4||0x0f);
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui)>>4 & 0x0f];
         hex_count++;
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui) & 0x0f];
         hex_count++;
      }
    delete[] digest;
    digest=NULL;
    HEX_format_Digest[hex_count]='\0';
static int pass=0;
pass++;
    printf("\nhere is the string %d:%s\n\n",pass,HEX_format_Digest);

    //At this point we have got the string format for the digest


    delete []st;
    st=NULL;
    ////////////////////

    // modulus in public key

    mp_int modulus;
    mp_init(&modulus);
    mp_read_radix(&modulus,signKey.modulus,10);

    //private key exponent
    mp_int Private_key;
    mp_init(&Private_key);
    mp_read_radix(&Private_key,signKey.private_exponent,10);

   //making the signing process pkcs.1.15 compatible
   char padding_SHA256[446]="1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff003031300d060960864801650304020105000420";
   strcat(padding_SHA256,HEX_format_Digest);
   //printf("\npadded sha digest :%s\n",padding_SHA256);
   mp_int hash_to_sign;
   mp_init(&hash_to_sign);
   mp_read_radix(&hash_to_sign,padding_SHA256,16);
   char aux_hash_ex[514];
   mp_to_hex(&hash_to_sign,aux_hash_ex,sizeof(aux_hash_ex));
   printf("\n   hash to sign ===\n%s\n\n",aux_hash_ex);

   ////Signing begins

   mp_int cipher;
   mp_init(&cipher);

   mp_exptmod(&hash_to_sign,&Private_key,&modulus,&cipher);

   char cipher_string_hex[513];
   mp_to_hex(&cipher,cipher_string_hex,sizeof(cipher_string_hex));
   printf("\n%s\n",cipher_string_hex);//this is the encrypted message


    base64Encoder(cipher_string_hex,strlen(cipher_string_hex),result);


mp_clear(&hash_to_sign);
mp_clear(&modulus);
mp_clear(&cipher);
mp_clear(&Private_key);

}

void pair_file_write(pair_set *ptr,int ptr_quant ,char *file_name , key Skey){


    char content[6000];//=(char*) malloc(sizeof(char)*6000);
    strcpy(content,"<content>\n");

    int i=0;
    while(i<ptr_quant){

    strcat(content,"<");
    strcat(content,ptr[i].tag);
    strcat(content,"=");
    strcat(content,ptr[i].value);
    strcat(content,">\n");
    i++;
    }
    strcat(content,"</content>\n");

    strcat(content,"\0");

    printf("\nINside pair_file_write %s \n",content);

    // now signing takes place
    char signature[2000];//=(char*) malloc(sizeof(char)*2000);
    signing_support( Skey,content,signature);
    content[int(strlen(content))-1]='\n';
    strcat(content,"<Sign>\n");
    strcat(content,"<Signature=");
    strcat(content,signature);
    strcat(content,">\n");
    strcat(content,"</Sign>");
    strcat(content,"\0");




    ////////////////////

    FILE *fptr;
    fptr=fopen(file_name,"w");
    fprintf(fptr, "%s", content);
    fclose(fptr);
    //free(content);
  //  free(signature);
}

//this function is for validating the txt file wrt  the given key

int inBase64(char *d){
    for(int i=0;i<64;i++){
        if(*d==BASE64_DIGITS[i]){
            return i;
        }
    }
}

void base64decoder(char *base64_ptr,char *hex){



    int ik=0;



    int aux1,aux2,aux3,num_bits;

    for(int i=0;base64_ptr[i]!='\0';i=i+4){
        aux1=0;num_bits=0;

        for(int j=0;j<4;j++){

        if(base64_ptr[i+j]!='='){
        aux1=aux1<<6;
        num_bits=num_bits+6;

        aux2=inBase64(&base64_ptr[i+j]);


        aux1=aux1|aux2;


        }else{
            aux1=aux1>>2;
            num_bits=num_bits-2;
        }

        }

        int count_bit=num_bits;
        while(num_bits!=0){
            num_bits=num_bits-4;
            aux3=(aux1>>num_bits) & 15;


           // printf("aux3:: %d ",aux3);
            hex[ik]=HEX_DIGITS[aux3];

          // printf("%c |",hex[ik]);
            ik++;
        }

       // printf("\n");


    }
    hex[ik]='\0';

    printf("length %d\n",ik);
   // printf("\n");
    int i=0;
    while(hex[i]!='\0'){
        printf("%c",hex[i]);
        i++;
    }

   // printf("\n");

}

char small_letter(char a){
   char uset[]="0123456789";
   for(int i=0;i<10;i++){
      if(uset[i]==a){
         return a;
      }
   }
   return a+32;
}
int Validating_File(char *file , key key){

    char tag1[30]="</content>";
    char tag2[30]="Signature";
    FILE *fptr;
    fptr=fopen(file,"r");
    char ch;
    char file_content[3000];//=(char*) malloc(sizeof(char)*3000);
    int i=0;

    while((ch=fgetc(fptr))!=EOF)
    {
        file_content[i]=ch;
        i++;
    }

    file_content[i]='\0';

    int index=isSubstring(tag1,file_content);
    index=index+int(strlen(tag1));



    unsigned char *st=new unsigned char[2000];
    unsigned int i_sha=0;
    //// assigning char to unsigned char, this has to be done to implement sha256
    int i_sha_a=0;
    while(i_sha_a<index+1)
    {

	   ch=file_content[i_sha_a];
      // printf("%c",ch);

       if(ch!='\n'){
        *(st+i_sha)=ch;
	//	printf("%c",*(st+i_sha));
		i_sha++;}
        i_sha_a++;

    }



    SHA256 sha;

	sha.update(st,i_sha);

	uint8_t * digest = sha.digest();
    int it=0;//just an iterator

    char HEX_format_Digest[65];

    // code for making hex string
    int hex_count=0;
    for(int yui=0;yui<32;yui++)
    {
        // printf("%d ",*(digest+yui)>>4||0x0f);
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui)>>4 & 0x0f];
         hex_count++;
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui) & 0x0f];
         hex_count++;
    }
    HEX_format_Digest[hex_count]='\0';
    static int pass=0;
    pass++;
    printf("\nhere is the     string %d: %s\n\n",pass,HEX_format_Digest);
   // delete []st;
    //st=NULL;

    char signature[500];//=(char*) malloc(sizeof(char)*500);

    fetch_tag(file_content,tag2,signature);
   // free(file_content);

    printf("signature\n%s\n",signature);

    // base64 decoder


   char hex[1000];//=(char*) malloc(sizeof(char)*(1000));

    base64decoder(&signature[0],hex);



   // printf("\n%s\n",hex);

    mp_int HEX1;
    mp_init(&HEX1);
    //mp_read_radix(&HEX1,"6676CB59FC89868EB6F2EF269CEF076E265C963779DE44B9E2E234A3391043B10E7667892400753214A9B1FD51AB7F48A429BD6AE73B0EC894785CCE3E0EFD735C4BBD54D2B9F7709629BC6C5A635F1AF52BBEBA1352D876154EDFA8EF4F1C58D4EFF9ADAB3EB81AF329B35595BA94B98505B67EBB814963B71C35312CA2904BA56CC2A4DDBD53D161BB900A74B5CF647531476A343293895433F70A0A35E7110EC220299A9F685BF6A98685925C3DA603BFC11EE0BF6E2216F47873DEF58EDB0CFB4CAE158F70E60E6233B09542CAA1F21722CAC5F24A8C09E4A32AFD34B879C8CA68E0DBD4CBE65F30A793333D5983B006EB91CC5FD86549939D526EBB3CE6",16);
    mp_read_radix(&HEX1,hex,16);
    //free(hex);
    //free(signature);
    mp_int modulus;
    mp_init(&modulus);
    mp_read_radix(&modulus,key.modulus,10);

    mp_int public_key;
    mp_init(&public_key);
    mp_read_radix(&public_key,"65537",10);

    mp_int decrypted;
    mp_init(&decrypted);
    mp_exptmod(&HEX1,&public_key,&modulus,&decrypted);

    char decrypted_hex[1000];//=(char*) malloc(sizeof(char)*1000);

    mp_to_hex(&decrypted,decrypted_hex,sizeof(decrypted_hex));

    printf("\n%s\n",decrypted_hex);//this is the encrypted message

    char useful_decrypted_hex[400];//=(char*) malloc(sizeof(char)*400);

    i=445;
    int i_counter=0;
    while(i<int(strlen(decrypted_hex))){
        useful_decrypted_hex[i_counter]=small_letter(decrypted_hex[i]);
        i++;
        i_counter++;
    }
    useful_decrypted_hex[i_counter]='\0';

    printf("heloo\n\n%s\n\n %s\n\n",useful_decrypted_hex,HEX_format_Digest);

  //  free(decrypted_hex);
   // free(useful_decrypted_hex);
   mp_clear(&HEX1);
   mp_clear(&modulus);
   mp_clear(&public_key);
   mp_clear(&decrypted);

    return 0;


}


// this function is for encrypting files that has to be remained inside the RFM
int encrypting_File(char *content, key key, char *fname)
{


int i=0;
char buf[50];
int aux;
mp_int public_key;
mp_init(&public_key);
mp_read_radix(&public_key,"65537",10);//
mp_int private_key;
mp_init(&private_key);
mp_read_radix(&private_key,key.private_exponent,10);//
mp_int modulus;
mp_init(&modulus);
mp_read_radix(&modulus,key.modulus,10);//

char encrypted_content[20000];//=(char*) malloc(sizeof(char)*20000);
mp_int aux_int;
mp_init(&aux_int);
mp_int aux_int_result;
mp_init(&aux_int_result);
char snum[10];
while(content[i]!='\0'){
    aux=int(content[i]);
    sprintf(snum, "%d",aux );
    mp_read_radix(&aux_int, snum,10);
  // printf("\n %c  %d \n",content[i],aux);
    mp_exptmod(&aux_int,&private_key,&modulus,&aux_int_result);

    mp_to_decimal(&aux_int_result,buf,sizeof(buf));
  //  printf("\n modulus product==\n%s\n\n",buf);
  //  break;

    if(i==0){
    strcpy(encrypted_content,buf);
    strcat(encrypted_content,";");
    }else{
        strcat(encrypted_content,buf);
        strcat(encrypted_content,";");
    }



    i=i+1;


}
strcat(encrypted_content,"\0");

FILE *fptr;

fptr=fopen(fname,"w");
fprintf(fptr,"%s",encrypted_content);
fclose(fptr);
mp_clear(&aux_int);
mp_clear(&aux_int_result);
mp_clear(&modulus);
mp_clear(&public_key);
mp_clear(&private_key);

}

// this function is for decrypting files that has to be remained inside the RFM
void decrypting_File(char *file , key key, char *result){
    FILE *fptr;
    fptr=fopen(file,"r");

    char buf[50];

    mp_int public_key;
    mp_init(&public_key);
    mp_read_radix(&public_key,"65537",10);//
    mp_int private_key;
    mp_init(&private_key);
    mp_read_radix(&private_key,key.private_exponent,10);//
    mp_int modulus;
    mp_init(&modulus);
    mp_read_radix(&modulus,key.modulus,10);//
    mp_int aux_int;
    mp_init(&aux_int);
    mp_int aux_int_result;
    mp_init(&aux_int_result);

    char tao;
    char decrypt;
    char buff2[30];
    int j=0,k=0,first_pass=0,di=0;

    while((tao=fgetc(fptr))!=EOF){

        if(tao!=';'){
            buff2[k]=tao;
            k++;
        }else{
            buff2[k]='\0';
            mp_read_radix(&aux_int, buff2,10);
            mp_exptmod(&aux_int,&public_key,&modulus,&aux_int_result);
            mp_to_decimal(&aux_int_result,buf,sizeof(buf));
      // printf("\n modulus product==\n%s\n\n",buf);
            decrypt = atoi(buf);
            result[di]=decrypt;
            di++;
            memset(buff2, 0, sizeof(buff2));
            k=0;
        }

    }
    result[di]='\0';
mp_clear(&aux_int);
mp_clear(&aux_int_result);
mp_clear(&modulus);
mp_clear(&public_key);
mp_clear(&private_key);

}



// this function is for fetching value of a particular tag in the given string

void fetch_tag(char *content, char *target, char *result){
    char TAG[30];
    strcpy(TAG,"<");
    strcat(TAG,target);

    int index=isSubstring(TAG,content);

    index=index+int(strlen(TAG))+1;
    int i=0;
    while(content[index]!='>'){
        result[i]=content[index];
        i++;
        index++;

    }
    result[i]='\0';


}


//function to generate DroneID.txt
void DroneIDcreation( ){

    char DroneID[]="ABBDDJEDNDJK";//any string
    char RFM_version[2]="0";
    char RPAS_category[10]="Small";


 // RFM_public_key: needs to be fetched from PublicPrivateInuse.txt(which is kept encrypted)
 // inside the rfm, decrypted using decrypting key.
    key RFM_private_key;
    char fname[60]="PublicPrivateInuse.txt";
    key Inside_RFM;
    strcpy(Inside_RFM.modulus,"180919775566931");
    strcpy(Inside_RFM.private_exponent,"32102716896161");

    char content[5000];//=(char*) malloc(sizeof(char)*5000);

    decrypting_File(fname , Inside_RFM, content);
   // printf("\n\nThe content of decrypted file is :\n\n%s\n",content);

    char target0[30]="Modulus";
    char target1[30]="PrivateKey";

    char value_modulus[2000];//=(char*) malloc(sizeof(char)*2000);
    char value_private[800];//=(char*) malloc(sizeof(char)*800);

    fetch_tag(content,target0,value_modulus);
    fetch_tag(content,target1,value_private);
    strcpy(RFM_private_key.private_exponent,value_private);
    strcpy(RFM_private_key.modulus,value_modulus);




    //printf("vadd %s",value_modulus);
   // free(content);

    char rfm_key0[30]="RFM_public_key_modulus";
    char rfm_key1[30]="RFM_public_key_exponent";



 //    this value is just for writing into the file
    char DigitalSky_public_key[1020]="MIIC8TCCAdmgAwIBAgIJAJRDnqfLydHvMA0GCSqGSIb3DQEBCwUAMA8xDTALBgNVBAMMBHRlc3QwHhcNMTkwMzI2MDcxMTQzWhcNMjkwMzIzMDcxMTQzWjAPMQ0wCwYDVQQDDAR0ZXN0MIIBIjANBgkqhkiG9w0BAQEFAAOCAQ8AMIIBCgKCAQEAq51cjR/mcgd0nWO33O3SM84yu3DRdaG8OMYSqzPixY5R+D8niOTVLZvOtaFROSneP1JmUAcaBn5sFhsFxgpJX8O6ee0m9PqLL+LKjexEs5dZ85IG8GqF+UJABaKfBeTPOgI5NAwoyZPBphzxsra1fH2OV2roaCf4ErMnYluuyey/VfFlHTVgC5+VX2wvO+o6pYUuzdNqCvgYwZrMEDCXm+08iZk/qpLgqgUCQTs8qGu/Y0d/EqwGmv9xN8tyxX+IbaeQM7uztN8PbMf8wY40OqdgNmgaVmMR4mfAO2XJiryR5Y8JACDGf3dhmcDrdtfmNjaHR109o2/wUPhSdWB/3QIDAQABo1AwTjAdBgNVHQ4EFgQUR4p2KJJXG5cZ8STI66RG6l2o7yowHwYDVR0jBBgwFoAUR4p2KJJXG5cZ8STI66RG6l2o7yowDAYDVR0TBAUwAwEB/zANBgkqhkiG9w0BAQsFAAOCAQEAV3uurlHMtyopefBpdGj59eLWCrpRYJLKbDtLFCj+tY1/uiwogUMNsEEHEBeEdwM+PIPuzWZ4tSYQ+SvdCCt4/6e9x+c2/1mZKhnRzL/s9o70RyWZXQO+Dz43B5aIIy/qARUhLxU2NVL42q90pInIh/ltT02IVkcibwDnsM4XJhsSyvQlRyYXdPzDeBjEOVYFpafLbC/7a5FBuNwfNKEMWhOj6AELnC8fWb3maNevhjSH5amGU2XrUp6yIdWUL2HuW7ReSer93Lg6iYujd/aaqk+pWE5bQsC+r2kHpNcpntHJLsd9E1cwzWCJiEM9zK4GXqKV/QDUdPC6FYfEf+ti9A==";

    char Firmware_version[4]="1.0";

    pair_set pairset[7];

    strcpy(pairset[0].tag,"DroneID");
    strcpy(pairset[0].value,DroneID);


    strcpy(pairset[1].tag,"RFM_version");
    strcpy(pairset[1].value,RFM_version);

    strcpy(pairset[2].tag,"RPAS_category");
    strcpy(pairset[2].value,RPAS_category);

    strcpy(pairset[3].tag,"Firmware_version");
    strcpy(pairset[3].value,Firmware_version);

    strcpy(pairset[4].tag,"DigitalSky_public_key");
    strcpy(pairset[4].value,DigitalSky_public_key);

    strcpy(pairset[5].tag,rfm_key0);
    strcpy(pairset[5].value,value_modulus);

    strcpy(pairset[6].tag,rfm_key1);
    strcpy(pairset[6].value,"65537");

    char fileName[20]="DroneID.txt";



  //  free(value_modulus);

    pair_file_write(pairset,7,fileName,RFM_private_key);





}

// this function  creates .txt file when an amendment takes
// place in the hardware(currently only gps)
void HardwareInuseCreation(int gps){
    char DroneID[]="ABBDDJEDNDJK";//any string
    char RFM_version[2]="0";
    char RPAS_category[10]="Small";
    char GPS_ID[20];
    sprintf(GPS_ID, "%d",gps);
 // RFM_public_key: needs to be fetched from PublicPrivateInuse.txt(which is kept encrypted)
 // inside the rfm, decrypted using decrypting key.
    key RFM_private_key;
    char fname[60]="PublicPrivateInuse.txt";
    key Inside_RFM;
    strcpy(Inside_RFM.modulus,"180919775566931");
    strcpy(Inside_RFM.private_exponent,"32102716896161");

    char content[5000];//=(char*) malloc(sizeof(char)*5000);

    decrypting_File(fname , Inside_RFM, content);
    //printf("\n\nThe content of decrypted file is :\n\n%s\n",content);

    char target0[30]="Modulus";
    char target1[30]="PrivateKey";

    char value_modulus[2000];//=(char*) malloc(sizeof(char)*2000);
    char value_private[800];//=(char*) malloc(sizeof(char)*800);

    fetch_tag(content,target0,value_modulus);
    fetch_tag(content,target1,value_private);
    strcpy(RFM_private_key.private_exponent,value_private);
    strcpy(RFM_private_key.modulus,value_modulus);
   /// RFM private key fetched



   // printf("vadd %s",value_modulus);
  //  free(content);




    pair_set pairset[4];

    strcpy(pairset[0].tag,"DroneID");
    strcpy(pairset[0].value,DroneID);


    strcpy(pairset[1].tag,"RFM_version");
    strcpy(pairset[1].value,RFM_version);

    strcpy(pairset[2].tag,"RPAS_category");
    strcpy(pairset[2].value,RPAS_category);



    strcpy(pairset[3].tag,"GPS_ID");
    strcpy(pairset[3].value,GPS_ID);



    char fileName[20]="HardwareInuse.txt";



   // free(value_modulus);

    pair_file_write(pairset,4,fileName,RFM_private_key);
}

void get_RFM_Key(key *key){
    char fname[60]="PublicPrivateInuse.txt";
    key Inside_RFM;
    strcpy(Inside_RFM.modulus,"180919775566931");
    strcpy(Inside_RFM.private_exponent,"32102716896161");

    char content[5000];//=(char*) malloc(sizeof(char)*5000);

    decrypting_File(fname , Inside_RFM, content);
   // printf("\n\nThe content of decrypted file is :\n\n%s\n",content);

    char target0[30]="Modulus";
    char target1[30]="PrivateKey";

    char value_modulus[2000];//=(char*) malloc(sizeof(char)*2000);
    char value_private[800];//=(char*) malloc(sizeof(char)*800);

    fetch_tag(content,target0,value_modulus);
    fetch_tag(content,target1,value_private);
    strcpy(key->private_exponent,value_private);
    strcpy(key->modulus,value_modulus);

}
int file_read(char *fname,char *tag,char *result,int file_type){

    key *RFM_key;
    key RFM_KEY;
    RFM_key=&RFM_KEY;

    get_RFM_Key(RFM_key);

    //file_type 0,1
    //0 = this is for file needing RFM key pair to get validated
    //1= this for file needing FMP key pair to get validated
    if (file_type==0){
        // RFM key pair needed to validate


    }

}
