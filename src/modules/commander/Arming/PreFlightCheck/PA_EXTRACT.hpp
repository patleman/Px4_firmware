#include<string.h>
#include<stdio.h>
#include "MP_INT.hpp"
#include "PreFlightCheck.hpp"
#include <drivers/drv_hrt.h>
#include <HealthFlags.h>
#include <lib/parameters/param.h>
#include <systemlib/mavlink_log.h>
#include <uORB/Subscription.hpp>
#include <uORB/topics/vehicle_gps_position.h>
#include<stdlib.h>
#include<time.h>

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

int aux_space_count;
struct stack{
    char list_att[40][40];
    int counter;
};

void updateStack_remove(char *ptr, stack *ss){
    int i=0;
    ss->counter=ss->counter-1;
    while(1){
    char t=ss->list_att[ss->counter][i];
    if (t=='\0'){
        break;
    }else{
        t=' ';
    }
    i=i+1;
    }

}

void updateStack_add(char *ptr,int set, stack *ss)
{
    int i=0;

    while(1){

        if(*(ptr+set+i+1)==' ' || *(ptr+set+i+1)=='>'){

            break;
        }

        ss->list_att[ss->counter][i]=ptr[set+i+1];
        i++;
    }
    ss->list_att[ss->counter][i]='\0';
   // printf("ha %s ha",ss->list_att[ss->counter]);
    ss->counter=ss->counter+1;


}

void replaceTag(char *ptr, stack *ss,char *result,int *ptr_count){
    char aux[40];
    int i=0;
    ss->counter=ss->counter-1;
    while(1){
        char aux_v=ss->list_att[ss->counter][i];
        if (aux_v=='\0' || aux_v==' '){
            if(aux_v=='\0'){
               // printf("endarray");
            }
            if(aux_v==' '){
               // printf("space");
            }
            break;
        }
        aux[i]=aux_v;
        i++;
    }

  //  printf("%c",'>');
    result[*(ptr_count)]='>';
    *(ptr_count)=*(ptr_count)+1;
  //  printf("%c",'<');
    result[*(ptr_count)]='<';
    *(ptr_count)=*(ptr_count)+1;

   // printf("%c",'/');
    result[*(ptr_count)]='/';
    *(ptr_count)=*(ptr_count)+1;


    for(int f=0;f<i;f++){

       // printf("%c",aux[f]);
        result[*(ptr_count)]=aux[f];
        *(ptr_count)=*(ptr_count)+1;

    }
}

int isSignature(char *ptr){
    /* if signature start ::0
       if signature end ::1
       if none  ::2
    */
    char check0[18]="<Signature xmlns=";

    char check1[14]="</Signature>";

    char check2[14]="</SignedInfo>";

    char check3[13]="<SignedInfo>";

    int check0_status=1;

    for(int i=0;i<(int)strlen(check0);i++)
    {
        if(*(ptr+i)!=check0[i])
        {
            check0_status=0;
            break;
        }
    }

    int check1_status=1;

    for(int i=0;i<(int)strlen(check1);i++)
    {
        if(*(ptr+i)!=check1[i])
        {
            check1_status=0;
            break;
        }
    }

    int check2_status=1;
    for(int i=0;i<(int)strlen(check2);i++){
        if(*(ptr+i)!=check2[i]){
            check2_status=0;
            break;
        }
    }


    int check3_status=1;
    for(int i=0;i<(int)strlen(check3);i++){
        if(*(ptr+i)!=check3[i]){
            check3_status=0;
            break;
        }
    }

    if(check0_status==1){
        return 0;
    }
    if(check1_status==1){
        return 1;
    }
    if(check2_status==1){
        return 2;
    }
    if(check3_status==1){
        return 3;
    }
    return 4;

}

void Reference_canon(char *file_name, char *res){

    //char space=' ';
    //char new_line='\n';

    FILE *fp;

    //char open_bracket='<';//  0
    //char close_bracket='>';//  1
    //char slash_bracket[3]="/>";//  2

    fp=fopen(file_name,"r");
    if(fp)
    {


	    //exit(1);

    char ch;
    int count_result=0;
    int *ptr_count_result;
    ptr_count_result=&count_result;
    int i=0;



    char *content=(char*) malloc(sizeof(char)*10000);

    while(1)
    {

        ch=fgetc(fp);
        if(ch==EOF){
            break;
        }

        content[i]=ch;
        i++;

    }
    stack *dd;
    stack dds;
    dds.counter=0;
    dd=&dds;


    int flag0=0,flag1=0;
    int aux_check_flag;
    int space_count=0;
    for(int u=0;u<i;u++)
    {

        if(content[u]=='<')
        {
            aux_check_flag=isSignature(content+u);
            if(aux_check_flag==0)
                {   aux_space_count=space_count;
                    flag0=1;
                }
            if(aux_check_flag==1){
                flag1=1;
                aux_space_count=space_count;
                u=u+12;

            }
            if(aux_check_flag==2){

            }
        }

        if((flag0==0 && flag1==0) || (flag0==1 && flag1==1) )
        {

            if(content[u]=='<')
            {
                if(content[u+1]!='/')
                {
                    updateStack_add(content,u,dd);
                    space_count=space_count+2;
                }

            }

            if(content[u]=='/' && content[u+1]=='>')
                {

                    replaceTag(content+u,dd,res,ptr_count_result);
                    space_count=space_count-2;

                }else{


                res[count_result]=content[u];
                count_result=count_result+1;
                }

            if(content[u]=='>' && content[u+1]=='<')
            {


                if(content[u+2]=='/')
                {
                    updateStack_remove(content+u,dd);

                     space_count=space_count-2;
                }





            }

        }

    }
    res[count_result]='\0';

    free(content);}
    else{
        printf("PA File does not exist");
    }
}


void SignedInfo_canon(char *file_name, char *res)
{

    FILE *fp;

    fp=fopen(file_name,"r");
    if(fp)
    {


    char *content=(char*) malloc(sizeof(char)*10000);

    char ch;
    char aux_copy[300];
    int count=0;
    int flag_ini=0;
    int flag_end=0;
     int *ptr_count_result;

    int res_count=0;

    ptr_count_result=&res_count;

    while(1)
    {

        ch=fgetc(fp);
        if(ch==EOF){
            break;
        }
        content[count]=ch;
        count++;
    }

    stack dde;
    stack *de;
    de=&dde;
    dde.counter=0;
    int space_count;
    int first_pass_flag=0;
    int in_count=0;//aux length
    int first_time_flag=0;
    for(int i=0;i<count;i++)

    {
        if(content[i]=='<')
        {

        int check_Signature;

        check_Signature=isSignature(content+i);
        if(check_Signature==0)
        {  // printf("\nfound start of signed info\n");
            int k=11;
            while(content[k+i]!='>')
            {
                aux_copy[k-11]=content[k+i];
                k++;
                in_count++;

            }
            aux_copy[k]='\0';
          //  printf("%s\n",aux_copy);
            flag_ini=1;
        }

        if(check_Signature==2)
        {
            //end()
           // printf("\njhghggg\n");
            flag_end=1;
            for(int f=0;f<13;f++)
            {
            res[res_count]=content[i+f];
            res_count++;
            }
        }

        if(check_Signature==3)
        {
           // content[i+12]=' ';
            space_count=aux_space_count;
            for(int o=0;o<11;o++)
            {
                res[res_count]=content[i+o];
              //  printf("%c",res[res_count]);
                res_count++;
            }
            i=i+12;
            res[res_count]=' ';
            res_count++;

            for(int j=0;j<in_count;j++)
            {
                res[res_count]=aux_copy[j];
              //  printf("%c",res[res_count]);
                res_count++;
            }
            res[res_count]='>';
            res_count++;
          //  res[res_count]='\n';
          //  res_count++;

            first_pass_flag=1;
        }
        }

        if(flag_ini==1 && first_pass_flag==1 && flag_end!=1)
        {
            if(first_time_flag==0){
                first_time_flag=1;
                space_count=space_count+4;
                for(int y=0;y<space_count;y++)
                {
                   // res[res_count]=space;
                   // res_count++;
                }
            }

            if(content[i]=='<')
            {
                if(content[i+1]!='/')
                {
                    updateStack_add(content,i,de);
                    space_count=space_count+2;
                }
                if(content[i+1]=='/' && content[i-1]!='>')
                {

                    updateStack_remove(content+i,de);
                   /*  printf("\nprinting stack::\n");
                    for(int i=0;i<de->counter;i++){
                        printf("%s ",de->list_att[i]);
                    }
                    printf("\n\n");*/
                    space_count=space_count-2;
                }
              //  printf("\npassed updatestack\n");

            }

            if(content[i]=='/' && content[i+1]=='>')
                {

                    replaceTag(content+i,de,res,ptr_count_result);
                    space_count=space_count-2;

                }else{

              //  printf("%c",content[i]);
                res[res_count]=content[i];
                res_count=res_count+1;
                }

            if(content[i]=='>' && content[i+1]=='<')
            {
               // printf("\n");//new line

               // res[res_count]=0x0a;//'\n';

               // res_count=res_count+1;


                if(content[i+2]=='/')
                {
                    updateStack_remove(content+i,de);
                  /*  printf("\nprinting stack::\n");
                    for(int i=0;i<de->counter;i++){
                        printf("%s ",de->list_att[i]);
                    }
                    printf("\n\n");*/
                     space_count=space_count-2;
                }





            }


        }

    }
    free(content);
    fclose(fp);
    }
    else{
        printf("PA file does not exist");

    }



}

//function get tag value
void getTagvalue(char *Tag, char *certificate,char *file_nam_perm){
   FILE *fp;
   fp=fopen(file_nam_perm,"r");
   // Digest DD;

    char buf[3000];

    char tagi[30];
    char tage[30];

    memset(tagi, 0, sizeof(tagi));
    strcpy(tagi, "<");
    strcat(tagi,Tag );
    strcat(tagi, ">");

    memset(tage, 0, sizeof(tage));
    strcpy(tage, "</");
    strcat(tage,Tag );
    strcat(tage, ">");

    char *ptr1,*ptr2;

    int c_flag_i=0,c_flag_e=0,line=0;

    int index[2][2]={{0,0},{0,0}};

    int  aux=0;

    int c=0;

   while(fscanf(fp, "%s", buf) != EOF )
    {
    line=line+1;
    if(isSubstring(tagi,buf)!=-1 || isSubstring(tage,buf)!=-1 ){
      //  printf(" %d ",line);
        if (isSubstring(tagi,buf)!=-1){       //for "<X509Certificate>"
          // printf("start");
           index[0][0]=line;
           index[0][1]=isSubstring(tagi,buf);
           c_flag_i=1;
           aux=1;

        }
        if(isSubstring(tage,buf)!=-1){                          //for "</X509Certificate>"
          // printf("  end");
           index[1][0]=line;
           index[1][1]=isSubstring(tage,buf);
          // printf("\n%d   %d\n",index[1][0],index[1][1]);
           c_flag_e=1;


        }
        }

        if(c_flag_i==1 && c_flag_e==0 ){
        if (aux==1){

        for(int u=index[0][1]+strlen(tagi);u<(int)strlen(buf);u++){
            certificate[c]=buf[u];

            c++;
        }

        aux=0;
        }
        else{ptr1=certificate+c;
            ptr2=&(buf[0]);
            memcpy(ptr1,ptr2,strlen(buf));
            c=c+strlen(buf);
        }
        }
        else if(c_flag_i==1 && c_flag_e==1 )
        {
            if(index[0][0]==index[1][0])//start line and end line is same
        {
            for(int u=index[0][1]+strlen(tagi);u<index[1][1];u++)
            {
            certificate[c]=buf[u];
            c++;
            }

        }else
        {  //end line is different than start line
           for(int u=0;u<index[1][1];u++)
            {
            certificate[c]=buf[u];
            c++;
            }

        }
        break;
        }else{continue;}

    }
    fclose(fp);

}
//static const char *const BASE64_DIGITS ="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

//static const char *const HEX_DIGITS = "0123456789abcdef";
int inBase64(char *d)
  {
        for(int i=0;i<64;i++){
            if(*d==BASE64_DIGITS[i]){
                return i;
            }
        }
	return 0;
  }

void base64decoder(char *base64_ptr, char *result){

   int length_str=0;
   while(*base64_ptr!='\0'){
            length_str++;
            base64_ptr++;
      }
   base64_ptr=base64_ptr-length_str;

      int size_of_hex;
      if((length_str*6)%8!=0){
            printf("invalid base64 to hex\n");
      }
      size_of_hex=(length_str*6)/8;

      char *hex=(char*) malloc(sizeof(char)*(2*size_of_hex+100));

      int ik=0;

      int aux1,aux2,aux3,num_bits;

      for(int i=0;i<length_str;i=i+4){
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

        // int count_bit=num_bits;
         while(num_bits!=0){
               num_bits=num_bits-4;
               aux3=(aux1>>num_bits) & 15;


           // printf("aux3:: %d ",aux3);
         hex[ik]=HEX_DIGITS[aux3];
         *(result+ik)=hex[ik];
         printf("%c",hex[ik]);
         ik++;
         }
      }
   result[ik]='\0';
   free(hex);
   printf("\n");
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



char* strip(char *str)
        {
        size_t len = strlen(str);
        memmove(str, str+1, len-2);
        str[len-2] = 0;
        return str;
        }




int Is_PA_VAlid(){
	char Tag_Digest_value0[12]="DigestValue";
        char Digest_Value0[100];// in base 64
        char Digest_Value_in_hex[100];// in hex

        char Tag_Signed_Value0[17]="SignatureValue";
        char Signature_Value0[550];// in base 64
        char Signature_Value_in_hex[500];// in hex
        char  Extracted_sha_from_Signature_value[100];

        char file_name[48]="./log/permission_artifact_breach.xml";//"permission_artifact_1.xml";



        char Sha_of_Reference[300];

            char Sha_of_SignedInfo[300];

            char SignedInfo_canonilized[2000];

            char Reference_canonilized[2000];// variable to hold canonicalized reference section


        // This function will extract the reference section out of the xml file (PA)
	    //and willl perform canonicalization step over it
        Reference_canon(file_name,Reference_canonilized);

	    //This section will extract the signed info section out of the xml file (PA)
	    //and will perform canonicalization step over it
	    SignedInfo_canon(file_name,SignedInfo_canonilized);

        printf("\n");


	    // SHA value calculated for extracted canonicalized SIgned info and reference section
	    //start
        char ch;
        unsigned char st_arr[3000];
	    unsigned char *st=&st_arr[0];
        unsigned int iff=0;
        // assigning char to unsigned char (this is done to have asci value in 8 bits rather than 7 bits)
        while(1){
	      ch=Reference_canonilized[iff];
	     if(ch=='\0'){
		   break;
	     }
         *(st+iff)=ch;
		printf("%c",*(st+iff));
		iff++;
        }


	printf("\n");
	printf("\n");

        SHA256 sha;

	sha.update(st,iff);

	uint8_t * digest = sha.digest();
        int it=0;//just an iterator
	//printing digest
	printf("\n");

        while(it<32){
		printf("%02x",*(digest+it));
		it++;
	}
       int aux0;
       int count_digest=0;
       it=0;
       while(it<32)
          {

            aux0=(*(digest+it)<<24)|(*(digest+it+1)<<16)|(*(digest+it+2)<<8)|*(digest+it+3);

            for(int h=28;h>=0;h=h-4){
               Sha_of_Reference[count_digest] =HEX_DIGITS[(aux0>>h)&0xf];
               count_digest++;
            }

		    it=it+4;
	   }

      Sha_of_Reference[count_digest]='\0';
      printf("\nSha of reference section that has been canonicalized :%s\n",Sha_of_Reference);

      printf("\n");


   char ch1;
   unsigned char st1_arr[3000];
   unsigned char *st1=&st1_arr[0];
   unsigned int i1=0;
    //// assigning char to unsigned char
    while(1){
	   ch1=SignedInfo_canonilized[i1];
	   if(ch1=='\0'){
		   break;
	   }
        *(st1+i1)=ch1;
		printf("%c",*(st1+i1));
		i1++;

    }


	printf("\n");
	printf("\n");

   SHA256 sha1;

	sha1.update(st1,i1);

	uint8_t * digest1 = sha1.digest();
   it=0;//just an iterator
	//printing digest
	printf("\n");

   while(it<32){
		printf("%02x",*(digest1+it));
		it++;
	}
   int aux1;
   count_digest=0;
   it=0;
   while(it<32)
      {

            aux1=(*(digest1+it)<<24)|(*(digest1+it+1)<<16)|(*(digest1+it+2)<<8)|*(digest1+it+3);

            for(int h=28;h>=0;h=h-4){
               Sha_of_SignedInfo[count_digest] =HEX_DIGITS[(aux1>>h)&0xf];
               count_digest++;
            }

		    it=it+4;
	   }

   Sha_of_SignedInfo[count_digest]='\0';
   printf("\nsha of SignedInfo section that has been canonicalized :%s\n",Sha_of_SignedInfo);

	printf("\n");


   // SHA value calculated for extracted canonicalized SIgned info and reference section
   // end

   //In coming functions :
   //1) extracting digest value
   //2)converting it from base64 to hex
   //same above two steps for Signature value.

    // start
    // storing Digest value
   getTagvalue(Tag_Digest_value0,Digest_Value0,file_name);
   //  printf("\nhello  %d\n",Digestv.len);


   printf("\n Digest value in the permission artefact is %s\n",Digest_Value0);

   printf("\n");
   base64decoder(Digest_Value0, Digest_Value_in_hex);
   printf("\n");
   printf("\n Digest value in hex format %s \n",Digest_Value_in_hex);


   // storing Signed value
   getTagvalue(Tag_Signed_Value0,Signature_Value0, file_name);



   printf("\n Signature value in the permission artefact is %s\n",Signature_Value0);
   printf("\n");

   base64decoder(Signature_Value0, Signature_Value_in_hex);
   printf("\n");
   printf("\n Signature value in hex format %s \n",Signature_Value_in_hex);

   //end

   // message: bignum object created for storing value of Signature value
   // now Performing modulus operation on bignumber integers

   mp_int message,modulus,public_key,Decrypted;
   mp_init_multi(&message,&modulus,&public_key,&Decrypted,NULL);
   mp_read_radix(&message,Signature_Value_in_hex,16);

   /// Public key of dgca: This has to be taken from a reserved file in directory./ firmware update

 //  mp_int ;
   char Modulus[513]="ab9d5c8d1fe67207749d63b7dcedd233ce32bb70d175a1bc38c612ab33e2c58e51f83f2788e4d52d9bceb5a1513929de3f526650071a067e6c161b05c60a495fc3ba79ed26f4fa8b2fe2ca8dec44b39759f39206f06a85f9424005a29f05e4cf3a0239340c28c993c1a61cf1b2b6b57c7d8e576ae86827f812b327625baec9ecbf55f1651d35600b9f955f6c2f3bea3aa5852ecdd36a0af818c19acc1030979bed3c89993faa92e0aa0502413b3ca86bbf63477f12ac069aff7137cb72c57f886da79033bbb3b4df0f6cc7fcc18e343aa76036681a566311e267c03b65c98abc91e58f090020c67f776199c0eb76d7e6363687475d3da36ff050f85275607fdd";
   mp_read_radix(&modulus,Modulus,16);

  // mp_int ;
   mp_read_radix(&public_key,"65537",10);

   //mp_int ;
   mp_exptmod(&message, &public_key,&modulus,&Decrypted);   //                       this part for decrypting encrypted text

   char message_string_hex[513];
   mp_to_hex(&Decrypted,message_string_hex,sizeof(message_string_hex));
   printf("\n Decrypted signature value : %s\n",message_string_hex);//this is the decrypted message
   // At this point message has been decrypted
   // Now having decrypted message in small letters
   mp_clear_multi(&message,&modulus,&public_key,&Decrypted,NULL);

   printf("\n\n");
   int k_cout=0;
   for(int i2=445;message_string_hex[i2]!='\0';i2++){
      Extracted_sha_from_Signature_value[k_cout]=small_letter(message_string_hex[i2]);
      k_cout++;
     // printf("%c",message_string_hex[i]);
   }
   //with Extracted SHA we are reffering to the expected SHA of SIgned info section
   Extracted_sha_from_Signature_value[k_cout]='\0';//
   printf("\n%s\n",Extracted_sha_from_Signature_value);

   if(strcmp(Extracted_sha_from_Signature_value,Sha_of_SignedInfo)==0){
      printf("SEE if next two vstring strings are same");
      printf("\n%s\n ",Sha_of_Reference);
      printf("%s\n",Digest_Value_in_hex);
      if(strcmp(Sha_of_Reference,Digest_Value_in_hex)==0){
         printf("\n PA is valid and non tampered \n");
	  int pa_done=1;
          param_set(param_find("PA_CHECK"),&pa_done);
	 return 1;
      }else{
         printf("\n PA is not valid\n");
         return 0;
      }

   }else{
      printf("\n PA is not valid \n");
    //  PX4_INFO("PA is valid???????????");
      return 0;
   }

/// at this point validation is complete

}


void Key_rotation_start(){
  /// NOW prime number generation begins for KEY pair creation (later on this can be added in separate module also)
  // FILE *key_rotation;
  // key_rotation= fopen("Key_rotation_required","r");

int key_rotation_buff = 0;
param_set(param_find("KEY_ROT"),&key_rotation_buff);


   //   fclose(key_rotation);

   mp_int  p1, q1;

   char    buf[4096];
   //mp_int z,r;
   mp_init(&p1);
   mp_init(&q1);

  // mp_prime_rand(&p1,1,90,MP_PRIME_SAFE);
   mp_prime_rand(&p1,1,1024,MP_PRIME_2MSB_ON);
   mp_to_decimal(&p1,buf,sizeof(buf));
   printf("\nprime1 == %s \n",buf);


   //mp_prime_rand(&q1,1,90,MP_PRIME_SAFE);
   mp_prime_rand(&q1,1,1024,MP_PRIME_2MSB_ON);
   mp_to_decimal(&q1,buf,sizeof(buf));
   printf("\n prime2 == %s \n",buf);

   //phi(n)=(p1-1)*(q1-1)
   mp_int s1,s2,product_p1_q1;
   mp_init_multi(&s1,&s2,&product_p1_q1);

   mp_mul(&p1,&q1,&product_p1_q1);//modulus n=p1*q1

   mp_sub_d(&p1, 1uL, &s1);
   mp_to_decimal(&s1, buf, sizeof(buf));
   printf("\n\ns1 == %s\n", buf);//p1-1

   mp_sub_d(&q1, 1uL, &s2);
   mp_to_decimal(&s2, buf, sizeof(buf));
   printf("\n\ns2 == %s\n", buf);//q1-1

   mp_int product_s1_s2,e,d;
   mp_init_multi(&e,&d,&product_s1_s2);
   mp_mul(&s1,&s2,&product_s1_s2);// phi(n)=(p1-1)*(q1-1)
   //mp_int e;
   mp_read_radix(&e,"65537",10);//
 //  mp_int d;
   mp_invmod(&e,&product_s1_s2,&d);

   mp_to_decimal(&e,buf,sizeof(buf));
   printf("\n Public key : e===\n%s\n\n",buf);
   mp_to_decimal(&d,buf,sizeof(buf));
   printf("\n Private key : d==\n%s\n\n",buf);
   mp_to_decimal(&p1,buf,sizeof(buf));
   printf("\n prime number 1: p1==\n%s\n\n",buf);
   mp_to_decimal(&q1,buf,sizeof(buf));
   printf("\n prime number 2: q1==\n%s\n\n",buf);


   mp_to_decimal(&product_p1_q1,buf,sizeof(buf));
   printf("\n modulus product==\n%s\n\n",buf);
   int key_rot_done=1;
   param_set(param_find("KEY_ROT"),&key_rot_done);
   mp_clear(&d);
   mp_clear(&p1);
   mp_clear(&q1);
   mp_clear(&e);
   mp_clear(&product_s1_s2);
   mp_clear(&s1);
   mp_clear(&s2);
   mp_clear(&product_p1_q1);





/// KEY pair generation ends here
}

struct Date_time{
       int Hours;
       int Minutes;
       int Seconds;
       int date;
       int Month;
       int Year;
};
struct GEO_DATE_TIME_XML{// for storing xml data
    Date_time start;
    Date_time end;
};

struct ParsedData{
   char start_time[50];
   char end_time[50];
   float long_lat_coords[20][10];
   int num_coords;

};

ParsedData parse_artifact()
{


   ParsedData result{};

   FILE *fp;
   fp = fopen("./log/permission_artifact_breach.xml", "r"); // read mode

   char buf[1000], start_time[200], end_time[200],long_lat_coords[20][20];
   rewind(fp);

   int coord_ind = 0;

   while(fscanf(fp, "%s", buf) != EOF )
		{
         int succ = 0;
         succ = sscanf(buf,"flightEndTime=%s",buf);
			if (succ != 0)strcpy(end_time,strip(buf));
         succ = sscanf(buf,"flightStartTime=%s",buf);
			if (succ != 0) strcpy(start_time,strip(buf));
         succ = sscanf(buf,"latitude=%s>",buf);
         if (succ != 0) {
            strcpy(long_lat_coords[coord_ind],strip(buf));
            coord_ind++;
         }
         succ = sscanf(buf,"longitude=%[^/>]s>",buf);
         if (succ != 0) {
            strcpy(long_lat_coords[coord_ind],strip(buf));
            coord_ind++;
         }
		}

      strcpy(result.start_time,start_time);
      strcpy(result.end_time , end_time);
      result.num_coords = static_cast<int>(coord_ind/2);
      for(int i  = 0; i <coord_ind;i+=2){
        result.long_lat_coords[i][0] = atof(long_lat_coords[i]);
        result.long_lat_coords[i][1] = atof(long_lat_coords[i+1]);
      }
   fclose(fp);
   return result;
}


bool In_Time(Date_time current, GEO_DATE_TIME_XML Xml){

	// full check
    printf("hello hello");
	if(Xml.start.Year<=current.Year && current.Year<=Xml.end.Year){
        printf("year is ok");
		if(Xml.start.Month<=current.Month && current.Month<=Xml.end.Month){
            printf("month is ok");
                       if(Xml.start.date<=current.date && current.date<=Xml.end.date){
                           printf("date is ok");

		       }
		       else{ printf("%d",current.Year);
			       return 0;
		       }
		}else{
			return 0;
		}
	}
	else{
        	return 0;
	}
	//time check
	if(Xml.start.Hours<=current.Hours||current.Hours<=Xml.end.Hours){
        printf("\nhours is ok");
        if(Xml.start.Hours==current.Hours||current.Hours==Xml.end.Hours){
            printf("both the hours are equal");
            if(Xml.start.Hours==current.Hours){
                printf("\nstart hours are equal");
                if(Xml.start.Minutes<=current.Minutes){
                    printf(" minutes are ok");
                    if(Xml.start.Minutes==current.Minutes){
                        printf("miutes are equal");
                    if(Xml.start.Seconds<=current.Seconds){
                        return 1;
                    }else{return 0;}
                }else{return 1;}
                }else{return 0;}
            }
            else{
                 if(Xml.end.Minutes>=current.Minutes){
                    if(Xml.end.Minutes==current.Minutes){
                    if(Xml.end.Seconds>current.Seconds){
                        return 1;
                    }else{return 0;}
                }else{return 1;}
                }else{return 0;}

            }
        }
}else{return 0;}
return 0;
}


int date_time_extract_and_check() {

	FILE * file;
	int vehi_gps_pos=orb_subscribe(ORB_ID(vehicle_gps_position));
        time_t timestamp;
        struct vehicle_gps_position_s raw;
        orb_copy(ORB_ID(vehicle_gps_position),vehi_gps_pos,&raw);
        timestamp=(raw.time_utc_usec)/1000000;//microsecond to seconds
	printf("\n%u\n",timestamp);
        struct tm  ts;
        Date_time DT;//current time
        ts = *localtime(&timestamp);
        DT.Year=ts.tm_year+1900;
        DT.Month=ts.tm_mon+1;
        DT.date=ts.tm_mday;
	DT.Hours=ts.tm_hour;
	DT.Minutes=ts.tm_min;
	DT.Seconds=ts.tm_sec;
	printf("\n%d/%d/%d current %d:%d:%d \n",DT.Year,DT.Month,DT.date,DT.Hours,DT.Minutes,DT.Seconds);

	file = fopen("./log/permission_artifact_breach.xml", "r");// ./log/ for posix /log/ for nuttx
	if (file){
	fclose(file);
	GEO_DATE_TIME_XML Rv;
	ParsedData data= parse_artifact();
        printf("\n  %s          %s \n ",data.start_time,data.end_time);
        char aux_buf1[60];
        // extracting start year
	printf("\n");
	int length_str=0;
	int i_iter=0;
	while(length_str<(int)strlen(data.start_time)){
		if (data.start_time[i_iter]=='.'){
         	break;
                }

      		if ( data.start_time[i_iter]!='T' && data.start_time[i_iter]!='-' && data.start_time[i_iter]!=':'){


      		aux_buf1[length_str]=data.start_time[i_iter];

      		if(i_iter==3){
         	aux_buf1[length_str+1]='\0';
           	printf(" %s ",aux_buf1);
         	Rv.start.Year=strtol(aux_buf1,NULL,10);
         	printf("year: %d ",Rv.start.Year);
         	memset(aux_buf1,'0',4);

      		}
      if(i_iter==6){
         aux_buf1[length_str+1]='\0';
           printf(" %s ",aux_buf1);
         Rv.start.Month=strtol(aux_buf1,NULL,10);
         printf("month: %d ",Rv.start.Month);
         memset(aux_buf1,'0',7);
      }
       if(i_iter==9){
         aux_buf1[length_str+1]='\0';
         printf(" %s ",aux_buf1);
         Rv.start.date=strtol(aux_buf1,NULL,10);
         printf("date: %d ",Rv.start.date);
         memset(aux_buf1,'0',10);
      }
       if(i_iter==12){
         aux_buf1[length_str+1]='\0';
         printf(" %s ",aux_buf1);
         Rv.start.Hours=strtol(aux_buf1,NULL,10);
         printf("Hours: %d ",Rv.start.Hours);
         memset(aux_buf1,'0',10);
      }
      if(i_iter==15){
         aux_buf1[length_str+1]='\0';
         printf(" %s ",aux_buf1);
         Rv.start.Minutes=strtol(aux_buf1,NULL,10);
         printf("Minutes: %d ",Rv.start.Minutes);
         memset(aux_buf1,'0',12);
      }
      if(i_iter==18){
         aux_buf1[length_str+1]='\0';
         printf(" %s ",aux_buf1);
         Rv.start.Seconds=strtol(aux_buf1,NULL,10);
         printf("Seconds: %d ",Rv.start.Seconds);
         memset(aux_buf1,'0',12);
      }

        length_str++;
      }
      i_iter++;

   }
   printf("\n\n");
   memset(aux_buf1,'0',60);
   // extracting end year
   printf("\n");
   length_str=0;
   i_iter=0;
   while(length_str<(int)strlen(data.end_time)){
      if (data.end_time[i_iter]=='.'){
         break;
      }

      if ( data.end_time[i_iter]!='T' && data.end_time[i_iter]!='-' && data.end_time[i_iter]!=':'){


      aux_buf1[length_str]=data.end_time[i_iter];

      if(i_iter==3){
         aux_buf1[length_str+1]='\0';
           printf(" %s ",aux_buf1);
         Rv.end.Year=strtol(aux_buf1,NULL,10);
         printf("year: %d ",Rv.end.Year);
         memset(aux_buf1,'0',4);

      }
      if(i_iter==6){
         aux_buf1[length_str+1]='\0';
           printf(" %s ",aux_buf1);
         Rv.end.Month=strtol(aux_buf1,NULL,10);
         printf("month: %d ",Rv.end.Month);
         memset(aux_buf1,'0',7);
      }
       if(i_iter==9){
         aux_buf1[length_str+1]='\0';
         printf(" %s ",aux_buf1);
         Rv.end.date=strtol(aux_buf1,NULL,10);
         printf("date: %d ",Rv.start.date);
         memset(aux_buf1,'0',10);
      }
       if(i_iter==12){
         aux_buf1[length_str+1]='\0';
         printf(" %s ",aux_buf1);
         Rv.end.Hours=strtol(aux_buf1,NULL,10);
         printf("Hours: %d ",Rv.end.Hours);
         memset(aux_buf1,'0',10);
      }
      if(i_iter==15){
         aux_buf1[length_str+1]='\0';
         printf(" %s ",aux_buf1);
         Rv.end.Minutes=strtol(aux_buf1,NULL,10);
         printf("Minutes: %d ",Rv.end.Minutes);
         memset(aux_buf1,'0',12);
      }
      if(i_iter==18){
         aux_buf1[length_str+1]='\0';
         printf(" %s ",aux_buf1);
         Rv.end.Seconds=strtol(aux_buf1,NULL,10);
         printf("Seconds: %d ",Rv.end.Seconds);
         memset(aux_buf1,'0',12);
      }

        length_str++;
      }
      i_iter++;

   }
   printf("\n\n");

   printf("\nfrom xml:: %d/%d/%d and  %d:%d:%d",Rv.end.Year,Rv.end.Month,Rv.end.date,Rv.end.Hours,Rv.end.Minutes,Rv.end.Seconds);
   printf("\n\n");
   printf("\nfrom xml:: %d/%d/%d and  %d:%d:%d",Rv.start.Year,Rv.start.Month,Rv.start.date,Rv.start.Hours,Rv.start.Minutes,Rv.start.Seconds);


	//int has_artifact = 1;
	//param_set(param_find("PERM_ARTIFACT"),&has_artifact);
	printf("Permission Artifact Found \n");
	bool in_time=In_Time(DT,Rv);
	if (in_time!=1){
		printf("\n permission denied\n");
		return 0;
	}else{
		return 1;
	}

	}

	else{
   	//file doesn't exists or cannot be opened (es. you don't have access permission)
	//int has_artifact = 0;
	//param_set(param_find("PERM_ARTIFACT"),&has_artifact);
	return 0;

	}
}
