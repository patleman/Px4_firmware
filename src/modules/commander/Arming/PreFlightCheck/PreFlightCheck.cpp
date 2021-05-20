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

#include <drivers/drv_hrt.h>
#include <HealthFlags.h>
#include <lib/parameters/param.h>
#include <systemlib/mavlink_log.h>
#include <uORB/Subscription.hpp>


using namespace time_literals;

/*
Class for SHA256 implementation

///////////////////////////////////////////

static const char *const BASE64_DIGITS ="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static const char *const HEX_DIGITS = "0123456789abcdef";

class SHA256 {

public:
	SHA256();
	void update(uint8_t * data, size_t length);
	uint8_t * digest();


private:
	uint8_t  m_data[64];
	uint32_t m_blocklen;
	uint64_t m_bitlen;
	uint32_t m_state[8]; //A, B, C, D, E, F, G, H

    uint32_t K[64] = {
		0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,
		0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
		0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
		0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
		0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,
		0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
		0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,
		0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
		0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
		0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
		0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,
		0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
		0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,
		0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
		0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
		0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
	};

	static uint32_t rotr(uint32_t x, uint32_t n);
	static uint32_t choose(uint32_t e, uint32_t f, uint32_t g);
	static uint32_t majority(uint32_t a, uint32_t b, uint32_t c);
	static uint32_t sig0(uint32_t x);
	static uint32_t sig1(uint32_t x);
	void transform();
	void pad();
	void revert(uint8_t * hash);

};

SHA256::SHA256(): m_blocklen(0), m_bitlen(0) {
	m_state[0] = 0x6a09e667;
	m_state[1] = 0xbb67ae85;
	m_state[2] = 0x3c6ef372;
	m_state[3] = 0xa54ff53a;
	m_state[4] = 0x510e527f;
	m_state[5] = 0x9b05688c;
	m_state[6] = 0x1f83d9ab;
	m_state[7] = 0x5be0cd19;
}
void SHA256::update( uint8_t * data, size_t length) {
	for (size_t i = 0 ; i < length ; i++) {
		m_data[m_blocklen++] = data[i];
	   // printf("%c",data[i]);
		if (m_blocklen == 64) {
			transform();

			// End of the block
			m_bitlen += 512;
			m_blocklen = 0;
		}
	}
}

uint8_t * SHA256::digest() {
	uint8_t * hash = new uint8_t[32];
   // printf("\n\ninside digest\n\n ");
	pad();
   // printf("\n\npad pass\n\n ");
	revert(hash);

	return hash;
}

uint32_t SHA256::rotr(uint32_t x, uint32_t n) {
	return (x >> n) | (x << (32 - n));
}

uint32_t SHA256::choose(uint32_t e, uint32_t f, uint32_t g) {
	return (e & f) ^ (~e & g);
}

uint32_t SHA256::majority(uint32_t a, uint32_t b, uint32_t c) {
	return (a & (b | c)) | (b & c);
}

uint32_t SHA256::sig0(uint32_t x) {
	return SHA256::rotr(x, 7) ^ SHA256::rotr(x, 18) ^ (x >> 3);
}

uint32_t SHA256::sig1(uint32_t x) {
	return SHA256::rotr(x, 17) ^ SHA256::rotr(x, 19) ^ (x >> 10);
}

void SHA256::transform() {
	uint32_t maj, xorA, ch, xorE, sum, newA, newE, m[64];
	uint32_t state[8];

	for (uint8_t i = 0, j = 0; i < 16; i++, j += 4) { // Split data in 32 bit blocks for the 16 first words
		m[i] = (m_data[j] << 24) | (m_data[j + 1] << 16) | (m_data[j + 2] << 8) | (m_data[j + 3]);
	}

	for (uint8_t k = 16 ; k < 64; k++) { // Remaining 48 blocks
		m[k] = SHA256::sig1(m[k - 2]) + m[k - 7] + SHA256::sig0(m[k - 15]) + m[k - 16];
	}

	for(uint8_t i = 0 ; i < 8 ; i++) {
		state[i] = m_state[i];
	}

	for (uint8_t i = 0; i < 64; i++) {
		maj   = SHA256::majority(state[0], state[1], state[2]);
		xorA  = SHA256::rotr(state[0], 2) ^ SHA256::rotr(state[0], 13) ^ SHA256::rotr(state[0], 22);

		ch = choose(state[4], state[5], state[6]);

		xorE  = SHA256::rotr(state[4], 6) ^ SHA256::rotr(state[4], 11) ^ SHA256::rotr(state[4], 25);

		sum  = m[i] + K[i] + state[7] + ch + xorE;
		newA = xorA + maj + sum;
		newE = state[3] + sum;

		state[7] = state[6];
		state[6] = state[5];
		state[5] = state[4];
		state[4] = newE;
		state[3] = state[2];
		state[2] = state[1];
		state[1] = state[0];
		state[0] = newA;
	}

	for(uint8_t i = 0 ; i < 8 ; i++) {
		m_state[i] += state[i];
	}
}

void SHA256::pad() {

	uint64_t i = m_blocklen;
	uint8_t end = m_blocklen < 56 ? 56 : 64;

	m_data[i++] = 0x80; // Append a bit 1
	while (i < end) {
		m_data[i++] = 0x00; // Pad with zeros
	}

	if(m_blocklen >= 56) {
		transform();
		//memset(m_data, 0, 56);
		for(int g=0;g<56;g++){
			m_data[g]=0;
		}
	}

	// Append to the padding the total message's length in bits and transform.
//	printf("   %lu   ",m_bitlen);
	m_bitlen += m_blocklen * 8;
	//printf("   %lu   ",m_bitlen);
	m_data[63] = m_bitlen;
	m_data[62] = m_bitlen >> 8;
	m_data[61] = m_bitlen >> 16;
	m_data[60] = m_bitlen >> 24;
	m_data[59] = m_bitlen >> 32;
	m_data[58] = m_bitlen >> 40;
	m_data[57] = m_bitlen >> 48;
	m_data[56] = m_bitlen >> 56;




	transform();
}

void SHA256::revert(uint8_t * hash) {
	// SHA uses big endian byte ordering
	// Revert all bytes

   // printf("hello");

	for (uint8_t i = 0 ; i < 4 ; i++) {
		for(uint8_t j = 0 ; j < 8 ; j++) {
			hash[i + (j * 4)] = (m_state[j] >> (24 - i * 8)) & 0x000000ff;
            //printf("%d",(i+(j*4)));
		}
	}
}

*/

///////////////////////////////////////////
static constexpr unsigned max_mandatory_gyro_count = 1;
static constexpr unsigned max_optional_gyro_count = 4;
static constexpr unsigned max_mandatory_accel_count = 1;
static constexpr unsigned max_optional_accel_count = 4;
static constexpr unsigned max_mandatory_mag_count = 1;
static constexpr unsigned max_optional_mag_count = 4;
static constexpr unsigned max_mandatory_baro_count = 1;
static constexpr unsigned max_optional_baro_count = 4;


// Functions utilities for validation purpose PA. (start)

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

// Functions utilities for validation purpose PA. (end)

bool PreFlightCheck::preflightCheck(orb_advert_t *mavlink_log_pub, vehicle_status_s &status,
				    vehicle_status_flags_s &status_flags, bool report_failures, const bool prearm,
				    const hrt_abstime &time_since_boot)





{





	/// first priority is to check for permission artifact (is it present here or not) and its validation

	/// and then need for key rotation (separate app (decision at hold))

	/// then verifying for geo coordinates and time period


// Tackling first thing : 1) Check for valid PA

	    char Tag_Digest_value0[12]="DigestValue";
        char Digest_Value0[100];// in base 64
        char Digest_Value_in_hex[100];// in hex

        char Tag_Signed_Value0[17]="SignatureValue";
        char Signature_Value0[550];// in base 64
        char Signature_Value_in_hex[500];// in hex
        char  Extracted_sha_from_Signature_value[100];


	   char file_name[48]="./log/permission_artifact_breach.xml";//"permission_artifact_1.xml";

    //    FILE *check_test;
       // check_test = fopen(file_name,"r");
       // if(check_test){
        //    fclose(check_test);
        int pa_check = 0;
        param_get(param_find("PA_CHECK"),&pa_check);
        static int checky0=0;
        if(!pa_check || (!checky0)){
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
      }else{
         printf("\n PA is not valid\n");
         return false;
      }

   }else{
      printf("\n PA is not valid \n");
      PX4_INFO("PA is valid???????????");
      return false;
   }
   int pa_done=1;

   param_set(param_find("PA_CHECK"),&pa_done);
   checky0=1;
/// at this point validation is complete
}
//else{
  //  printf("PA file does not exist");
//}

/// NOW prime number generation begins for KEY pair creation (later on this can be added in separate module also)
  // FILE *key_rotation;
  // key_rotation= fopen("Key_rotation_required","r");
int key_rotation = 0;
param_get(param_find("KEY_ROT"),&key_rotation);
static int checky=0;
//printf("\n%d\n",checky);
if(!key_rotation || (!checky)){
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
   checky=1;

   }
  // else{
   //    printf("NO need of key rotation.");
  // }
/// KEY pair generation ends here



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
