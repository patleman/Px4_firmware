#include"print_hello.hpp"

void printing_hello_world(){
printf("\n\nhello world\n\n");
}

void s_mp_zero_digs(mp_digit *d, int digits)
{
   while (digits-- > 0) {
      *d++ = 0;
   }
}
