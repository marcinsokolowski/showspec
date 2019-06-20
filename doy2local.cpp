#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "bg_date.h"

void usage()
{
   printf("ux2local YEAR[YYYY] DOY\n");
   exit(-1);
}

int main(int argc,char* argv[])
{
   // time_t doy2dttm_local( int year, int doy, int hour, int min, int sec, int& local_date );   
   int year = atol(argv[1]);
   int doy = atol(argv[2]);
   
   int local_date;
   time_t uxtime = doy2dttm_local(year,doy,12,0,0,local_date);
   printf("%d (%d)\n",local_date,uxtime);
}
