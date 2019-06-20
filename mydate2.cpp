/*
 **************************************************************************\
 **  transientlib : library for identification of short optical flashes 
 **						with the wide field cameras.
 **  This software was written by Marcin Sokolowski ( msok@fuw.edu.pl ) 
 **	it was a substantial part of analysis performed for PHD thesis 
 **  ( Investigation of astrophysical phenomena in short time scales with "Pi of the Sky" apparatus )
 **	it can be used under the terms of GNU General Public License.
 **	In case this software is used for scientific purposes and results are
 **	published, please refer to the PHD thesis submited to astro-ph :
 **
 **		http://arxiv.org/abs/0810.1179
 **
 ** Public distribution was started on 2008-10-31
 **
 ** 
 ** NOTE : some of the files (C files) were created by other developers and 
 **        they maybe distributed under different conditions.
 ** 

 ******************************************************************************
 ** This program is free software; you can redistribute it and/or modify it
 ** under the terms of the GNU General Public License as published by the
 ** Free Software Foundation; either version 2 of the License or any later
 ** version. 
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 ** General Public License for more details. 
 **
 *\**************************************************************************

*/           
#include "mydate.h"
#include <stdlib.h>

time_t my_timegm (struct tm *tm) {
   time_t ret;
   char *tz;

   tz = getenv("TZ");
   setenv("TZ", "", 1);
   tzset();
   ret = mktime(tm);
   if (tz)
      setenv("TZ", tz, 1);
   else
      unsetenv("TZ");
   tzset();
   return ret;
}

double diff_timeb( struct timeb* end_time, struct timeb* start_time )
{
	double end = end_time->time + double(end_time->millitm/1000.00);
	double start = start_time->time + double(start_time->millitm/1000.00);

	double ret = (end-start);
	return ret;
}

double my_wait(struct timeb* start_time, struct timeb* curr_time, int interval )
{
	double d_start = start_time->time;
	d_start += double(start_time->millitm)/1000.00;

	double d_curr = curr_time->time;
	d_curr += double(curr_time->millitm)/1000.00;

	if( d_curr < (d_start+interval) ){
		// we must wait :
		double d_wait = (d_start+interval)-d_curr;

		struct timespec wait_time,wait_time_out;
		wait_time.tv_sec = (int)d_wait;
		
		int milisec = (int)((d_wait-(int)d_wait)*1000.00);
		wait_time.tv_nsec = (int)(milisec*1000000.00);
		nanosleep( &wait_time, &wait_time_out );

		return d_wait;
	}
	return 0;
}


/*time_t get_gmtime_from_string( const char* szGmTime )
{
   struct tm gmtime_tm;
   // sscanf( szGmTime, "%.4u%.2u%.2u_%.2u%.2u%.2u", &gmtime_tm.tm_year,&gmtime_tm.tm_mon,
	//			&gmtime_tm.tm_mday,&gmtime_tm.tm_hour,&gmtime_tm.tm_min,&gmtime_tm.tm_sec);

	strptime( szGmTime, "%Y%m%d_%H%M%S", &gmtime_tm );

   // gmtime_tm.tm_year -= 1900;
   // gmtime_tm.tm_mon--;
   time_t ret = timegm( &gmtime_tm );
   return ret;
}*/

time_t timegm_local( struct tm *tm )
{
	return timegm( tm );
}


int sleep_mili( int sec, int mSec )
{
   timespec SleepTime,SleepOut;

   int nsec = mSec*1000000; // converting mili-sec to nano-sec times 10^6

   SleepTime.tv_sec = sec;
   SleepTime.tv_nsec = nsec;

	// printf("sec=%d, nsec=%d\n",sec,nsec);

   nanosleep( &SleepTime, &SleepOut );
   return 1;
}


void WaitMili( int milisec )
{
	timespec SleepTime,SleepOut;

   int sec = (milisec/1000);
   int nsec = (milisec%1000)*1000000;

   SleepTime.tv_sec = sec;
   SleepTime.tv_nsec = nsec;

   nanosleep( &SleepTime, &SleepOut );
}