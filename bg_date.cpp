#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/timeb.h>
#include <string.h>

#include <string>

#include <mydate.h>

using namespace std;

/*time_t get_dttm()
{  
  long gm_time;
  time( &gm_time );
  return gm_time;
}*/
      

time_t get_unixtime_from_local_string2( const char* szDTM, const char* format )
{
   struct tm local_time_tm;
   memset( &local_time_tm, '\0', sizeof(struct tm));
       
   // temporary correction due to fact that strptime does not fill fields :
   // tm_isdst = 1, tm_gmtoff = -10800,  tm_zone = 0x85e85a8 "CLST"
   // thus not working exactly good ... , but this is now 
   // filling current values of this field , which may not work for past 
   // and future dates ...
   time_t ut_time=get_dttm();
   localtime_r( &ut_time , &local_time_tm );
                            
   // 11:29:59 11.07.2012
   strptime( szDTM, format, &local_time_tm );
   time_t ret = mktime( &local_time_tm );
   return ret;    
}


/*string get_gmtime_string( time_t ut_time )
{
   struct tm gmtime_tm;
   string szRet;
   if(gmtime_r( &ut_time, &gmtime_tm )){
      char tempstring[64];

      // bug ??? first %.2u -> %.4u ???
      sprintf(tempstring,"%.2u%.2u%.2u_%.2u%.2u%.2u",
                           gmtime_tm.tm_year+1900,(gmtime_tm.tm_mon+1),gmtime_tm.tm_mday,
                           gmtime_tm.tm_hour,gmtime_tm.tm_min,gmtime_tm.tm_sec);
      szRet=tempstring;
   }  
   return szRet;
}*/
string get_gmtime_string_bg( time_t ut_time )
{
   mystring tmp = get_gmtime_string(ut_time);
   string ret = tmp.c_str();
   
   return ret;
}

string get_gmtime_string_bg()
{
   return get_gmtime_string_bg( get_dttm() );
}

/*string get_localtime_string( time_t ut_time )
{
   char szDate[40];
   struct tm _tm;
   string szRet;   
       
   localtime_r(&ut_time , &_tm);
             
   int year = _tm.tm_year+1900;
   sprintf(szDate,"%d%02d%02d_%02d%02d%02d",year,_tm.tm_mon+1,_tm.tm_mday,
                   _tm.tm_hour,_tm.tm_min,_tm.tm_sec);

   szRet = szDate;                       
   return szRet;  
}*/
string get_localtime_string( time_t ut_time )
{
   mystring tmp = get_date_time_string( ut_time );
   string ret = tmp.c_str();
   
   return ret;
}

/*int sleep_mili( int sec, int mSec )
{
   timespec SleepTime,SleepOut;
                     
   int nsec = mSec*1000000; // converting mili-sec to nano-sec times 10^6
                                             
   SleepTime.tv_sec = sec;
   SleepTime.tv_nsec = nsec;
                                                                    
   // printf("sec=%d, nsec=%d\n",sec,nsec);
                                                                              
   nanosleep( &SleepTime, &SleepOut );
   return 1;
}*/

                                                                                        
                                                                                                                                

time_t doy2dttm_local( int year, int doy, int hour, int min, int sec, int& local_date )
{
   struct tm _tm;
   memset( &_tm, '\0', sizeof(struct tm));
   if( year>1900 ){
      _tm.tm_year = year-1900;
   }


   _tm.tm_hour = 1;
   _tm.tm_min = 0;
   _tm.tm_sec = 0;
   _tm.tm_mday = 1;
   _tm.tm_mon = 0;
   time_t mid_new_year = timegm( &_tm );
   int day = 24*3600;

   time_t dttm = mid_new_year + day*(doy-1);


   struct tm _tm2;
   gmtime_r( &dttm  , &_tm2 );

   _tm2.tm_hour = hour;
   _tm2.tm_min = min;
   _tm2.tm_sec = sec;
   time_t ut_time = timegm( &_tm2 );
   
   struct tm gmtm;
   time_t _gm = ut_time-12*3600; // not to have filename change after midnight :
   localtime_r( &_gm , &gmtm );
   int year_full=1900+gmtm.tm_year; // in gmtm.tm_year years since 1900
             
   char szTmp[128];
   sprintf( szTmp, "%.4u%.2u%.2u" , year_full,(gmtm.tm_mon+1),gmtm.tm_mday );                    
   local_date = atol( szTmp );
                             
   return ut_time;
}


double get_local_hour_decimal( time_t ut_time )
{
   struct tm gmtm;
   time_t _gm = ut_time;
   localtime_r( &_gm , &gmtm );

   double sec = gmtm.tm_min*60 + gmtm.tm_sec;
   double ret = gmtm.tm_hour + sec/3600.00;
   return ret;
}

double get_ut_hour_decimal( time_t ut_time )
{
   struct tm gmtm;
   time_t _gm = ut_time;
   gmtime_r( &_gm , &gmtm );

   double sec = gmtm.tm_min*60 + gmtm.tm_sec;
   double ret = gmtm.tm_hour + sec/3600.00;
   return ret;
}

                     