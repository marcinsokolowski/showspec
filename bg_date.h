#ifndef _BG_DATE_H__
#define _BG_DATE_H__

#include <string>
using namespace std;

// date/time :
time_t get_dttm(); // in baselib/mydate.h 
string get_gmtime_string_bg();
string get_gmtime_string_bg( time_t ut_time );
string get_localtime_string( time_t ut_time );
int sleep_mili( int sec, int mSec );

time_t get_unixtime_from_local_string2( const char* szDTM, const char* format="%H:%M:%S %m.%d.%Y" );

// included externlly from baselib :
// void get_ymd_hms( time_t ut_time, int& year, int& month, int& day, int& hour, int& minute, int& sec );

time_t doy2dttm_local( int year, int doy, int hour, int min, int sec, int& local_date );                  

double get_local_hour_decimal( time_t ut_time );
double get_ut_hour_decimal( time_t ut_time );


#endif
