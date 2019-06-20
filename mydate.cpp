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
// MyDate.cpp: implementation of the CMyDate class.
//
//////////////////////////////////////////////////////////////////////
#define _XOPEN_SOURCE /* glibc2 needs this */

#include <time.h>
#include <sys/time.h>
#include <sys/timeb.h>

#include "mydate.h"
#include "cexcp.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

// global functions :

CTimer gTimer;

void show_time( BOOL_T bNewLine )
{
	mystring szDTM = get_date_time_string();
	printf("[ %s ]",szDTM.c_str());
	if( bNewLine ){
		printf("\n");
	}	
	fflush(0);
}

void get_night_date_gmt( mystring& szNIGHT_DATE )
{
	get_night_date_gmt( get_dttm(), szNIGHT_DATE );
}

void get_night_date_gmt( time_t ut_time , mystring& szNIGHT_DATE )
{
	struct tm gmtm;
   time_t _gm = ut_time-12*3600; // not to have filename change after midnight :
   gmtime_r( &_gm , &gmtm );
   int year=1900+gmtm.tm_year; // in gmtm.tm_year years since 1900
	
	char szTmp[128];
	sprintf( szTmp, "%.4u%.2u%.2u" , year,(gmtm.tm_mon+1),gmtm.tm_mday );
	
	szNIGHT_DATE = szTmp;	
}

mystring get_night_date_local()
{
	mystring szNIGHT_DATE;
	get_night_date_local( get_dttm(), szNIGHT_DATE );
	return szNIGHT_DATE;
}

void get_night_date_local( mystring& szNIGHT_DATE )
{
	get_night_date_local( get_dttm(), szNIGHT_DATE );
}

void get_ymd_hms( time_t ut_time, int& year, int& month, int& day,
                  int& hour, int& minute, int& sec )
{
	struct tm gmtm;
	localtime_r( &ut_time , &gmtm );
	year = 1900 + gmtm.tm_year; // + 1900 
   month = (gmtm.tm_mon+1);
   day = gmtm.tm_mday;
	hour = gmtm.tm_hour;
	minute = gmtm.tm_min;
	sec = gmtm.tm_sec;
}

void get_ymd_hms_ut( time_t ut_time, int& year, int& month, int& day,
                  int& hour, int& minute, int& sec )
{
	struct tm gmtm;
	gmtime_r( &ut_time , &gmtm );
	year = 1900 + gmtm.tm_year; // + 1900
   month = (gmtm.tm_mon+1);
   day = gmtm.tm_mday;
	hour = gmtm.tm_hour;
	minute = gmtm.tm_min;
	sec = gmtm.tm_sec;
}


void get_night_ymd( time_t ut_time, int& year, int& month, int& day )
{
   struct tm gmtm;
   time_t _gm = ut_time-12*3600; // not to have filename change after midnight :
   localtime_r( &_gm , &gmtm );

   char szTmp[128];
   year = 1900 + gmtm.tm_year; // + 1900
   month = (gmtm.tm_mon+1);
   day = gmtm.tm_mday;
}


void get_night_date_local( time_t ut_time , mystring& szNIGHT_DATE )
{
	struct tm gmtm;
   time_t _gm = ut_time-12*3600; // not to have filename change after midnight :
   localtime_r( &_gm , &gmtm );
	
	int year=1900+gmtm.tm_year; // in gmtm.tm_year years since 1900 
	
	char szTmp[128];
	sprintf( szTmp, "%.4u%.2u%.2u" ,year,(gmtm.tm_mon+1),gmtm.tm_mday );
	
	szNIGHT_DATE = szTmp;	
}

int get_night_date_local( time_t ut_time )
{
	struct tm gmtm;
   time_t _gm = ut_time-12*3600; // not to have filename change after midnight :
   localtime_r( &_gm , &gmtm );
	int year=1900+gmtm.tm_year; // in gmtm.tm_year years since 1900
	
	char szTmp[128];
	sprintf( szTmp, "%.4u%.2u%.2u" , year,(gmtm.tm_mon+1),gmtm.tm_mday );
	
	return atol( szTmp );	
}



time_t get_uttime_rounded_to_5( time_t ut_time )
{
	struct tm gmtm;
   time_t _gm = ut_time;
   localtime_r( &_gm , &gmtm );
	
	char szTmp[128];
	sprintf( szTmp, "%.2u%.2u",gmtm.tm_hour,gmtm.tm_min);

	int hm = atol(szTmp);	
	int hm_r = round_hm( hm );

	printf("hm=%d, hm_r=%d, gm=%d:%d\n",hm,hm_r,gmtm.tm_hour,gmtm.tm_min);

	gmtm.tm_hour = (hm_r/100);
	gmtm.tm_min  = (hm_r%100);

	printf("after gm=%d:%d\n",gmtm.tm_hour,gmtm.tm_min);

	time_t ret = mktime( &gmtm );
	return ret;
}

int get_local_hour( time_t ut_time )
{
	struct tm gmtm;
   time_t _gm = ut_time;
   localtime_r( &_gm , &gmtm );
	
	char szTmp[128];
	sprintf( szTmp, "%.2u%.2u",gmtm.tm_hour,gmtm.tm_min);

	return atol(szTmp);
}

void get_short_night_date_local( time_t ut_time , mystring& szNIGHT_DATE )
{
	struct tm gmtm;
   time_t _gm = ut_time-12*3600; // not to have filename change after midnight :
   localtime_r( &_gm , &gmtm );
	
	char szTmp[128];
	sprintf( szTmp, "%.2u%.2u%.2u" , (gmtm.tm_year-100),(gmtm.tm_mon+1),gmtm.tm_mday );
	
	szNIGHT_DATE = szTmp;	
}


void convert_date_from_fits_header_format( struct tm* _tm )
{
	// in header years are in format YYYY :
	_tm->tm_year = (_tm->tm_year - 1900);

	// and months from 1 :
	_tm->tm_mon--;
}

int get_time_zone()
{
	time_t dttm;// = get_dttm();
	tm TimeInfo;
	localtime_r( &dttm, &TimeInfo );
	
	tzset();
	return timezone;
}


mystring get_gmtime_string()
{
	time_t ut = get_dttm();
	mystring szRet = get_gmtime_string( ut );
	return szRet;
}

mystring get_gmtime_string( time_t ut_time )
{
	struct tm gmtime_tm;
	mystring szRet;
	if(gmtime_r( &ut_time, &gmtime_tm )){
		char tempstring[64];

		// bug ??? first %.2u -> %.4u ???
		sprintf(tempstring,"%.2u%.2u%.2u_%.2u%.2u%.2u",
									gmtime_tm.tm_year+1900,(gmtime_tm.tm_mon+1),gmtime_tm.tm_mday,
									gmtime_tm.tm_hour,gmtime_tm.tm_min,gmtime_tm.tm_sec);
		szRet << tempstring;
	}	
	return szRet;
}

time_t get_unixtime_from_local_string( const char* szDTM )
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

   strptime( szDTM, "%Y%m%d_%H%M%S", &local_time_tm );
	time_t ret = mktime( &local_time_tm );
	return ret;		
}

time_t get_unixtime_from_local_string2( const char* szDTM )
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
   strptime( szDTM, "%H:%M:%S %m.%d.%Y", &local_time_tm );
	time_t ret = mktime( &local_time_tm );
	return ret;		
}


time_t get_gmtime_from_string( const char* szGmTime )
{
   struct tm gmtime_tm;
   // sscanf( szGmTime, "%.4u%.2u%.2u_%.2u%.2u%.2u", &gmtime_tm.tm_year,&gmtime_tm.tm_mon,
   // &gmtime_tm.tm_mday,&gmtime_tm.tm_hour,&gmtime_tm.tm_min,&gmtime_tm.tm_sec);

   strptime( szGmTime, "%Y%m%d_%H%M%S", &gmtime_tm );

   // gmtime_tm.tm_year -= 1900;
   // gmtime_tm.tm_mon--;
   time_t ret = timegm_local( &gmtime_tm );
	// time_t ret = timegm( &gmtime_tm );
   return ret;
}


BASELIB_EI time_t local_to_dttm(struct tm* _localtime)
{
	return mktime(_localtime);	
}

BASELIB_EI tm dttm_to_local(time_t _dttm)
{
	struct tm ret = *localtime( &_dttm );
	return ret;
}

BASELIB_EI time_t get_dttm()
{
	long gm_time;
	time( &gm_time );
	return gm_time;
}



int get_sod()
{
	time_t t = get_dttm();
	time_t mid = get_gmmidnight_ut( t );
	return (t-mid);
}

int get_sod( time_t t )
{
	time_t mid = get_gmmidnight_ut( t );
	return (t-mid);
}



mystring get_date()
{
	mystring ret;
	ret << (long)get_mydate();
	return ret;
}


MYDATE_T get_mydate( time_t ut_time )
{
	return dttm_to_mydate( ut_time );
}

BASELIB_EI MYDATE_T get_mydate()
{
	time_t tt = get_dttm();	
	return dttm_to_mydate( tt );
}

BASELIB_EI MYDATE_T dttm_to_mydate(time_t _dttm)
{
	struct tm* _tm;
	_tm = gmtime(&_dttm);
	return ((_tm->tm_year+1900)*10000+(_tm->tm_mon+1)*100+_tm->tm_mday);
}

int dttm2grbdate( time_t _dttm, char* grbdate ){
	struct tm* _tm;
   _tm = gmtime(&_dttm);

	sprintf(grbdate,"%.2u%.2u%.2u",(_tm->tm_year-100),(_tm->tm_mon+1),_tm->tm_mday);
	return atol( grbdate );
}

BASELIB_EI mystring dttm2dbstring(time_t _dttm)
{
	mystring szRet;

	if (_dttm!=0){
		char szDate[40];
		struct tm* _tm;
		_tm = localtime( &_dttm );
		int year = _tm->tm_year+1900;
	
		sprintf(szDate,"%d%.2d%.2d%.2d%.2d%.2d",year,_tm->tm_mon+1,_tm->tm_mday,
		                               _tm->tm_hour,_tm->tm_min,_tm->tm_sec);
		Assert(strlen(szDate)<sizeof(szDate),"Date buffer size exceeded");
		szRet << szDate;
	}else{
		szRet << "NULL";
	}
	return szRet;
}

mystring get_date_time_string( time_t ut_time )
{
	char szDate[40];
	struct tm _tm;


	localtime_r(&ut_time , &_tm);
	
	int year = _tm.tm_year+1900;
	sprintf(szDate,"%d%02d%02d_%02d%02d%02d",year,_tm.tm_mon+1,_tm.tm_mday,
		                               _tm.tm_hour,_tm.tm_min,_tm.tm_sec);
	Assert(strlen(szDate)<sizeof(szDate),"Date buffer size exceeded");
	return mystring(szDate);	

}

BASELIB_EI mystring get_date_time_string()
{
	char szDate[40];
	time_t _dttm;
	struct tm _tm;

	time(&_dttm);
	// _tm = gmtime(&_dttm);

	localtime_r( &_dttm , &_tm );
	
	int year = _tm.tm_year+1900;
	sprintf(szDate,"%d%02d%02d_%02d%02d%02d",year,_tm.tm_mon+1,_tm.tm_mday,
		                               _tm.tm_hour,_tm.tm_min,_tm.tm_sec);
	Assert(strlen(szDate)<sizeof(szDate),"Date buffer size exceeded");
	return mystring(szDate);	
}

BASELIB_EI const char* get_dttm_db_string()
{
	static char szDate[40];
	time_t _dttm;
	struct tm* _tm;

	time(&_dttm);
	_tm = gmtime(&_dttm);
	
	int year = _tm->tm_year+1900;
	sprintf(szDate,"%d%02d%02d%02d%02d%02d",year,_tm->tm_mon+1,_tm->tm_mday,
		                               _tm->tm_hour,_tm->tm_min,_tm->tm_sec);
	Assert(strlen(szDate)<sizeof(szDate),"Date buffer size exceeded");
	return szDate;
}


// returns unix_time of LOCAL MIDNIGHT !!!
time_t get_midnight_ut( time_t ut_time )
{
	struct tm gmtm;
   time_t _gm = ut_time;

	localtime_r( &_gm , &gmtm );
	gmtm.tm_hour = 0;
	gmtm.tm_min = 0;
	gmtm.tm_sec =0;
	time_t midt = mktime( &gmtm );		

	return midt;
}

time_t get_midday_ut( time_t ut_time )
{
	struct tm gmtm;
   time_t _gm = ut_time;

	localtime_r( &_gm , &gmtm );
	gmtm.tm_hour = 12;
	gmtm.tm_min = 0;
	gmtm.tm_sec =0;
	time_t midt = mktime( &gmtm );		

	return midt;
}

time_t get_gmmidnight_ut( time_t ut_time )
{
	struct tm gmtm;
   time_t _gm = ut_time;

	gmtime_r( &_gm , &gmtm );
	gmtm.tm_hour = 0;
	gmtm.tm_min = 0;
	gmtm.tm_sec =0;
	time_t midt = timegm( &gmtm );		

	return midt;
}


BASELIB_EI int sec_to_local_midnight(time_t _dttm)
{
	struct tm* _tm;
	_tm = localtime( &_dttm );
	int sec = (24-_tm->tm_hour)*60*60+(60-_tm->tm_min)*60+(60-_tm->tm_sec);
	return sec;
}

BASELIB_EI int sec_to_gm_midnight(time_t _dttm)
{
	struct tm* _tm;
	_tm = gmtime( &_dttm );
	int sec = (24-_tm->tm_hour)*60*60+(60-_tm->tm_min)*60+(60-_tm->tm_sec);
	return sec;
}

double msec2sec( int msec )
{
	double sec = ((double)(msec))/1000.00;
	return sec;
}

CMyDate::CMyDate()
{

}

CMyDate::~CMyDate()
{

}

const char* CMyDate::getdate(const char* fmt)
{
	static mystring szDate;
	struct tm* _tm;
	time_t _dttm;


	time(&_dttm);
	_tm = localtime( &_dttm );
	int year = _tm->tm_year+1900;
	szDate = "";
	szDate << _tm->tm_mon+1<< "-" << _tm->tm_mday << "-" 
			 << year << " " << _tm->tm_hour << ":";
	if (_tm->tm_min<10)
		szDate << "0";			 
	szDate << _tm->tm_min << ":";
	if (_tm->tm_sec<10)
		szDate << "0";
	szDate << _tm->tm_sec;
	return szDate.c_str();
}

mystring get_clock_in_sec_string( long cl)
{
	double r = double(cl)/double(CLOCKS_PER_SEC);
	long msec = (long)(r*100.00);
	long sec = (msec/100);
	long l_msec = (msec%100);
	mystring szStr;
	szStr << cl << " ticks (=" << sec << ".";
	if(l_msec<10)
		szStr << "0";
	szStr << l_msec << " sec)";
	return szStr;
}


CTimer::CTimer()
{
	start_t = clock();
}

void CTimer::StartTimer()
{
	start_t = clock();
}

mystring CTimer::GetTimer()
{
	clock_t t2=clock();
   mystring msg=get_clock_in_sec_string( t2-start_t );
	return msg;
}



time_t CMyDate::getTime( const char* szDT, const char* fmt )
{
	struct tm _tm;
	// memset(&_tm,'\0',sizeof(_tm));
	time_t _t = get_dttm();
	_tm = *::localtime( &_t );
	strptime( szDT	, fmt, &_tm );
	time_t ret = mktime( &_tm );	
	return ret;
}

void CMyDate::SetCLT()
{
   setenv("TZ","/usr/share/zoneinfo/Chile/Continental",1);
   tzset();
}


time_t get_runupto_dtm( const char* szDTM, const char* fmt )
{
	time_t upto=0;
	if(strstr(szDTM,"_")){
		upto = CMyDate::getTime( szDTM, fmt );
	}else{
		// only time provided :	
		mystring szCurrDate = get_date();
		mystring szDateTime;
		szDateTime << szCurrDate << "_" << szDTM;
		time_t currTime = get_dttm();
		time_t tt =  CMyDate::getTime( szDateTime.c_str(), fmt );
		if(tt<currTime){
			// next day :
			tt += (3600*24);
		}	
		upto = tt;					
	}
	return upto;
}


int round_hm( int hm )
{
	int ret = hm;
	int m = (hm % 100);
	int h = (hm / 100);
	if( m % 5 > 0 ){
		int new_m = (m - (m % 5))+5;
		if( new_m >= 60 ){
			h++;
			new_m -= 60;			
		}
		ret = h*100+new_m;
	}
	return ret;
}


time_t doy2dttm( int year, int doy, int hour, int min, int sec )
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

	for( int i=0;i<=367;i++){
		// DOY is counted from 1-365 :
		time_t dttm = mid_new_year + day*i;


		struct tm _tm2;
		gmtime_r( &dttm  , &_tm2 );
		
		int as_doy = _tm2.tm_yday+1;

		if( as_doy == doy ){
			_tm2.tm_hour = hour;
			_tm2.tm_min = min;
			_tm2.tm_sec = sec;
			time_t ret = timegm( &_tm2 );
			return ret;
		}
	}

	return 0;
}

int get_doy()
{
	time_t dttm=get_dttm();
	struct tm res_tm;
	localtime_r( &dttm, &res_tm);
	return (res_tm.tm_yday+1);
}

double get_doyd( time_t dttm ){
	struct tm res_tm;
	gmtime_r( &dttm, &res_tm);

	double sec = res_tm.tm_hour*3600.00+res_tm.tm_min*60.00+res_tm.tm_sec;	
	double doyd = res_tm.tm_yday + sec/(24.00*3600.00);

	return doyd;
	
}

double get_epoch( time_t dttm )
{
	struct tm res_tm;
	gmtime_r( &dttm, &res_tm);

	double sec = res_tm.tm_hour*3600.00+res_tm.tm_min*60.00+res_tm.tm_sec;	
	double doyd = res_tm.tm_yday + sec/(24.00*3600.00);

	int year = 1900 + res_tm.tm_year;
	int year_mod = (year % 1000);

	double epoch = year_mod*1000+doyd;
	return epoch;	
}

double get_epoch_new( time_t dttm )
{
	struct tm res_tm;
	gmtime_r( &dttm, &res_tm);

	double sec = res_tm.tm_hour*3600.00+res_tm.tm_min*60.00+res_tm.tm_sec;	
	double doyd = res_tm.tm_yday + sec/(24.00*3600.00);

	int year = 1900 + res_tm.tm_year;

	double epoch = year + (doyd/365.00);
	return epoch;	
}

