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
// MyDate.h: interface for the CMyDate class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MYDATE_H__C881F1AA_AE7A_11D5_B636_382E07C10000__INCLUDED_)
#define AFX_MYDATE_H__C881F1AA_AE7A_11D5_B636_382E07C10000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <time.h>
#include <sys/timeb.h>
#include <stdio.h>
#include "mytypes.h"
#include "mystring.h"
#include "basedefines.h"

#define DAY_IN_SEC (60*60*24)

// DOY - day of the year calculations :
int get_doy();
double get_doyd( time_t dttm );
double get_epoch( time_t dttm );
double get_epoch_new( time_t dttm );

int round_hm( int hm );

void show_time( BOOL_T bNewLine=TRUE );

void convert_date_from_fits_header_format( struct tm* _tm );

time_t doy2dttm( int year, int doy, int hour, int min, int sec );

time_t my_timegm (struct tm *tm);

time_t timegm_local( struct tm *tm );

void WaitMili( int milisec );

double my_wait(struct timeb* start_time, struct timeb* curr_time, int interval );

double diff_timeb( struct timeb* end_time, struct timeb* start_time );

BASELIB_EI mystring get_date_time_string();

BASELIB_EI mystring get_date_time_string( time_t ut_time );

BASELIB_EI time_t local_to_dttm(struct tm* _localtime);

BASELIB_EI tm dttm_to_local(time_t _dttm);

BASELIB_EI time_t get_dttm();

int get_sod();

int get_sod( time_t t );

BASELIB_EI MYDATE_T get_mydate();

int dttm2grbdate( time_t _dttm, char* grbdate );

MYDATE_T get_mydate( time_t ut_time );

mystring get_date();

BASELIB_EI MYDATE_T dttm_to_mydate(time_t _dttm);

BASELIB_EI const char* get_dttm_db_string();

BASELIB_EI mystring dttm2dbstring(time_t _dttm);

BASELIB_EI int sec_to_local_midnight(time_t _dttm);

BASELIB_EI int sec_to_gm_midnight(time_t _dttm);

// returns unix_time of LOCAL MIDNIGHT !!!
time_t get_midnight_ut( time_t ut_time );
time_t get_midday_ut( time_t ut_time );

// returns unix_time of GM MIDNIGHT !!!!
time_t get_gmmidnight_ut( time_t ut_time );


BASELIB_EI mystring get_clock_in_sec_string( long cl);

mystring get_gmtime_string( time_t ut_time );

mystring get_gmtime_string();

time_t get_gmtime_from_string( const char* szGmTime );
time_t get_unixtime_from_local_string( const char* szDTM );
time_t get_unixtime_from_local_string2( const char* szDTM );

int get_time_zone();


void get_night_date_gmt( mystring& szNIGHT_DATE );
void get_night_date_local( mystring& szNIGHT_DATE );
mystring get_night_date_local();

time_t get_uttime_rounded_to_5( time_t ut_time );
int get_local_hour( time_t ut_time );
void get_night_ymd( time_t ut_time, int& year, int& month, int& day );
void get_short_night_date_local( time_t ut_time , mystring& szNIGHT_DATE );
void get_night_date_local( time_t ut_time , mystring& szNIGHT_DATE );
int get_night_date_local( time_t ut_time );
void get_night_date_gmt( time_t ut_time , mystring& szNIGHT_DATE );

void get_ymd_hms( time_t ut_time, int& year, int& month, int& day,
						int& hour, int& minute, int& sec );

void get_ymd_hms_ut( time_t ut_time, int& year, int& month, int& day,
						int& hour, int& minute, int& sec );

time_t get_runupto_dtm( const char* szDTM, const char* fmt="%Y%m%d_%H%M" );

int sleep_mili( int sec, int mSec );

double msec2sec( int msec );

class BASELIB_EI CMyDate  
{
public:
	CMyDate();
	virtual ~CMyDate();
	static const char* getdate(const char* fmt="MM/DD/YYYY HH:MI:SS");
	static time_t getTime( const char* szDT, const char* fmt );
	static void SetCLT();
protected:
	mystring m_szDate;
	time_t m_gmtime;
};

class CTimer
{
public :
	CTimer();
	void StartTimer();
	void StopTimer();
	mystring GetTimer();

	clock_t start_t;
};

extern CTimer gTimer;

#endif // !defined(AFX_MYDATE_H__C881F1AA_AE7A_11D5_B636_382E07C10000__INCLUDED_)
