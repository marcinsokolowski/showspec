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
#ifndef _MY_UTIL_H__
#define _MY_UTIL_H__

#include "mytypes.h"
#include "mystring.h"
#include <math.h>
#include <vector>
#include "mymacros.h"
#include "basestructs.h"

using namespace std;

const char* GetOK( BOOL_T bOK );

float non_zero( float val, float val_default );

void my_printf_now( const char* msg );

double MAX_func( double x, double y );
double MIN_func( double x, double y );


void my_insertion_sort( LONG_T* tab, LONG_T cnt );
void my_qsort( LONG_T* tab, LONG_T cnt );
void my_qsort_int( int* tab, LONG_T cnt );
int get_median( LONG_T* tab, LONG_T cnt );
int get_no_stars_median( LONG_T* tab, LONG_T cnt );


void my_sort_float( double* ftab, long cnt );

BOOL_T my_table_overlap_check( LONG_T* tab1, LONG_T cnt1, 
			 						    LONG_T* tab2, LONG_T cnt2 );

BOOL_T my_sorted_table_overlap_check( LONG_T* tab1, LONG_T cnt1, 
			 				 				     LONG_T* tab2, LONG_T cnt2 );

LONGLONG_T my_find_max_value( LONGLONG_T* table, LONG_T cnt );

LONG_T my_find_max_value_long( LONG_T* table, LONG_T cnt );

double find_max_value( double* values, int cnt, int& pos );

BOOL_T check_sort( LONG_T* tab, LONG_T size );

void DumpTable( LONG_T* tab, LONG_T size );

/*inline LONG_T my_round(double x)
{
	register int x_l = (int)floor(x);

	register double h = (x-x_l);
	if(h>0.5)
		return (x_l+1);
	else
		return (x_l);
}*/

inline int my_round(double x)
{
	if( x>=0 ){
		return (int)(x+0.50);
	}else{
		return (int)(x-0.50);
	}
}

const char* my_stristr( const char* in_string, const char* substr );

LONG_T long_sum( LONG_T* tab, LONG_T cnt, LONG_T start_from=0 );

int find_value( LONG_T* tab, int tab_size, LONG_T val );

int qfind_value( LONG_T* tab, int tab_size, LONG_T val );

int find_value( vector<int>& tab, int value );

int find_string( const char** tab, const char* szValue );

int find_min_value( int* tab, int count, int& min_pos );

int do_print_usage( int argc, char* argv[], int min_argc_req );

int is_default( const char* param );

mystring get_uname();

inline int mystrlen( const char* str )
{ if(!str)return 0; else return strlen(str); }

inline void DecZero( int& value )
{
	value--;
	if(value<0)
		value=0;
}

BOOL_T compare_double( double l, double r );

int safe_atol( const char* szValue );
double safe_atof( const char* szValue );
const char* safe_string( const char* szValue );

mystring get_param_value( const char* szParam );
BOOL_T check_param( const char* param, const char* name, mystring& szValue );
BOOL_T check_param( const char* param, const char* name, int& nValue );
BOOL_T check_param( const char* param, const char* name, double& fValue );
BOOL_T check_param( const char* param, const char* name, float& fValue );
BOOL_T ParseCommand( const char* szCmd, sCommand& cmd );

int exec_cmd( const char* szCmd, char* szOutput, int size );

// average / RMS 
double calc_average_rms( double* values, int cnt, double& rms );
double calc_rms( int mean, LONG_T* values, int cnt );

// cluster :
void CopyCluster( LONG_T* tab_out, LONG_T* tab_in, int cnt );

void print_line();

void dump_table( LONG_T* tab ,int cnt );

int check_zero( int value );

void my_calc_mean_and_sigma( int* tab, int count, double& mean, double& rms, int min_val=-10000000, int max_val=10000000 );
void my_calc_mean_and_sigma_elem( ELEM_TYPE* tab, int count, double& mean, double& rms, int min_val=-10000000, int max_val=10000000 );

#endif
