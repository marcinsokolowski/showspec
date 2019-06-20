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
#ifndef _BASELIB_DEFS_H__
#define _BASELIB_DEFS_H__

#include "mytypes.h"
#include <stdlib.h>


// defines
#define DEFAULT_SIZE 100
#define NOT_FOUND -1
#define BUFF_SIZE 100

#define FILE_BUFF_SIZE 2048

// #define ELEM_TYPE short
#define ELEM_TYPE unsigned short

// NEW BIG_ELEM_TYPE changed long -> int !!! 20050930 !!!
// in case of changes, template in myfitslib must be modified also 
// #define BIG_ELEM_TYPE LONGLONG_T
#define BIG_ELEM_TYPE int

// samples with negative values are kept in :
#define ELEM_SAMPLE_TYPE short

#define DEFAULT_CFG_FILE "ccd.cfg"


#define MAX_PIXELS_IN_LAPLACE 20
#define MAX_LAPLACE_DEFINED 20 

#ifdef _UNIX
	#define BASELIB_EI
   #define STRCASECMP strcasecmp
#else
   #ifdef _BASELIB_DEF_
	    #define BASELIB_EI __declspec(dllexport)
   #else
	   #define BASELIB_EI __declspec(dllimport)	
   #endif
   #define STRCASECMP stricmp
#endif


#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define MAX(x,y) ((x)>=(y) ? (x) : (y))


// constants :
// #define PI 3.14159265358979323846

// structures :
class ImageStat {
public :
	ImageStat();
	~ImageStat();
	ImageStat( const ImageStat& right );	
	void AllocDistribTab( long maxValue );
	LONG_T GetMaxOccurValue();

   double Average;
   double RMS;
   double Min;
   double MinNonZero;
   double Max;
   double Sum;
   double MeanFit;
   double SigmaFit;

	double MaxX;
	double MaxY;		

	long* m_pValDistrib;
	long  m_ValCount;
	
	BOOL_T bFitOK;
};

inline BOOL_T CheckRange( int x, int y, int xSize, int ySize )
{
	return (x>=0 && x<xSize && y>=0 && y<ySize );
}

inline const char* GetResultString( int res ){
	if( res>0 ){
		return "OK";
	}else{
		return "FAILED";
	}
}

const char* GetFlipDesc( eDriverReverseImage_T fliptype );

#endif
