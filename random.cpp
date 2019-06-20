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
#ifdef _WINDOWS
#include <afxwin.h>
#endif
#include <stdlib.h>
#include "random.h"
#include "mydate.h"
#include "mycmnglobals.h"
#include <math.h>
#include <myfits.h>

BOOL_T CRandom::m_bInitialized=FALSE;
unsigned long CRandom::m_seed=DEFAULT_SEED;

CRandom::CRandom()
{
	Initialize();
}

CRandom::~CRandom()
{
}

void CRandom::Initialize(unsigned long seed)
{
	if (!m_bInitialized){
		m_seed = seed+get_dttm();
		if( gChangeSeed ){
			srand(m_seed);
		}else{
			if( gExternalSeed>0 ){
            srand( gExternalSeed );
         }else{
				srand(seed); // just temporary to make things more predictable ... :-)
			}
		}
		m_bInitialized = TRUE;
	}
}

double  CRandom::GetRandom()
{
	Initialize();
	double r = double(rand())/double(RAND_MAX+1.0); 
	return r;
}

ULONG_T CRandom::GetRandomInteger(ULONG_T lowerbound,ULONG_T upperbound)
{
	Initialize();
	double r = GetRandom(); 
	ULONG_T ret = lowerbound+(ULONG_T)((upperbound - lowerbound)*r);
	return ret;
}

double CRandom::GetGauss( double sigma, double mean )
{
	Initialize();
	
	double mean50 = mean*50;
	double norm=1.0/(sqrt(2*3.1415)*sigma);
	double sigma2 = sigma*sigma;
	int i=0;
	int max_iter=10000;
	while( i<=max_iter ){
		double x = ( GetRandom()-0.5 )*mean50;
		double y = GetRandom()*norm;
		double gauss_value = norm*exp( - ( ((x-mean)*(x-mean))/(2.00*sigma2) ) );

		if( y <= gauss_value ){
			return x;
		}

		i++;
	}

	printf("ERROR : value from gauss distribution could not be obtained !\n");
	return -1000;	
}

double CRandom::GetFastGauss( double sigma, double mean )
{
	return CMyFit::GetGaussFast( sigma, mean );
}