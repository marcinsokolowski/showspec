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
#ifndef _BASESTRUCTS_H__
#define _BASESTRUCTS_H__ 

#include "mystring.h"

class  CEnvVar {
public:
	mystring szName;
	mystring szValue;
	mystring szComment;

	CEnvVar(){};		
	CEnvVar( const char* sName, const char* sValue, const char* sComment="" );
	CEnvVar(const CEnvVar& right);
	CEnvVar& operator=(const CEnvVar& right);


	friend bool operator > ( const CEnvVar& left, const CEnvVar& right ){
		return (strcmp(left.szName.c_str(),right.szName.c_str())>0);
	}
	friend bool operator < ( const CEnvVar& left, const CEnvVar& right ){
		return (strcmp(left.szName.c_str(),right.szName.c_str())<0);
	}
	friend bool operator == ( const CEnvVar& left, const CEnvVar& right ){
		return (strcmp(left.szName.c_str(),right.szName.c_str())==0 &&
		        strcmp(left.szValue.c_str(),right.szValue.c_str())==0);
	}
};

struct sCommand
{
   int command;
   mystring szParam1;
	mystring szParam2;
};
         

/*bool CompareEnvVar( const CEnvVar& left, const CEnvVar& right )
{
	return (strcmp(left.szName.c_str(),right.szName.c_str())==1);
}*/


#endif
