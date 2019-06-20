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
#ifndef _MYKEYTAB_H__
#define _MYKEYTAB_H__

#include "mystring.h"
#include "basestructs.h"
#include <vector>


using namespace std;



class CKeyTab : public vector<CEnvVar>
{
public:
	CKeyTab(){};
	~CKeyTab(){};
	void Add( const CEnvVar& elem );		
	void Add(const char* key,const char* val, const char* comment="");
	void Add(const char* key,int val,const char* comment="");
	void Add(const char* key,double val, const char* comment="");
	void Delete( const char* key);
	CEnvVar* Find(const char* key);
	CEnvVar* Find(const char* key, int pos);

	void Set( const char* key,const char* val, const char* comment="" );
	void Set( const char* key,double val, const char* comment="" );
	void Set( const char* key,int val, const char* comment="" );
};


#endif
