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
#ifndef _MY_STRING_TABLE_HH
#define _MY_STRING_TABLE_HH

#include "mytypes.h"
#include "basedefines.h"
#include "mystring.h"
#include <vector>

using namespace std;

BASELIB_EI const char* mygetenv(const char* varname);

class BASELIB_EI  CMyStrTable : public vector<mystring>
{
public:
   CMyStrTable(){};
	CMyStrTable(const CMyStrTable& right);   
   
   CMyStrTable& operator=(const CMyStrTable& right);
   
	void Add(const char* str);
	void AddInt( int value );
	const char* Find(const char* str);
	int FindPos(const char* str);
	int Remove(const char* str);

	mystring& get(LONG_T pos);
	LONG_T GetCount(){ return size(); };
	void Sort();

	void SortByDateFormat( const char* fmt );
	
	void DumpToFile( const char* szFName );
};


#endif
