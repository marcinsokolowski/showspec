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
#ifndef _MY_STRING_HH
#define _MY_STRING_HH

#include "mytypes.h"
#include "basedefines.h"
#include "basestring.h"


BASELIB_EI const char* mygetenv(const char* varname);

class BASELIB_EI mystring : public BaseString
{
public:
	mystring(){};
	virtual ~mystring();
	mystring(const char* str):BaseString(str){};
	mystring(int n);
	mystring(const mystring& str);
	mystring(int size,int size2);
	// mystring(double db);

	/*BOOL_T operator==(const mystring& right);
	BOOL_T operator>(const mystring& right);		
	BOOL_T operator<(const mystring& right);		*/

	static mystring double2string( double x, const char* fmt="%.2f" );

	static mystring get_number( mystring& szStr );
	static int get_item( const char* pLine, int start_pos, int end_pos, mystring& szItem );
	static int get_item( const char* pLine, int start_pos, int end_pos, mystring& szItem, 
								int& end_out, int& last_space );
	static int get_item_safe( const char* pLine, int start_pos, int end_pos, mystring& szItem );
	static const char* skip_white( const char* line );

	mystring getupper() const;
	mystring getlower() const;
	mystring fromgreat();
	mystring& replace_char(char old,char nnew);
	static void replace_char( char* ptr, char old,char nnew);
	
	void splitpath(BaseString& drv,BaseString& dir,
		 	       BaseString& fname,BaseString& ext) const;

	mystring& SkipApostrophs();
	mystring& SkipChar( char z );

	inline int Strcmp( const char* str ){ return strcmp( m_pchData, str ); }

	const char* env2str();
	
	static char get_first_non_white( const char* szString );
	static char get_last_non_white( const char* szString );

	const char* Fgets();	
	
	mystring TrimApostro( char z='\'' );
	mystring&  TrimApostrophs( char z='\'' );
	static BOOL_T IsInApostrophs( const char* szString );

	void SkipTrailSpaces();
	static void SkipTrailSpaces( char* ptr );

	// saving string to file :
	void DumpToFile( const char* szFileName, const char* mode="a+" );
};



#endif
