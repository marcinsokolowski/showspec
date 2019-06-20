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
#ifndef _MYPARSER_H__
#define _MYPARSER_H__

#include "mystring.h"
#include "basedefines.h"
#include "mytypes.h"

class CMyStrTable;
class CFraction;
class CEnvVar;

class BASELIB_EI MyParser : public mystring
{
public:
	MyParser();
	MyParser(const char* szString);
	const char* GetNextItem(const char* upto=",");
	const char* GetNextItemNew(const char* upto=",");
	const char* GetNextLine( BOOL_T bAddEndLine=TRUE );
	const char* GetBetween( char b,char e);
	const char* GetItem(int pos);

	BOOL_T GetVarAndValue(mystring& varname,mystring& value, const char* sep="=");
	BOOL_T GetVarAndValue( CEnvVar& paramDesc, const char* sep="=" );
	BOOL_T GetVarAndValueFromLine( mystring& varname,mystring& value );		
	
	void SkipWhite();
	void SkipWhiteAndSeps(const char* seps=",");
	void SkipNonNumbers();
	const char* GetItemUpTo(const char* upto,BOOL_T bReadNull=FALSE);
	const char* GetItemUpToString(const char* upto,BOOL_T bReadNull=FALSE);
	const char* GetItemUpToWithNull(const char* upto);
	int GetItems( CMyStrTable& tab, const char* sep="," );
	int GetItemsNew( CMyStrTable& tab, const char* sep="," );
	void Reset();
	void trimleft();
	void operator++();
	void PlusPlus();
	const char* GetCurrent(){ return &( (*this)[m_pos] ); }

	static BOOL_T IsWhite(char c);
	static BOOL_T IsNumber(char c);
	static BOOL_T IsNumerical(const char* str);
	
	static BOOL_T GetFraction( const char* szFraction, CFraction& fract  );
	static BOOL_T ParseSaveArea( const char* szVal, int& x0, int& y0, int& x1, int& y1 );
protected:
	int m_pos;
	int m_linepos;
	mystring m_szItem;
	mystring m_szLine;
};

#endif
