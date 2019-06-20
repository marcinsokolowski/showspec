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
#include "myparser.h"
#include <string.h>
#include "mystrtable.h"
#include "myfract.h"
#include "basestructs.h"

MyParser::MyParser():mystring()
{
}


MyParser::MyParser(const char* szString) : mystring(szString)
{
	m_pos = 0;	
}

void MyParser::Reset()
{
	m_pos=0;
}


BOOL_T MyParser::IsNumber(char c)
{
	char last='9';
	char first='0';
	return (c>=first && c<=last);
}

BOOL_T MyParser::IsNumerical(const char* str)
{
	int len = strlen(str);
	for(int i=0;i<len;i++){
		if(!IsNumber(str[i]))
			return FALSE;	
	}
	return TRUE;	
}

BOOL_T MyParser::IsWhite(char c)
{
	return (c==' ' || c=='\t');
}

void MyParser::SkipWhite()
{
	while(m_pos<length() && IsWhite((*this)[m_pos]))
		m_pos++;
}

void MyParser::SkipWhiteAndSeps(const char* seps)
{
	while(m_pos<length() && ( IsWhite((*this)[m_pos]) ||
		  strchr(seps,(*this)[m_pos]) ) )
		m_pos++;
}

BOOL_T MyParser::GetVarAndValue( CEnvVar& paramDesc, const char* sep )
{
	return GetVarAndValue( paramDesc.szName, paramDesc.szValue, sep );
}

BOOL_T MyParser::GetVarAndValueFromLine( mystring& varname,mystring& value )
{
	const char* pVarName;
	const char* pVarVal;
	BOOL_T bRet = TRUE;	

	if (pVarName = GetItemUpTo("=")){
		varname = pVarName;

		PlusPlus();
		if (pVarVal = GetItemUpTo("\n")) 
			value = pVarVal;
		else
			value = "";
		return TRUE;
	}
	return FALSE;

}

BOOL_T MyParser::GetVarAndValue(mystring& varname,mystring& value,const char* sep)
{
	const char* pVarName;
	const char* pVarVal;
	BOOL_T bRet = TRUE;	

	if (pVarName = GetItemUpTo(sep)){
		varname = pVarName;

		if (pVarVal = GetItemUpTo("=\n")) 
			value = pVarVal;
		else
			value = "";
		return TRUE;
	}
	return FALSE;
}
	

const char* MyParser::GetItemUpTo(const char* upto,BOOL_T bReadNull)
{
	m_szItem="";
	if (!bReadNull){
		SkipWhiteAndSeps(upto);
	}else{
		if (m_pos<length()){
			if(strchr(upto,(*this)[m_pos]))
				m_pos++;
		}
	}
	if (m_pos<length() && (*this)[m_pos]!='\n'){
		int i=m_pos;
		while(i<length() && (*this)[i]!='\n' && 
			  !strchr(upto,(*this)[i]))
			i++;
		char* pTmp = new char[i-m_pos+1];
		strncpy(pTmp,&((*this)[m_pos]),i-m_pos);
		pTmp[i-m_pos]='\0';
		m_szItem = pTmp;
		delete [] pTmp;
		m_pos = i;
		return (m_szItem.length() ? m_szItem.c_str() : NULL);
	}
	return NULL;
}

const char* MyParser::GetItemUpToString(const char* upto,BOOL_T bReadNull)
{
	int upto_len = strlen( upto );
	m_szItem="";
	if (!bReadNull){
		SkipWhiteAndSeps(upto);
	}else{
		if (m_pos<length()){
			if(strchr(upto,(*this)[m_pos]))
				m_pos++;
		}
	}
	int len = length();
	if (m_pos<(len-upto_len) && (*this)[m_pos]!='\n'){
		const char* curr_pos = &( (*this)[m_pos] );
		const char* ret = strstr( curr_pos , upto );
		if( ret ){
			int l = ( ret - curr_pos + 1 );
			char* tmp = new char[l];
			strncpy( tmp, curr_pos, l );
			m_szItem = tmp;
			delete [] tmp;
			m_pos += (l+upto_len-1);
			return m_szItem.c_str();
		}
	}
	return NULL;
}


const char* MyParser::GetItemUpToWithNull(const char* upto)
{
	return GetItemUpTo(upto,TRUE);
}


const char* MyParser::GetNextItem(const char* upto)
{
	mystring szUpTo;
	szUpTo << " " << "\t" << upto;
	return GetItemUpTo(szUpTo.c_str());
}

const char* MyParser::GetNextItemNew(const char* upto)
{
	return GetItemUpTo( upto );
}



const char* MyParser::GetBetween( char b,char e)
{
	long pos1;
	pos1 = find(b,m_pos);
	if (pos1!=NOT_FOUND){
		unsigned pos2;
		pos2 = find(e,pos1+1);
		if (pos2!=NOT_FOUND){
			m_szItem = substr(pos1+1,pos2-pos1-1).c_str();
			return m_szItem.c_str();
		}
	}
	return NULL;
}

const char* MyParser::GetNextLine(BOOL_T bAddEndLine)
{
	m_szLine="";
	while( m_pos<length() && (*this)[m_pos]!='\n'){
		m_szLine += (*this)[m_pos];
		m_pos++;
	}
	
	// NEW MS - 2009-06-27
	// if added due to error when parsing null terminated strings - 
	// multiply lines strings :
	if( m_pos<=length() ){ 
		if ((*this)[m_pos]=='\n' && bAddEndLine)
			m_szLine += (*this)[m_pos];
	}
	m_pos++;
	return (m_szLine.length() ? m_szLine.c_str() : NULL);
}


void MyParser::trimleft()
{
	m_pos=0;
	BaseString szTmp;
	SkipWhite();
	szTmp = this->substr(m_pos,length()-m_pos);
	(*this) = szTmp.c_str();
}


const char* MyParser::GetItem(int pos)
{
	MyParser parser;
	const char* ret;
	parser = (*this);

	int i=0;
	Reset();
	while(i<=pos){
		ret = GetNextItem();
		if (!ret)
			break;
		i++;
	}
	return m_szItem.c_str();
}

void MyParser::operator++(){ 
	int len = length();
	if (m_pos+1 < len)
		m_pos++;
}

// note this function is different then operator++ : it makes plus plus also
// if at last character of string :
void MyParser::PlusPlus(){ 
	int len = length();
	if (m_pos < len)
		m_pos++;
}

void MyParser::SkipNonNumbers()
{
	mystring szTmp;
	int i = 0;
	while(i<length()){
		char c = (*this)[i];
		if (IsNumber(c) || c=='.')
			szTmp << c;
		i++;
	}
	(*this) = szTmp.c_str();
}
	

int MyParser::GetItems( CMyStrTable& tab, const char* sep/*=","*/ )
{
	const char* item;
	tab.clear();
	while( item = GetNextItem(sep) ){
		tab.Add( item );
	}
	return tab.size();
}

int MyParser::GetItemsNew( CMyStrTable& tab, const char* sep/*=","*/ )
{
	const char* item;
	tab.clear();
	while( item = GetNextItemNew(sep) ){
		tab.Add( item );
	}
	return tab.size();
}



BOOL_T MyParser::GetFraction( const char* szFraction, CFraction& fract  )
{	
	MyParser pars;
	pars = szFraction;

	fract.down = 1;
	if (strstr(szFraction,"/")){
		const char* up = pars.GetItemUpTo("/");
		const char* down = "1";
		down = pars.GetNextItem();
		
		LONG_T _up = (LONG_T)( atof(up) );
		LONG_T _down = (LONG_T)( atof(down) );		
		fract.Set( _up, _down );		
	}else{
		LONG_T _up = (LONG_T)( atof(szFraction) );
		fract.Set( _up, 1 );
	}
	return TRUE;
}

BOOL_T MyParser::ParseSaveArea( const char* szVal, 
											int& x0, int& y0, int& x1, int& y1 )
{
	x0 = 0;
	y0 = 0;
	x1 = 0;
	y1 = 0;
	if( szVal && szVal[0] ){
   	CMyStrTable items;
      MyParser szPars = szVal;
      szPars.GetItems(items);
      if( items.size()>=2 ){
      	if( atol( items[0].c_str() )>0 )
         	x0 = atol( items[0].c_str() );
	      if( atol( items[1].c_str() )>0 )
   	      y0 = atol( items[1].c_str() );
      }
      if( items.size()>=4 ){
      	if( atol( items[2].c_str() )>0 )
         	x1 = atol( items[2].c_str() );
	      if( atol( items[3].c_str() )>0 )
   	      y1 = atol( items[3].c_str() );
      }
		
		return TRUE;
	}
	return FALSE;
}
