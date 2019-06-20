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
#include "mystring.h"
#include "cexcp.h"
#include "myenv.h"
#include <stdio.h>
#include <stdarg.h>
#ifdef _UNIX
#include <ctype.h>
#endif
#include "myfile.h"

#include "myparser.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

const char* mygetenv(const char* varname)
{
	return CMyEnv::mygetenv(varname);
}


mystring::~mystring()
{
}
	
mystring::mystring(const mystring& str)
{
	(*this) = str;
}

mystring::mystring(int n)
{
	char buffer[BUFF_SIZE];
	if (sprintf(buffer,"%d",n)>=BUFF_SIZE)
		AssertNULL(NULL);
	(*this) << buffer;
}

mystring::mystring(int size,int size2)
: BaseString(size,size2) 
{

}

/* mystring::mystring(double db)
{
	(*this) = db;	
}*/
		
void mystring::replace_char( char* ptr, char old,char nnew)
{
	int i,len = strlen(ptr);
	i=0;
	while(i<len){
		if (ptr[i]==old)
			ptr[i]=nnew;
		i++;
	}
}

mystring& mystring::replace_char(char old,char nnew)
{
	int i,len = length();
	i=0;
	while(i<len){
		if ((*this)[i]==old)
			(*this)[i]=nnew;
		i++;
	}
	return (*this);
}



void mystring::splitpath(BaseString& drv,BaseString& dir,
					     BaseString& fname,BaseString& ext) const
{
#ifndef _UNIX
	char szDrive[_MAX_DRIVE];
	char szDir[_MAX_DIR];
	char szFName[_MAX_FNAME];
	char szExt[_MAX_EXT];
	_splitpath(c_str(),szDrive,szDir,szFName,szExt);
	drv = szDrive;
	dir = szDir;
	fname = szFName;
	ext = szExt;	
#else
	drv = "";
	int pos=length()-1;
	mystring szTmp,szTmp2;
	while(pos>=0 && m_pchData[pos]!='\\' && m_pchData[pos]!='/'){
		szTmp2 << m_pchData[pos];
		pos--;		
	}
	for(int kk=szTmp2.length()-1;kk>=0;kk--){
		szTmp<< szTmp2[kk];
	}
	int i;
	for(i=0;i<pos;i++){
		dir << m_pchData[i];
	}
	if(i>0)
		dir << "/";
	const char* ptr;
	if (ptr = strstr(szTmp.c_str(),".")){
		fname = szTmp.substr(0,ptr-szTmp.c_str());
		ext = ptr+1;		
	}else{
		dir << "/" << szTmp;
	}		
#endif
}


const char* mystring::env2str()
{
	unsigned pos=0;
	int len=0;
	unsigned pos1=0;
	while( (pos1=BaseString::find("$(",pos1+len)) !=NOT_FOUND){
		unsigned pos2=find(")",pos1+2);
		if (pos2!=NOT_FOUND)
		{
			BaseString szVarName;			
			szVarName = substr(pos1+2,pos2-pos1-2);
			const char* szEnvValue;
			szEnvValue = mygetenv(szVarName.c_str());
			if(szEnvValue!=NULL){
				BaseString::replace(pos1,pos2-pos1+1,szEnvValue);
				len = strlen(szEnvValue);
			}else{
				len = pos2 - pos1 + 1;
			}
		}	
	}
	return c_str();
}




mystring& mystring::SkipChar( char z )
{
	mystring ret;
   register int len = length();
   for(register int i=0;i<len;i++){
		if( (*this)[i]!=z ){
			ret << (*this)[i];
		}
	}
	(*this)=ret;
	return (*this);
}



mystring& mystring::SkipApostrophs()
{
	mystring ret;
   register int len = length();
   for(register int i=0;i<len;i++){
		if( (*this)[i]!='\'' ){
			ret << (*this)[i];
		}
	}
	(*this)=ret;
	return (*this);
}

mystring mystring::getupper() const
{
	mystring ret;
	ret = (*this);
	int len=ret.length();
	for(register int i=0;i<len;i++)
		ret[i] = toupper(ret[i]);
	return ret;
}


mystring mystring::getlower() const
{
	mystring ret;
	ret = (*this);
	int len=ret.length();
	for(int i=0;i<len;i++)
		ret[i] = tolower(ret[i]);
	return ret;
}
	
	
mystring mystring::fromgreat()
{
	mystring ret;
	ret = (*this);
	if (ret.length()){
		ret[0] = toupper(ret[0]);
	}
	return ret;
}
	

mystring mystring::double2string( double x, const char* fmt )
{
	char szTmp[1024];
	sprintf(szTmp,fmt,x);
	mystring szRet = szTmp;
	return szRet;
}


const char* mystring::Fgets()
{
	ReAllocBuffer(2048);
	char* ptr = ::fgets(m_pchData, 2048, stdin);
	return m_pchData;
}

/*BOOL_T mystring::operator==(const mystring& right)
{
	return (strcmp(m_pchData,right.m_pchData)==0);
}

BOOL_T mystring::operator>(const mystring& right)
{
	return (strcmp(m_pchData,right.m_pchData)==1);
}

BOOL_T mystring::operator<(const mystring& right)
{
	return (strcmp(m_pchData,right.m_pchData)==-1);
}*/



mystring mystring::get_number( mystring& szStr )
{
	mystring szOut;
	int len = szStr.length();
	for(int i=0;i<len;i++){
		if(MyParser::IsNumber( szStr[i] ))
			szOut << szStr[i];		
	}
	return szOut;
}

char mystring::get_first_non_white( const char* szString )
{
	register int len = strlen(szString);
	if(len>0){
		for(register int i=0;i<len;i++){
			if(szString[i]!=' ' && szString[i]!='\t'){
				return szString[i];
			}		
		}
	}

	return '\0';
}

char mystring::get_last_non_white( const char* szString )
{
	register int len = strlen(szString);
	if(len>0){
		for(register int i=(len-1);i>=0;i--){
			if(szString[i]!=' ' && szString[i]!='\t'){
				return szString[i];
			}		
		}
	}

	return '\0';
}

BOOL_T mystring::IsInApostrophs( const char* szString )
{
	if( get_first_non_white(szString)=='\'' && get_last_non_white(szString)=='\'')
		return TRUE;
	return FALSE;
}


mystring&  mystring::TrimApostrophs( char z )
{
	mystring ret="";
	int i=0;
	int len0 =strlen( m_pchData );
	while( i<len0 && m_pchData[i] && (m_pchData[i]==' ' || m_pchData[i]=='\t' || m_pchData[i]==z) ){
		i++;
	}
	ret = (char*)(&(m_pchData[i]));

	int len = strlen(ret.m_pchData);
	int pos=len-1;
	while(pos>=0 && (ret.m_pchData[pos]==' ' || ret.m_pchData[pos]=='\t' || ret.m_pchData[pos]==z )){
		pos--;
	}
	Assert(pos>=-1 && pos<len,"Error in TrimApostro for string %s",m_pchData);
	ret.m_pchData[pos+1]='\0';

	(*this) = ret;
	return (*this);
}

mystring mystring::TrimApostro( char z )
{
	mystring ret="";
	int i=0;
	int len0 =strlen( m_pchData );
	while( i<len0 && m_pchData[i] && (m_pchData[i]==' ' || m_pchData[i]=='\t' || m_pchData[i]==z) ){
		i++;
	}
	ret = (char*)(&(m_pchData[i]));

	int len = strlen(ret.m_pchData);
	int pos=len-1;
	while(pos>=0 && (ret.m_pchData[pos]==' ' || ret.m_pchData[pos]=='\t' || ret.m_pchData[pos]==z )){
		pos--;
	}
	Assert(pos>=-1 && pos<len,"Error in TrimApostro for string %s",m_pchData);
	ret.m_pchData[pos+1]='\0';

	return ret;
}

void mystring::DumpToFile( const char* szFileName, const char* mode )
{
	MyOFile out( szFileName , mode );
	out.Printf("%s\n",m_pchData);	
}

void mystring::SkipTrailSpaces()
{
	SkipTrailSpaces( m_pchData );
	m_Length = strlen( m_pchData );
}

void mystring::SkipTrailSpaces( char* ptr )
{
	int len = strlen(ptr);
	int i=(len-1);
	while( i>=0 && ptr[i]==' ' ){
		ptr[i]='\0';
		i--;
	}
}

int mystring::get_item( const char* pLine, int start_pos, int end_pos, mystring& szItem )
{
	szItem = "";
	int len = strlen(pLine);
	int pos=start_pos;
	while(pos<len && pos<=end_pos ){
		szItem << pLine[pos];
		pos++;
	}
	return strlen(szItem.c_str());
}

int mystring::get_item( const char* pLine, int start_pos, int end_pos, 
								mystring& szItem, int& end_out, int& last_space )
{
	szItem = "";
	int len = strlen(pLine);
	int pos=start_pos;
	int nGood=0;
	BOOL_T bLastGood=FALSE;

	while(pos<len && ( pos<=end_pos || bLastGood ) ){
		if( pLine[pos] != ' ' && pLine[pos] != '\t' ){
			end_out = pos;
			nGood++;
			bLastGood = TRUE;
		}else{
			last_space = pos;
			bLastGood = FALSE;
		}
		szItem << pLine[pos];
		pos++;
	}

	if( nGood == 0 )
		end_out = pos-1;

	return nGood;
}



int mystring::get_item_safe( const char* pLine, int start_pos, int end_pos, mystring& szItem )
{
	szItem = "";
	int len = strlen(pLine);
	int pos=start_pos;
	
	if( pos>0 ){
		if( pLine[pos-1] != ' ' ){
			int i=pos-1;
			while(i>=0 && pLine[i]!=' '){
				szItem << pLine[i];
				i--;
			}
		}
	}


	BOOL_T bWasEnd=FALSE;	
	while(pos<len && ( pos<=end_pos || pLine[pos]!=' ' ) ){
		if( pos>end_pos && bWasEnd ){
			break;
		}
		if( pLine[pos] == ' ' || pLine[pos] == '\t' ){
			bWasEnd = TRUE;
		}else{
			bWasEnd = FALSE;
		}
		szItem << pLine[pos];
		pos++;
	}
	return strlen(szItem.c_str());
}


const char* mystring::skip_white( const char* line )
{
	int len = strlen(line);
	int add=0;
	while( add<len && ( line[add]==' ' || line[add]=='\t' )){
		add++;
	}
	return (line+add);
}
