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
#include <stdio.h>
#include <string.h>
#include "basestring.h"
#include "cexcp.h"



#ifdef _WINDOWS
#ifdef _DEBUG
#include <afxwin.h>
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif

int my_strlen(const char* str)
{
	return (str ? strlen(str) : 0);
}

BaseString::~BaseString()
{
	delete [] m_pchData;
}


BaseString::BaseString() : m_Size(-1), m_Length(0), m_pchData(NULL)
{	
	AllocBuffer();
}

BaseString::BaseString(const BaseString& szString): m_Size(-1), m_Length(0), m_pchData(NULL)
{
	AllocBuffer();
	(*this) = szString;
}
	

BaseString::BaseString(const char* szString) : m_Size(-1), m_pchData(NULL)
{
	int len=my_strlen(szString);
	AllocBuffer(len>0 ? len : DEFAULT_SIZE);
	if (len)
		strcpy(m_pchData,szString);
	else
		m_pchData[0] = '\0';
	m_Length = len;
}

BaseString::BaseString(const char z)
:m_Size(-1), m_Length(0), m_pchData(NULL)
{
	AllocBuffer();
	sprintf(m_pchData,"%c",z);	
	m_Length=2;
}
	
BaseString::BaseString(long number)
:m_Size(-1), m_Length(0), m_pchData(NULL)
{
	char buffer[BUFF_SIZE];
	if (sprintf(buffer,"%d",(int)number)>=BUFF_SIZE)
		AssertNULL(NULL);
	(*this) << buffer;
}

BaseString::BaseString(long size,long size2)
:m_Size(-1), m_Length(0), m_pchData(NULL)
{
	AllocBuffer( size );
}
	
void BaseString::AllocBuffer(int size/*=DEFAULT_SIZE*/)
{
	if (size>=m_Size){
		char* pTmp=m_pchData;
		m_Size = size+1;
		m_pchData = new char[m_Size];	
		if (my_strlen(pTmp))
			strcpy(m_pchData,pTmp);
		else
			m_pchData[0]='\0';
		if (pTmp)
			delete [] pTmp;
	}
}

void BaseString::ReAllocBuffer(int size/*=DEFAULT_SIZE*/)
{
	if (size>=m_Size){
		m_Size = size+1;
		delete [] m_pchData;
		m_pchData = new char[m_Size];	
	}
}

char& BaseString::operator[](int pos)
{
	Assert(pos>=0 && pos<=m_Length,"Wrong index");
	return m_pchData[pos];
}

const char& BaseString::operator[](int pos) const
{
	Assert(pos>=0 && pos<=m_Length,"Wrong index");
	return m_pchData[pos];
}
	
BaseString& BaseString::operator +=(const char* pszString)
{
	// NEW - just to be safe, NULL is ignored 
	if( pszString ){
		int len=my_strlen(pszString);
		AllocBuffer(m_Length+len);
		if(pszString)
			strcat(m_pchData,pszString);
		m_Length += len;
	}
	return (*this);
}

BaseString& BaseString::operator +=(const BaseString& szString)
{
	AllocBuffer(m_Length+szString.length());
	strcat(m_pchData,szString.m_pchData);
	m_Length += szString.length();
	return (*this);
}

void BaseString::Assign(const char* pszString,int len)
{
	m_Length = len;
	ReAllocBuffer(m_Length);
	strcpy(m_pchData,pszString);	
}


BaseString& BaseString::operator =(const char* pszString)
{
	Assign( pszString, my_strlen(pszString));
	return (*this);
}

BaseString& BaseString::operator =(const BaseString& szString)
{
	Assign(szString.m_pchData,szString.length());
	return (*this);
}

BaseString& BaseString::operator=(int number)
{
	clear();
	(*this) << number;
	return (*this);
}
	
BaseString BaseString::operator+(const char* pszString)
{
	BaseString tmp = (*this);
	tmp += pszString;
	return tmp;
}

BaseString BaseString::operator+(const BaseString& szString)
{
	BaseString tmp = (*this);
	tmp += szString;
	return tmp;
}
	
	

void BaseString::MakeReplacement(const char* pOld,
							   const char* old,
							   const char* nnew)
{
	int lenOld = strlen(old);
	int lenNew = strlen(nnew);
	char* pTmp=new char[ m_Length - lenOld + lenNew+1];
	strncpy( pTmp , m_pchData , pOld-m_pchData);
	pTmp[pOld-m_pchData]='\0';
	strcat(pTmp,nnew);
	strcat(pTmp,pOld+lenOld);
	Assign( pTmp, m_Length - lenOld + lenNew);
	delete [] pTmp;
}

int BaseString::replace(const char* old,const char* nnew)
{
	return replace(old,nnew,TRUE);
}
	

int BaseString::replace(const char* old,const char* nnew,BOOL_T bAll)
{
	char* ptr = m_pchData;
	char* pOld;
	int count=0;
	while(pOld = strstr(ptr,old)){
		MakeReplacement(pOld,old,nnew);
		ptr = m_pchData + strlen(nnew);
		count++;
		if (!bAll)
			break;
	}
	return count;
}

void BaseString::replace(int pos,int len,const char* nnew)
{
	if(m_pchData){
		Assert(pos>=0 && pos<m_Length,"Position out of range");
		Assert(len>=0 && pos+len-1<m_Length,"Position range exceeded");
		m_pchData[pos] = '\0';
		BaseString szBegin = m_pchData;
		BaseString szEnd = m_pchData+pos+len;
		
		// calculating length of new string :
		int begin_len = szBegin.length();
		int end_len = szEnd.length();
		int nnew_len = strlen(nnew);		            		
		m_Length = begin_len + nnew_len + end_len; // lenght of new string

		// always use AllocBuffer not new char[] , because it also updates m_Size !!! 
		// which is later used to check if ReAlloc is needed, the problem with option -cat was 
		// that gCCDParams.m_szASASStarCatalog variable was not assign any value in constructor, so m_Size=101 buffer was allocated
		// but later there was :
		// 1) env2str which changed m_Lenght to 41 and allocated 41 bytes to m_pchData buffer, but left m_Size=101 which was absured, because m_pchData had only 41 bytes !!!
		// 2) assigned of 51 characters string, which called AllocBuffer but exited because m_Size=101 > 51 so no allocation was needed ! , which caused core dump at exit
		//    becase memory was written in not allocated region ( 51 > 41 ) !!! 
		// 
		// CONCLUSION : always call AllocBuffer which makes it in a wise way , do not do it manually with m_pchData allocation and m_Lenght assignment !!!
		AllocBuffer( m_Length ); // allocates space for m_Length string ( + end character )

		strcpy(m_pchData,szBegin.c_str());
		strcat(m_pchData,nnew);
		m_pchData[begin_len+nnew_len]='\0';
		strcat(m_pchData,szEnd.c_str());
		m_pchData[m_Length]='\0';
	}
}
	
long BaseString::find(char z,long pos)
{
	BaseString szTmp = z;
	return find(szTmp.c_str(),pos);
}
	

long BaseString::find(const char* txt,long pos)
{
	Assert(pos>=0,"Position out of range");
	if(pos>=m_Length)
		return NOT_FOUND;
	long ret = NOT_FOUND;
	const char* pStartPos = m_pchData + pos;
	if(pStartPos && pStartPos[0]){
		const char* ptr = strstr(pStartPos,txt);
		if(ptr){
			ret = (ptr-m_pchData);
		}
	}
	return ret;
}

BaseString BaseString::substr(long pos,long len)
{
	Assert(pos>=0 && pos<m_Length,"Position out of range");
	Assert(len>=0,"Position range exceeded");
	BaseString tmp = (*this);
	if(len)
		tmp.m_pchData[pos+len]='\0';
	BaseString szRet = tmp.m_pchData+pos;
	return szRet;
}

int BaseString::operator==(const BaseString& szString) const
{
	return (strcmp(m_pchData,szString.m_pchData)==0);
}

int BaseString::operator==(const char* szString) const
{
	return (strcmp(m_pchData,szString)==0);
}

int BaseString::operator<(const BaseString& szString) const
{
	return (strcmp(m_pchData,szString.m_pchData)<0);
}

int BaseString::operator<(const char* szString) const
{
	return (strcmp(m_pchData,szString)<0);
}

int BaseString::operator>(const BaseString& szString) const
{
	return (strcmp(m_pchData,szString.m_pchData)>0);
}

int BaseString::operator>(const char* szString) const
{
	return (strcmp(m_pchData,szString)>0);
}

BaseString& BaseString::operator<<=(const char* str)
{
	clear();
	return ( (*this)<<str );
}

BaseString& BaseString::operator<<(const char* str){
	(*this)+=str;
	return (*this);
}

BaseString& BaseString::operator<<(const BaseString& str){
	(*this)+=str;
	return (*this);
}	

BaseString& BaseString::operator<<(const int n){
	BaseString tmp((long)n);
	(*this) += tmp.c_str();
	return (*this);
}

BaseString& BaseString::operator<<(const long n){
	BaseString tmp((long)n);
	(*this) += tmp.c_str();
	return (*this);
}

BaseString& BaseString::operator<<(INTN_T n)
{
	BaseString tmp((long)n);
	(*this) += tmp.c_str();
	return (*this);
}

BaseString& BaseString::operator<<(const char z)
{
	char tmp[2];
	sprintf(tmp,"%c",z);
	(*this) << tmp;
	return (*this);
}


BaseString& BaseString::operator<<(LONGLONG_T num)
{
	char buf[100];
	sprintf(buf,"%ld",(long int)num);
	(*this) << buf;
	return (*this);
}

BaseString& BaseString::operator<<(double num)
{
	char buf[256];

	// cannot be like 20 - all aplication have long values :
	if( sprintf(buf,"%.2f",num)>=256 ){
		printf("ERROR in function BaseString::operator= , buffer size 256 exceeded, num=%.2f\n",num);
		exit(0);
	}

	(*this) << buf;
	return (*this);
}

BaseString& BaseString::operator =(double num)
{
	char buf[256];

	// cannot be like 20 - all aplication have long values :
	if( sprintf(buf,"%.2f",num)>=256 )
	{
		printf("ERROR in function BaseString::operator= , buffer size 256 exceeded, num=%.2f\n",num);
		exit(0);
	}
	(*this) = buf;
	return (*this);
}
	

	


