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
#ifndef _BaseString_H__
#define _BaseString_H__

/* class BaseString */
#include "basedefines.h"
#include "mytypes.h"
#include <stdio.h>
#include <string.h>

int my_strlen(const char* str);

class BASELIB_EI BaseString {
protected:
	// member variables
	char* m_pchData;
	int m_Size;
	int m_Length;

	// member functions 
	void ReAllocBuffer(int size=DEFAULT_SIZE);
	void Assign(const char* pszString,int len);
	void MakeReplacement(const char* pOld,const char* old,
					     const char* nnew);
public:
	BaseString();
	~BaseString();
	BaseString(long number);
	BaseString(long size,long size2);
	BaseString(const char z);
	BaseString(const char* szString);
	BaseString(const BaseString& szString);
	BaseString operator +(const char* pszString);
	BaseString operator +(const BaseString& szString);
	BaseString& operator +=(const char* pszString);
	BaseString& operator +=(const BaseString& szString);
	BaseString& operator =(int number);
	BaseString& operator =(const char* pszString);
	BaseString& operator =(const BaseString& szString);
	BaseString& operator =(double);
	void AllocBuffer(int size=DEFAULT_SIZE);
	
	int operator ==(const BaseString& szString) const;
	int operator ==(const char* szString) const;
	int operator <(const BaseString& szString) const;
	int operator <(const char* szString) const;
	int operator >(const BaseString& szString) const;
	int operator >(const char* szString) const;

	BaseString& operator<<(const char* str);
	BaseString& operator<<(const BaseString& str);
	BaseString& operator<<(const int n);
	BaseString& operator<<(const long n);
	BaseString& operator<<(INTN_T n);
	BaseString& operator<<(const char z);
	BaseString& operator<<(LONGLONG_T);
	BaseString& operator<<(double);

	// same as << but clears buffer before << 
	BaseString& operator<<=(const char* str);


	char& operator[](int pos);
	const char& operator[](int pos) const;
	inline operator const char*(){	return m_pchData; }
	inline int length() const { return m_Length; } 
	inline const char* c_str() const { return m_pchData; }
	inline const char* c_str(){ return m_pchData; }
	inline char* str(){ return m_pchData; }
	inline void clear() { m_pchData[0]='\0'; }
	int replace(const char* old,const char* nnew,BOOL_T bAll);
	int replace(const char* old,const char* nnew);
	void replace(int pos,int len,const char* nnew);
	void Reset();
	long find(const char* txt,long pos=0);
	long find(char z,long pos=0);
	BaseString substr(long pos,long len=0);
	BOOL_T empty() const { return ((!m_pchData) || (!m_pchData[0]));}
};


#endif
