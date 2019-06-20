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
#ifndef _C_EXCP_H__
#define _C_EXCP_H__
#include "mytypes.h"
#include "basedefines.h"
#include "mystring.h"

#define ASSERT_FILE_NAME "$(TRACEDIR)/lib.trace"

class BASELIB_EI CExcBasic 
{
protected:
	BaseString m_Msg;
	BaseString m_szCode;
	eERRCODE_T m_Errcode;
	int m_errno;
public:	
	CExcBasic(const char* pMsg, const char* szCode);
	CExcBasic(const char* pMsg=NULL,eERRCODE_T errcode=ERR_UNKNOWN);
	CExcBasic(int _errno,const char* pMsg=NULL,eERRCODE_T errcode=ERR_UNKNOWN);	
	const char* GetMsg() const;
	eERRCODE_T GetErrcode() const;
	int GetErrno() const{ return m_errno;}
};

class BASELIB_EI CExcFile : public CExcBasic
{
public:
	CExcFile(int _errno,const char* pMsg=NULL,eERRCODE_T errcode=ERR_UNKNOWN)
		: CExcBasic(_errno,pMsg,errcode) {};

};

class BASELIB_EI CMemAlloc
{
   long m_size;	
   int m_errno;
	public :
		CMemAlloc(int _errno,int size) : m_errno(_errno),m_size(size){}
	   mystring GetReport();
};

void PrintToTraceFile( const char* fmt,...);
void BASELIB_EI AssertNULLFunc(const void* param,BaseString str);
void BASELIB_EI AssertFunc(BOOL_T bAssert,const char* fmt,...);

#define AssertNULL(param) AssertNULLFunc(param,BaseString()<<__FILE__<<":"<<__LINE__)
#define Assert AssertFunc


#endif
