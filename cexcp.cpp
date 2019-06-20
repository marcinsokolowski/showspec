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
#include <stdlib.h>
#include <stdarg.h>
#include "cexcp.h"
#include "myfile.h"
#include "mydate.h"
#include "mymacros.h"
#include "mycmnglobals.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void PrintToTraceFile( const char* fmt,...)
{
		mystring szFmt;
		szFmt << CMyDate::getdate() << " : " << fmt << "\n";
		MyOFile o_file(ASSERT_FILE_NAME,"a+");
		va_list plist;
		va_start(plist,fmt);
		o_file.VPrintf(szFmt.c_str(),plist);
		char buff[1000];
		o_file.Close();
		if (vsprintf(buff,szFmt.c_str(),plist)>=1000){
			o_file.Printf("Cirtical error buffer in Assert to small exiting");
			exit(0);
		}
		printf("%s\n",buff);
}

void AssertNULLFunc(const void* param,BaseString str) 
{
	if (!param){
		BaseString szTmp;
		szTmp << CMyDate::getdate() << " : Unexpected NULL pointer at : " << str.c_str() << "\n";
		printf("%s\n",szTmp.c_str());
		MyOFile o_file(ASSERT_FILE_NAME,"a+");
		o_file.Printf(szTmp.c_str());
		throw CExcBasic(szTmp.c_str());
	}
}


void AssertFunc(BOOL_T bAssert,const char* fmt,...)
{
	if(!bAssert){
		mystring szFmt;
		szFmt << CMyDate::getdate() << " : " << fmt << "\n";
		va_list plist;
      va_start(plist,fmt);

		// pritnf to STDOUT
		char buff[1000];
		if (vsprintf(buff,szFmt.c_str(),plist)>=1000){
			printf("Cirtical error buffer in Assert to small exiting");
			exit(-1);
		}
		printf("%s\n",buff);

		// to FILE :
		MyOFile o_file(ASSERT_FILE_NAME,"a+");
		o_file.VPrintf(szFmt.c_str(),plist);
		o_file.Close();
		exit(-1);
		// throw CExcBasic(buff);
	}
}

CExcBasic::CExcBasic(const char* pMsg, const char* szCode)
: m_Errcode( ERR_UNKNOWN ), m_errno(0)
{
	if(gPrintfLevel>=0){
		_TRACE_PRINTF_0("Exception %s occured due to : %s\n",szCode,pMsg);
		fflush(0);
	}

	m_Msg = pMsg;
	m_szCode = szCode;
}

CExcBasic::CExcBasic(int _errno,const char* pMsg,eERRCODE_T errcode)
{
	if(gPrintfLevel>=0){
		_TRACE_PRINTF_0("Exception errno=%d occured due to : %s\n",_errno,pMsg);
		_TRACE_PRINTF_0("Errno desc : %s\n",strerror(_errno));
		fflush(0);
	}

	m_Msg = pMsg;
	m_Errcode = ERR_UNKNOWN;
	m_errno = _errno;
}
	

CExcBasic::CExcBasic(const char* pMsg,eERRCODE_T errcode)
{
	m_Msg = pMsg;
	m_Errcode = ERR_UNKNOWN;
}


const char* CExcBasic::GetMsg() const
{
	return m_Msg.c_str();
}

eERRCODE_T CExcBasic::GetErrcode() const
{
	return m_Errcode;
}

mystring CMemAlloc::GetReport()
{
	mystring szMsg;
	szMsg << "Allocating : " << (int)m_size << " bytes\n";
	szMsg << "Errno      : " << m_errno;
	// szMsg << "Message    : " << m_Msg;
	return szMsg;
}
