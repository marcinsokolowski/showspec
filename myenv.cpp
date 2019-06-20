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
// myenv.cpp: implementation of the CMyEnv class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _UNIX
#include "stdafx.h"
#endif
#include "myenv.h"
#include "myfile.h"
#include "myparser.h"
#include "basestructs.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

BOOL_T CMyEnv::m_bInitialized=FALSE;
CKeyTab CMyEnv::m_EnvList;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMyEnv::CMyEnv()
{

}

CMyEnv::~CMyEnv()
{

}

void CMyEnv::InitEnv()
{
	if(!m_bInitialized){
		MyIFile envfile(ENV_FILE_NAME,FALSE);
		if (envfile.IsOpened()){
			const char* pLine;
			while(pLine = envfile.GetLine(TRUE)){
				MyParser parser = pLine;
				CEnvVar tmp;
				parser.GetVarAndValue(tmp.szName,tmp.szValue);
				m_EnvList.push_back(tmp);
			}
		}
		m_bInitialized = TRUE;
		vector<CEnvVar>::iterator i;
		for(i=m_EnvList.begin();i!=m_EnvList.end();i++){
			mystring szTmp = i->szValue.c_str();
			szTmp.env2str();
			i->szValue = szTmp.c_str();
		}
	}
}


const char* CMyEnv::mygetenv(const char* name)
{
	vector<CEnvVar>::iterator i;
	InitEnv();
	for(i=m_EnvList.begin();i!=m_EnvList.end();i++){
		if (STRCASECMP(i->szName.c_str(),name)==0){
			return i->szValue.c_str();
		}
	}
	const char* szRet = ::getenv(name);
	return szRet;
}

void CMyEnv::GetEnvList(mystring& szEnv) 
{ 
	vector<CEnvVar>::iterator i;
	InitEnv();
	szEnv = "";
	for(i=m_EnvList.begin();i!=m_EnvList.end();i++){
		szEnv << i->szName << "=" << i->szValue << "\n";
	}
}
