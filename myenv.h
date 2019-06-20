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
// myenv.h: interface for the CMyEnv class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MYENV_H__9D2077ED_B738_11D5_B636_382E07C10000__INCLUDED_)
#define AFX_MYENV_H__9D2077ED_B738_11D5_B636_382E07C10000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <vector>
#include "mystring.h"
#include "basedefines.h"
#include "mykeytab.h"

using namespace std;

#define ENV_FILE_NAME "env.cfg"


class BASELIB_EI CMyEnv  
{
public:
	CMyEnv();
	virtual ~CMyEnv();
	static const char* mygetenv(const char* name);
	static void GetEnvList(mystring& szEnv);
protected :
	static void InitEnv();
	static BOOL_T m_bInitialized;
	static CKeyTab  m_EnvList;
};

#endif // !defined(AFX_MYENV_H__9D2077ED_B738_11D5_B636_382E07C10000__INCLUDED_)
