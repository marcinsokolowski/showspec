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
#ifndef _MY_LOCK_H__
#define _MY_LOCK_H__

#include <pthread.h>
#include "mytypes.h"


class CMyMutex
{
public :
	CMyMutex();
	~CMyMutex();

	void Lock();
	void UnLock();
	int IsLocked(){ return is_locked; }
	
protected :	
	pthread_mutex_t mutex;	
	int id;
	int is_locked;
};


class CMutexLock
{
protected:
	CMyMutex* m_pLock;
public:
	CMutexLock( CMyMutex* pLock ): m_pLock(pLock) { 
		m_pLock->Lock(); 
	}
	~CMutexLock(){ m_pLock->UnLock(); }
};


class CMutexUnLock
{
protected:	
   CMyMutex* m_pLock;
   BOOL_T m_bDoUnLock;
public:
	CMutexUnLock( CMyMutex* pLock, BOOL_T bDoUnlock=TRUE ): m_pLock(pLock), m_bDoUnLock(bDoUnlock) {}
	~CMutexUnLock(){ 
		if( m_bDoUnLock && m_pLock ){
			m_pLock->UnLock(); 
		}			
	}
};

#endif
