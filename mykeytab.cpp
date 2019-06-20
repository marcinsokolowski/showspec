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
#include "mykeytab.h"


void CKeyTab::Add(const char* key,const char* val, const char* comment)
{
	CEnvVar tmp;
	tmp.szName = key;
	tmp.szValue = val;
	tmp.szComment = comment;
	push_back(tmp);
}

void CKeyTab::Add(const char* key,int val, const char* comment)
{
	CEnvVar tmp;
	tmp.szName = key;
	tmp.szValue << val;
	tmp.szComment = comment;
	push_back(tmp);
}

void CKeyTab::Add( const char* key,double val, const char* comment )
{
	CEnvVar tmp;
	tmp.szName = key;
	tmp.szComment = comment;
	char szTmp[64];
	sprintf(szTmp,"%.8f",val);
	tmp.szValue = szTmp;
	push_back(tmp);
}

void CKeyTab::Add( const CEnvVar& elem )
{
	push_back( elem );
}

void CKeyTab::Delete( const char* key )
{
	vector<CEnvVar>::iterator i;
	for(i=begin();i!=end();i++){
		if(strcmp(i->szName.c_str(),key)==0){
			erase( i );	
			return;
		}
	}
}

CEnvVar*  CKeyTab::Find(const char* key)
{
	vector<CEnvVar>::iterator i;
	for(i=begin();i!=end();i++){
		if(strcmp(i->szName.c_str(),key)==0)
			return &(*i);
	}	
	return NULL;
}

CEnvVar*  CKeyTab::Find(const char* key, int pos)
{
	int count=size();
	if( pos>=0 && pos<count ){
		for(int i=pos;i<count;i++){
			if( strcmp( (*this)[i].szName.c_str(),key)==0 ){
				return &( (*this)[i] );
			}			
		}		
	}
	return CKeyTab::Find( key );
}


void CKeyTab::Set( const char* key,double val, const char* comment )
{
	char szTmp[64];
   sprintf(szTmp,"%.8f",val);

	Set( key, szTmp, comment );
}

void CKeyTab::Set( const char* key,int val, const char* comment )
{
	char szTmp[64];
   sprintf(szTmp,"%d",val);

	Set( key, szTmp, comment );
}

void CKeyTab::Set( const char* key,const char* val, const char* comment )
{
	CEnvVar* pKey = Find( key );
	if(pKey){
		pKey->szValue = val;
		if(comment && comment[0]){
			pKey->szComment = comment;
		}
	}else{
		Add( key, val, comment );
	}
}