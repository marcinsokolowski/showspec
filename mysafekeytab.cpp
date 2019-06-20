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
#include <math.h>
#include "basestructs.h"
#include "mykeytab.h"
#include "myfile.h"
#include "myparser.h"
#include "mysafekeytab.h"
#include "mystrtable.h"

CKeyTab* m_Table;

CKeyTab& CSafeKeyTab::getKeyTab()
{
	return (*m_Table);
}

CSafeKeyTab::CSafeKeyTab()
{
	m_Table = new CKeyTab();
}

CSafeKeyTab::CSafeKeyTab(const CSafeKeyTab& right )
{
	m_Table = new CKeyTab();
	(*this) = right;
}

CSafeKeyTab& CSafeKeyTab::operator=(const CSafeKeyTab& right )
{
	(*m_Table) = (*right.m_Table);
	return (*this);
}

CSafeKeyTab::~CSafeKeyTab()
{
	delete m_Table;
}


void CSafeKeyTab::Add( const char* key, const char* val,const char* comment)
{
	m_Table->Add(key,val,comment);
}

void CSafeKeyTab::Add( const char* key, int val, const char* comment )
{
	m_Table->Add(key,val,comment);
}

void CSafeKeyTab::Add( const char* key, double val, const char* comment )
{
	m_Table->Add(key,val,comment);
}


void CSafeKeyTab::Add( const CEnvVar& newelem )
{
	m_Table->push_back( newelem );	
}

void CSafeKeyTab::Set( const char* key, double val, const char* comment )
{
	char szTmp[64];
//   sprintf(szTmp,"%.8f",val);
	if(fabs(val) < 0.001 || val > 9999999.){
		sprintf(szTmp,"%12.8e",val);
	}else{
		sprintf(szTmp,"%10.8f",val);
	}	                    

   Set( key, szTmp, comment );	
}

void CSafeKeyTab::Set( const char* key, int val, const char* comment )
{
	char szTmp[64];
   sprintf(szTmp,"%d",val);

   Set( key, szTmp, comment );	
}

void CSafeKeyTab::Set( const char* key, const char* val, const char* comment )
{
	CEnvVar* pKey = Find( key );
	if(pKey){
		pKey->szValue = val;
		if( comment && comment[0] ){
			pKey->szComment = comment;
		}
	}else{
		Add( key, val, comment );
	}
}

void CSafeKeyTab::Delete( const char* key )
{
	m_Table->Delete( key );
}

CEnvVar*  CSafeKeyTab::Find(const char* key)
{
	return m_Table->Find(key);
}

const char* CSafeKeyTab::getKeyVal( const char* key )
{
	CEnvVar* pVar = m_Table->Find(key);
	if(pVar){
		return pVar->szValue.c_str();
	}
	return NULL;
}

void CSafeKeyTab::Clear()
{
	m_Table->clear();
}

long CSafeKeyTab::GetCount()
{
	return m_Table->size();
}

CEnvVar& CSafeKeyTab::operator[](long i){
	return (*m_Table)[i];
}


BOOL_T CSafeKeyTab::ReadFromFile( const char* fname )
{
	mystring szName = fname;
	szName.env2str();
	if(!MyFile::DoesFileExist( szName.c_str() )){
		return FALSE;
	}
	const char* pLine;
	MyIFile in( szName.c_str() );
	while(pLine=in.GetLine()){
		MyParser pars=pLine;
		CMyStrTable items;
		pars.GetItems(items);
		if(items.size()>=2){
			CEnvVar tmp( items[0].c_str(), items[1].c_str() );
			if(items.size()>=3){
				tmp.szComment = items[2];
			}
			m_Table->push_back( tmp );			
		}
	}
	return TRUE;
}

BOOL_T CSafeKeyTab::Get( const char* key, int& val )
{
	CEnvVar* pVar = m_Table->Find(key);
	if( pVar ){
		val = atol( (pVar->szValue).c_str() );
		return TRUE;
	}
	return FALSE;
}

BOOL_T CSafeKeyTab::Get( const char* key, double& val )
{
	CEnvVar* pVar = m_Table->Find(key);
	if( pVar ){
		val = atof( (pVar->szValue).c_str() );
		return TRUE;
	}
	return FALSE;
}
