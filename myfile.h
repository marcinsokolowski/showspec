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
#ifndef _MYFILE_H__
#define _MYFILE_H__

#include <stdio.h>
#include <stdarg.h>
#include "mystring.h"
#include "mytypes.h"
#include "basedefines.h"

#define MAX_OUT_BUFFER_SIZE 20000
#define DUMP_WHEN           10000


class CMyStrTable;
class CSafeKeyTab;

// BASELIB_EI const char* getbasename(const mystring& name);
BASELIB_EI const char* getbasename_new(const mystring& name,mystring& out);

BASELIB_EI mystring getfname( const mystring& name);

mystring getfname( const mystring& name, mystring& szExt);

BASELIB_EI mystring getfname_with_ext( const mystring& name);

//const char* BASELIB_EI getbasename(const string& name);

class BASELIB_EI MyFile 
{
protected:
	BOOL_T m_bOutBuffEnabled;
	long m_nBuffUsed;
	char m_szOutBuffer[MAX_OUT_BUFFER_SIZE];
public :
	MyFile(const char* filename=NULL,const char* attr="r");
	~MyFile();
	virtual BOOL_T OpenFile(const char* filename=NULL,BOOL_T bExcp=TRUE);
	virtual void PrepareFile( const char* filename ){};
	BOOL_T Open(const char* filename=NULL,const char* attr="r",BOOL_T bExcp=TRUE);
	char* GetLine(BOOL_T bNoNewLine=FALSE);
	BOOL_T GetLine(mystring& szLine);
	void Close();
	const char* GetBaseName();
	void DumpBuffer();
	void PrintfNow(const char* fmt,...);
	void Printf(const char* fmt,...);
	void VPrintf(const char* fmt,va_list args);
	static void CreateDir(const char* path);
	static void CreateDirBase(const char* path);		
	static BOOL_T DoesFileExist(const char* fname);
	static BOOL_T Delete( const char* fname);
	static BOOL_T IsReadable(const char* fname);
	BOOL_T IsOpened() const { return (m_pFile!=NULL); };
	void Flush();
	const char* GetFileName() { return m_FileName.c_str(); }
	BOOL_T Rewind();
	
	static const char* GetCWD( mystring& szDir );
	static int GetFileSize( const char* filename );
	static time_t GetLastModif( const char* filename );
	static int ReadFile( const char* filename, char*& ptr );
	static int WriteFile( const char* filename, char* ptr, int size );
	static void Clear( const char* filename );

	// file counters :
	static void init_day_file_counter( const char* filename, int& counter );
	static void save_day_file_counter( const char* filename, int counter );
protected:
	FILE* m_pFile;
	char* m_pBuffer;
	int m_Size;
	mystring m_FileName;
	mystring m_Attr;
	mystring m_szBaseName;
};


class BASELIB_EI MyOFile : public MyFile
{
protected:
public:
	MyOFile(const char* filename=NULL,const char* attr="w");		
	virtual BOOL_T OpenFile(const char* filename=NULL,BOOL_T bExcp=TRUE);
	virtual void PrepareFile( const char* filename );
	MyOFile& operator <<(const char* str);
	MyOFile& operator <<(int num);
};

class BASELIB_EI MyIFile : public MyFile
{
public:
	MyIFile(const char* filename=NULL,BOOL_T bExcp=TRUE);
};

class BASELIB_EI CListFile : public MyIFile
{
protected:
	CMyStrTable* m_List;
public :
	CListFile();
	CListFile(const char* filename);	
	~CListFile();
	int ReadListFile( const char* filename );
	LONG_T GetCount();
	inline CMyStrTable& GetListTable() { return (*m_List); }
	const char* operator[](int i);

	int UpdateListFile( const char* filename );
};

class BASELIB_EI CDescFile
{
protected:
	CSafeKeyTab* m_KeyTab;
public :
	CDescFile();
	CDescFile(const char* filename);	
	void Init(const char* filename);
	~CDescFile();
	LONG_T GetCount();
	
	CSafeKeyTab& GetDescTab(){ return (*m_KeyTab); }
};



#endif
