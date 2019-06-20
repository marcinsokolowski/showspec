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
#ifndef _UNIX
#include <direct.h>
#else
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#define DIR_MODE (S_IRWXU | S_IRWXG | S_ISGID | S_IROTH | S_IXOTH)
#define FILE_MODE ( S_IRWXU | S_IRWXG | S_IROTH )
#endif
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include "myfile.h"
#include "cexcp.h"
#include "mystrtable.h"
#include "mysafekeytab.h"
#include "myparser.h"
#include "basestructs.h"
#include "basedefines.h"
#include "mydate.h"

// some useful globals :
/*BASELIB_EI const char* getbasename(const mystring& name)
{
	static mystring basename;
	int i=0;
	basename.clear();
	while(name[i]!='.' && name[i]!='\0'){
		basename += name[i];
		i++;
	}
	basename += '\0';
	return basename.c_str();
}*/

BASELIB_EI const char* getbasename_new(const mystring& name,mystring& out)
{
	int i=0;
	out.clear();
	while(name[i]!='.' && name[i]!='\0'){
		out += name[i];
		i++;
	}
	out += '\0';
	return out.c_str();
}


mystring getfname_with_ext( const mystring& name)
{
	mystring szDrv,szDir,szFName,szExt;
	mystring szTmp = name;
	szTmp.splitpath( szDrv,szDir,szFName,szExt );
	szFName << "." << szExt.c_str();
	return szFName;
}

mystring getfname( const mystring& name)
{
	mystring szDrv,szDir,szFName,szExt;
	mystring szTmp = name;
	szTmp.splitpath( szDrv,szDir,szFName,szExt );
	return szFName;
}

mystring getfname( const mystring& name, mystring& szExt)
{
	mystring szDrv,szDir,szFName;
	mystring szTmp = name;
	szExt = "";
	szTmp.splitpath( szDrv,szDir,szFName,szExt );
	return szFName;
}


MyFile::MyFile(const char* filename,const char* attr) :
m_pFile(NULL), m_pBuffer(NULL), m_Size(FILE_BUFF_SIZE), m_bOutBuffEnabled(FALSE),
m_nBuffUsed(0)
{
	if (filename && attr && strlen(filename) && strlen(attr))
		Open( filename , attr );
}

MyFile::~MyFile()
{
	if (m_pBuffer)
		delete [] m_pBuffer;
	Close();
}

BOOL_T MyFile::Rewind()
{
	if(m_pFile){
		rewind( m_pFile );	
		return TRUE;
	}
	return FALSE;
}

BOOL_T MyFile::OpenFile(const char* filename,BOOL_T bExcp)
{
	return MyFile::Open( filename, "r", bExcp );
}

void MyFile::Clear( const char* filename )
{
	MyOFile out( filename, "w" );	
}

BOOL_T MyFile::Open(const char* filename,const char* attr,BOOL_T bExcp)
{
	m_FileName = filename;
	m_Attr = attr;
	m_FileName.env2str();
	if (!m_FileName.empty()){
		if(m_pFile)
			Close();
		PrepareFile( m_FileName.c_str() );

		//if( strlen( m_FileName.c_str() ) ){
		//	printf("Cannot open file with empty file name !!!\n");
		//	return FALSE;
		//}

		m_pFile = fopen(m_FileName.c_str(),attr);
		if (!m_pFile && bExcp){
			mystring szMsg;
			szMsg << "Could not open file (encosed in :) :" << m_FileName.c_str() << ": in mode :" << attr << ":";
			PrintToTraceFile( szMsg.c_str() );
			throw CExcFile(errno,szMsg.c_str());
		}
	}
	return (m_pFile!=NULL);
}

void MyFile::Flush()
{
	if (m_pFile)
		fflush( m_pFile );
}

BOOL_T MyFile::GetLine( mystring& szLine )
{
	szLine="";

	char* line = GetLine( TRUE );
	if( line ){
		szLine = line;
		return TRUE;		
	}
	return FALSE;
}

char* MyFile::GetLine(BOOL_T bNoNewLine/*=FALSE*/)
{
	AssertNULL(m_pFile);
	if (!m_pBuffer){
		m_pBuffer = new char[m_Size];
	}
	m_pBuffer[0]='\0';	
	char* ret = fgets(m_pBuffer,m_Size,m_pFile);
	if(strlen(m_pBuffer)==(m_Size-1) && m_pBuffer[m_Size-1]!='\n'){
		mystring szFullLine;
		szFullLine << m_pBuffer;
		// Assert(FALSE,"Buffer to small to read line of file %s",m_FileName.c_str());		
		while(strlen(m_pBuffer)==(m_Size-1) && m_pBuffer[m_Size-1]!='\n'){
			char* ptr = fgets(m_pBuffer,m_Size,m_pFile);
			szFullLine << m_pBuffer;
		}
		m_Size = strlen(szFullLine.c_str())+1;
		char* pNewBuffer = new char[m_Size];
		strcpy(pNewBuffer,szFullLine.c_str());
		delete [] m_pBuffer;
		m_pBuffer = pNewBuffer;
	}

//	if( strlen(m_pBuffer)<=0 ){
//	  printf("CRITICAL ERROR !!!! - MEMORY ACCESS BUG !!!\n");
//	}

	if( m_pBuffer ){
    	int buff_len=strlen(m_pBuffer);
    	if( buff_len>0 ){ // otherwise already it is : m_pBuffer[0]='\0' 
    	  if ( bNoNewLine && m_pBuffer[buff_len-1]=='\n')
          m_pBuffer[buff_len-1]='\0';
      }
   }
	return (ret ? m_pBuffer : ret);
}

const char* MyFile::GetBaseName()
{
	if (m_szBaseName.empty()){
		int i=0;
		while(m_FileName[i]!='.'){
			m_szBaseName += m_FileName[i];
			i++;
		}
	}
	return m_szBaseName.c_str();
}



void MyFile::Close()
{
	if (m_pFile)
		fclose(m_pFile);
	m_pFile = NULL;
}

void MyFile::DumpBuffer()
{
	fprintf(m_pFile,"%s",m_szOutBuffer);
	m_nBuffUsed = 0;
}

void MyFile::PrintfNow(const char* fmt,...)
{
	if(strcmp(m_Attr.c_str(),"a+")==0 || strcmp(m_Attr.c_str(),"w")==0) 
		CreateDir(m_FileName.c_str());
	if (!m_pFile)Open(m_FileName.c_str(),m_Attr.c_str());
	va_list plist;
	va_start(plist,fmt);
	vfprintf(m_pFile,fmt,plist);
	va_end(plist);
}

void MyFile::Printf(const char* fmt,...)
{
	if(strcmp(m_Attr.c_str(),"a+")==0 || strcmp(m_Attr.c_str(),"w")==0) 
		CreateDir(m_FileName.c_str());
	if (!m_pFile)Open(m_FileName.c_str(),m_Attr.c_str());
	va_list plist;
	va_start(plist,fmt);
	if(!m_bOutBuffEnabled){
		vfprintf(m_pFile,fmt,plist);
	}else{
		char* pos = (m_szOutBuffer+m_nBuffUsed);
		m_nBuffUsed += vsprintf(pos,fmt,plist);
		if(m_nBuffUsed>=MAX_OUT_BUFFER_SIZE){
			Assert(FALSE,"Internal buffer size exceed !!! - disable file buffering");
		}
		if(m_nBuffUsed>=DUMP_WHEN)
			DumpBuffer();
	}
	va_end(plist);
}

void MyFile::VPrintf(const char* fmt,va_list args)
{
	if (!m_pFile)Open(m_FileName.c_str(),m_Attr.c_str());
	vfprintf(m_pFile,fmt,args);	
}


BOOL_T MyFile::Delete( const char* fname)
{
	unlink( fname );
	return TRUE;
}

BOOL_T MyFile::DoesFileExist(const char* fname)
{
	BOOL_T bRet=FALSE;
	if(fname && fname[0] ){
		mystring szFile=fname;
		szFile.env2str();

		if( access( szFile.c_str(), F_OK ) == 0 ) {
			bRet = TRUE;
		}
	}
	return bRet;
}


BOOL_T MyFile::IsReadable(const char* fname)
{
	if(DoesFileExist(fname)){
		if( access( fname, R_OK ) == 0 ) {
			return TRUE;
		}
	}
	return FALSE;
}

void MyFile::CreateDir(const char* path)
{
	mystring szPath = path;
	szPath.env2str();
#ifndef _UNIX
	szPath.replace_char('/','\\');
#endif
	mystring ext,dir,drv,fname;
	mystring pth;
	mode_t mode=DIR_MODE;

	szPath.splitpath(drv,dir,fname,ext);
	int pos=-1;
	mystring szSep="\\";
#ifdef _UNIX
	szSep = "/";	
	pos++;
#endif
	while((pos=dir.find(szSep.c_str(),pos+1))!=NOT_FOUND){
		pth = drv;
		pth << dir.substr(0,pos);
#ifndef _UNIX
		int res = _mkdir(pth.c_str());
#else
		int res = mkdir(pth.c_str(),mode);
#endif
		if (res==0 || (res==-1 && errno==EEXIST))
		{ /* OK */ }
		else
		{		
			mystring szError;
			szError << "Could not create dir : " << path;	
			throw CExcFile(errno,szError.c_str());
		}

	}

}

void MyFile::CreateDirBase(const char* path)
{
	mystring szPath = path;
	szPath.env2str();
#ifndef _UNIX
	szPath.replace_char('/','\\');
#endif
	mode_t mode=DIR_MODE;

	int pos=-1;
	const char* szSep="\\";
#ifdef _UNIX
	szSep = "/";	
	pos++;
	if(!strstr(szPath.c_str(),"/"))
		szPath << "/";
#endif
	while((pos=szPath.find(szSep,pos+1))!=NOT_FOUND){
		mystring pth;

		pth << szPath.substr(0,pos);
#ifndef _UNIX
		int res = _mkdir(pth.c_str());
#else
		int res = mkdir(pth.c_str(),mode);
#endif
		if (res==0 || (res==-1 && errno==EEXIST))
		{ /* OK */ }
		else
		{		
			mystring szError;
			szError << "Could not create dir : " << path;	
			throw CExcFile(errno,szError.c_str());
		}

	}

}


MyIFile::MyIFile(const char* filename,BOOL_T bExcp)
{ 
	if (filename && filename[0]){
		Open (filename,"r",bExcp);
	} 
}


MyOFile::MyOFile(const char* filename,const char* attr)
{
	if (filename && filename[0]){
		CreateDir(filename);
		Open( filename, attr );
	}
}

void MyOFile::PrepareFile( const char* filename )
{
	if (filename && filename[0]){
		CreateDir(filename);
	}	
}

BOOL_T MyOFile::OpenFile(const char* filename,BOOL_T bExcp)
{	
	BOOL_T bRet=FALSE;
	if (filename && filename[0]){
		CreateDir(filename);
		bRet = Open( filename, "w" );
	}
	return bRet;
}

MyOFile& MyOFile::operator <<(const char* str)
{
	Printf("%s",str);
	return (*this);
}

MyOFile& MyOFile::operator <<(int num)
{
	Printf("%d",num);
	return (*this);
}

CListFile::CListFile()
{
	m_List = new CMyStrTable();
}

CListFile::CListFile(const char* filename):MyIFile(filename,FALSE)
{
	Assert(m_pFile!=NULL,"Could not open file %s",filename);
	m_List = new CMyStrTable();
	ReadListFile( filename );
}

int CListFile::ReadListFile( const char* filename )
{
	const char* pLine;
	if(!IsOpened()){
		Open( filename );
	}
	m_List->clear();
	while(pLine = GetLine(TRUE)){
		if(strlen(pLine) && pLine[0]!='#')
			m_List->Add(pLine);
	}
	Close();	
	return m_List->size();
}

int CListFile::UpdateListFile( const char* filename )
{
	const char* pLine;
	if(!IsOpened()){
		Open( filename );
	}
	int new_count=0;
	while( pLine = GetLine(TRUE) ){
		if(strlen(pLine) && pLine[0]!='#'){
			if( !m_List->Find( pLine ) ){
				m_List->Add(pLine);
				new_count++;
			}
		}
	}
	Close();	
	return new_count;	
}


CListFile::~CListFile()
{
	if(m_List)
		delete m_List;
}

const char* CListFile::operator[](int i)
{
	Assert((i>=0 && i<m_List->GetCount()),"CListFile::operator[] index %d out of range [0-%d]",i,(m_List->GetCount()-1));
	return (*m_List)[i].c_str();
}

LONG_T CListFile::GetCount()
{
	if(!m_List)
		return 0;
	return m_List->GetCount();
}


CDescFile::CDescFile()
{
	m_KeyTab = new CSafeKeyTab();	
}

CDescFile::CDescFile(const char* filename)
{
	m_KeyTab = new CSafeKeyTab();
	Init( filename );
}

void CDescFile::Init(const char* filename)
{
	if(m_KeyTab->GetCount())
		return;
	MyIFile infile(filename);
   const char* pLine;
   while(pLine = infile.GetLine(TRUE)){
		if (strlen(pLine)==0 || pLine[0]=='\n')
	      continue;
		// skip comments :
		if (pLine[0]=='#' || strncmp(pLine,"//",2)==0)
			continue;

		MyParser parser = pLine;
		CEnvVar tmp;
		parser.GetVarAndValue(tmp.szName,tmp.szValue," \t");
		tmp.szValue.env2str();

		if(strlen(tmp.szName.c_str())==0 && strlen(tmp.szValue.c_str())==0)
			continue;

		m_KeyTab->Add(tmp);		
   }
   infile.Close();	
	
}

CDescFile::~CDescFile()
{
	if(m_KeyTab)
		delete m_KeyTab;
}

LONG_T CDescFile::GetCount()
{
	return m_KeyTab->GetCount();
}


const char* MyFile::GetCWD( mystring& szDir )
{
	char* buff = NULL;
 	const char* pDir = NULL;

	szDir = "";

	int init_size = 2048;

	// printf("%d, %s\n",errno,strerror(errno));fflush(0);
	errno = 0;
	while(errno==ERANGE || init_size==2048){
		buff = new char[init_size];
		pDir = getcwd( buff, init_size );
		init_size = init_size*2;
		// printf("%d, %s\n",errno,strerror(errno));fflush(0);
	}

	if(pDir){
		szDir = pDir;
		delete [] buff;
		return szDir.c_str();
	}

	delete [] buff;
	return NULL;
}



int MyFile::GetFileSize( const char* filename )
{
	if( !DoesFileExist( filename ) ){
		return -1;
	}

	struct stat buf;
	stat( filename, &buf );
	return (int)buf.st_size;
}

time_t MyFile::GetLastModif( const char* filename )
{
   if( !DoesFileExist( filename ) ){
      return 0;
   }
   
   struct stat buf;
   stat( filename, &buf );
   return (int)buf.st_mtime;                  
}

int MyFile::WriteFile( const char* filename, char* ptr, int size )
{
	int ret=0;
	if( filename && strlen(filename)){
		CreateDir( filename );
		int fd = open( filename, O_CREAT | O_WRONLY,  FILE_MODE  );		
		if( fd>0 ){
			ret=write( fd, ptr, size );
			close(fd);
			if(ret!=size){
				printf("error while writing to file : %s\n",filename);
				printf("written %d bytes, while input was %d\n",ret,size);
				ret=0;	
			}
		}else{
			printf("could not open file : %s\n",filename);
			printf("due to error : %s\n",strerror(errno));
		}
	}
	return ret;
}

int MyFile::ReadFile( const char* filename, char*& ptr )
{
	ptr = NULL;
	if( filename && strlen(filename)){
		if( DoesFileExist( filename ) ){
			int size = GetFileSize( filename );
			if(size>=0){
				ptr = new char[size];
				mystring szFileName=filename;
				szFileName.env2str();
				int fd = open( szFileName.c_str(), O_RDONLY );
				if( fd>0 ){
					int ret=read( fd, ptr, size );
					close(fd);
					if(ret!=size){
						printf("error while reading file : %s\n",szFileName.c_str());
						delete [] ptr;
						return 0;
					}else{
						return ret;
					}
				}
			}
		}
	}
	return 0;
}

void MyFile::init_day_file_counter( const char* filename, int& counter )
{
	mystring szFile = filename;
	szFile.env2str();
	BOOL_T bFromFile=FALSE;
	if( DoesFileExist( szFile.c_str() ) ){
		mystring szDT;					
		get_night_date_local( szDT );
		MyIFile in( szFile.c_str() );
      const char* pLine = in.GetLine();
      if(pLine){
      	MyParser pars = pLine;
         CMyStrTable items;
         pars.GetItems(items);
         if( items.size()>=2 ){
	         if(strcmp(items[0].c_str(),szDT)==0){
   	      	counter = atol( items[1].c_str() );
      	      bFromFile=TRUE;
         	}
	      }
		}
	}
	if( !bFromFile )
		counter = 0;
}

void MyFile::save_day_file_counter( const char* filename, int counter )
{
	mystring szFile = filename;
   szFile.env2str();
   mystring szDT;
   get_night_date_local( szDT );


   MyOFile out( szFile.c_str() );
   out.Printf("%s %d\n",szDT.c_str(),counter);
}