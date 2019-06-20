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
#ifndef _POINTS_H__
#define _POINTS_H__
#include "mytypes.h"
#include "mystring.h"
#include <vector>


using namespace std;


enum eStarIdentType_T { eStarIdent_OK=0, eStarIdent_ReUsedOld };

class CTransformMatrixInTime;

class CPoint 
{
public :
	int m_FrameIndex;
	double x;
	double y;
	time_t frame_time;
	double m_Value;
	eStarIdentType_T m_IdentType;		
	CPoint(double _x,double _y,int frame_index=0,time_t _time=0,
			 double val=0,eStarIdentType_T identType=eStarIdent_OK);
	CPoint(const CPoint& right);
	CPoint();	

	void Set( LONG_T x, LONG_T y );
	
	CPoint& Add( CPoint& elem1, CPoint& elem2 );
	CPoint& Subtract( CPoint& elem1, CPoint& elem2 );

	static double dist( CPoint& elem1, CPoint& elem2 );
	static double dist( double x1, double y1, double x2, double y2 );
	static const char* get_ident_desc( eStarIdentType_T identType, BOOL_T bNotEnd );	

	static double calc_angle( CPoint& left_limit, CPoint& right_limit );

	static int find_transform( CPoint& f1_s1, CPoint& f1_s2, 
										CPoint& f2_s1, CPoint& f2_s2,
									   int size_x, int size_y,
							         CTransformMatrixInTime& matrix, int nLevel=0 );

	static BOOL_T IsBetween( double x0, double x1, double x );
	static BOOL_T IsBetweenBig( double x0, double x1, double x, double min_dist=3 );
};

class CPointDesc : public CPoint
{
public :
	mystring m_Desc;
	double m_mag;
	int m_sample_size;
	
	CPointDesc(double _x,double _y,const char* szDesc,
				  double mag=-1,int sample_size=2) : CPoint( _x, _y ), m_Desc(szDesc), m_mag(mag), m_sample_size(sample_size) {}
	CPointDesc(const CPointDesc& right) : CPoint( right.x, right.y ), m_Desc(right.m_Desc), m_mag(right.m_mag), m_sample_size(right.m_sample_size) {} 
};

class CRectangle
{
	public:
		CRectangle()
		 :start_x(0),start_y(0),end_x(0),end_y(0),color(0){};
		CRectangle(long _s_x,long _s_y,long _e_x, long _e_y, long col=0 )
		 :start_x(_s_x),start_y(_s_y),end_x(_e_x),end_y(_e_y),color(col){};
		void Init(long _s_x,long _s_y,long _e_x, long _e_y, long col=0 );
		
		long start_x;
		long start_y;
		long end_x;
		long end_y;
		long color;
};

class CLongPoint 
{
public :
	LONG_T x;
	LONG_T y;
	CLongPoint(LONG_T _x,LONG_T _y);
	CLongPoint(const CLongPoint& right);
	CLongPoint();	
};

class CSector 
{
public :
	CPoint m_Begin;
	CPoint m_End;

	// for soring additional information - maybe a little bit not good
	// place for this - but I need it fast ..
	int m_FrameIndex;
	
	CSector(	double x1, double y1, double x2, double y2, int frame_index=0 );	
	CSector( const CPoint& begin, const CPoint& end, int frame_index=0 );
	
	CPoint& BegPoint(){ return m_Begin; }
	CPoint& EndPoint(){ return m_End; } 
};

class CPointList : public vector<CPoint>
{
public :
	CPointList(){};	
	void Add( int x, int y );
	void Add( const CPoint& new_elem );
	void Clear();	
	CPointList& operator+=(CPointList& newElems);
	CPoint* FindPoint( long x, long y );		
	CPoint* FindByFrame( int frame_index );
	
	int ReadFromFile( const char* filename );
};

class CPointDescList : public vector<CPointDesc>
{
public:
	CPointDescList(){};
	void Add( const CPointDesc& new_elem );		
	CPointDesc* FindPoint( long x, long y );	
};

class CSectorList : public vector<CSector>
{
public :
	CSectorList(){};
	void Add( double x1, double y1, double x2, double y2 );
	void Add( const CPoint& begin, const CPoint& end );
	void Add( const CSector& sec );
	void Clear();		
}; 

class CLongList : public vector<LONG_T>
{
public:
	CLongList();
	CLongList(CLongList& newElems);
	void Add(const LONG_T new_elem);		
	void AddUnique(const LONG_T new_elem);								
	void MergeUnique(CLongList& newElems);
	void RemoveAll(){ clear(); }
	LONG_T GetCount(){ return size(); }
	BOOL_T FindPoint(LONG_T pos);
	void Clear();		
	CLongList& operator+=(CLongList& newElems);		
	CLongList& operator=(CLongList& newElems);
	CLongList& operator=(const CLongList& newElems);
	void Sort();
	void RemoveLast();
	
	// debug :
#ifdef _DEBUG
	void dump();
#endif
};

class CVector3D 
{
public :
	CVector3D() : x(0), y(0), z(0) {}

	double x;
	double y;
	double z;

};

#endif
