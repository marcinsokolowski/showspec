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
#ifndef _CTYPES_H__
#define _CTYPES_H__

// TRACK TYPES ETC :
enum eTrackCheckType_T { eNormalTrack=0, ePlaneTrack, eSingleCamTrack,
                         eRejIfMoreTrack, eTrackOnSumedFrame,
                         eAddTo3Points };
                        
// wheather track contains events found in TLE file                          
enum eTrackEventsType_T { eTrackEvents_NoTLE=0, eTrackEvents_HasTLE=1 };

// fitted line type : streight line, parabola etc ...
enum eTrackFitLineType_T { eTrackFitSLine=0, eTrackFitParabola=1 };
                         

enum eHistoVariableType_T { eCVariableInfo=0, eHistoIfMore, eChi2ToOld, 
									 eHistoCoicRADEC, eHistoCoicPix, eMinChi2_On3Points,
									 eHistoAfterTv, eHistoFrameEvents, eMinDistStar,
									 eMinSatDist, 
									 eHistoVXRatioToOld, eHistoVYRatioToOld,
									 eHistoVXRatioTo3, eHistoVYRatioTo3,
									 eChi2OnCurrFrame, eChi2OnCurrFrame3Points,
									 eHistoVXRatioToOldSum, eHistoVYRatioToOldSum,
									 eHistoVXRatioTo3Sum, eHistoVYRatioTo3Sum,
									 eChi2ToOldSum,
									 eHistoIfMoreAfterCoic,
									 eRXall, eRYall, eVXvsVX, eVYvsVY };


enum eDriverReverseImage_T { eReverseImageNone=0, eReverseImageFull,
									  eReverseImageHor, eReverseImageVert, 
									  eReverseImageHorVert };


enum eLaplaceType_T { eSinglePoint=0, eTwoPoints, eFourPoints, eFivePoints, 
							 eFivePlusFourMin, eNineEightFive, eNineEightSevenVeryBig,
                      eEightFour, eEightTen, eFiveEight, eRawS, 
                      eFourTwelve, eFourTwelveFar, eClusterLaplace,
                      eLastEnumElem };
                      
                      
enum eFitType_T { eFitHorizontalLine=0, eFitLine, eFitGauss };


typedef unsigned long ULONG_T;
typedef unsigned char BOOL_T;
typedef long DATE_T;//date in format YYYYMMDD
typedef long TIME_T;//time in format HHMMSS
typedef long DATETIME_T; // dttm format
typedef long MYDATE_T;//date in format YYYYMMDD
typedef long YEAR_T;//year YYYY
typedef ULONG_T INTN_T;
typedef long MONEY_T;
typedef long PERCENT_T;
typedef long LONG_T;
typedef long long LONGLONG_T;

#define FALSE 0
#define TRUE  1
#define INTN_NULL 0

#define BZERO(obj,class_name) memset(obj,'\0',sizeof(class_name))


// enums :
enum eERRCODE_T {ERR_UNKNOWN=0,ERR_DBOPEN,ERR_SELECT};
//enum eFIELDTYPE_T {F_NUMBER=0,F_STRING,F_MONEY,F_DATE,F_ENUM};
enum eMYBOOL_T { MYFALSE=0, MYTRUE, MYERROR };
enum eFIELDTYPE_T {	MY_STRING=0,MY_NUMBER,MY_DATE,MY_ENUM,
			 		MY_DATETIME,MY_PASSWORD,MY_MONEY};
enum eOPERATION_T { eBorrow=0,eOrder };


// constans :
#define DEFAULT_DATE_FORMAT "%H:%M:%S %d/%m/%Y"
#define NOT_DEFINED -1

int is_number(const char* string);
int is_digit( char znak );

#endif
