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
#ifndef _MY_FITS_H__
#define _MY_FITS_H__

/* Changes Log 
	2010-03-12 : instead of using CalcChi2 I use CalcChi2_Dist which 
					takes into account fact that x,y error is not only in dY but 
					also in dX !
*/

#include <mytypes.h>
#include <math.h>
#include <stdlib.h>

extern double gMinVeloToCheck;
extern double gMaxVeloRelativeError;
extern int gMaxFrameDistFor3PointSeed;

#define DEFAULT_MAX_DIFF_TIME 3600
#define DEFAULT_CCD_IDX_FOR_CHECK_COND 0

// enum eHistoType_T { eHistoTypeUnknown=0, eChi2Of3Events, eChi2Of3SumEvents };

enum eCheckVelocityReason_T { eCheckVeloDefault=0, eCheckVeloGoodChi2=1, eCheckVeloTryNewPoint=2 };

struct sEventDesc
{
	double x;
	double y;
	int frame;
	int timeUT;
};

void sort_event_list( sEventDesc* list, int cnt );

void show_points( const char* cmt, sEventDesc* list, int cnt );

class CMyFit
{
 protected:

public:
	CMyFit();
	~CMyFit();

	 static BOOL_T (*m_pFitGauss)( double minValue, double& rms, double& center,
	 										 double &max, int binNo, double binWidth,
	 										 int* pCountTab, int& nstep, int first_bin ); 
	 static BOOL_T (*m_pFillFunc)( void* event_info, int type, int ccd_idx, void* event_info2  );

	// very valnorable for initital parameters - especially 
	// when to small number of zero points in Gauss tail
	static BOOL_T FitGauss( double minValue, double& rms, double& center, 
									double &max, int binNo, double binWidth, 
									int* pCountTab, int& nstep,
									int first_bin=0);

	static BOOL_T FitParabola( double* x_values, double* y_values, int count,
	                           double& a, double& b, double& c, double& chi2 );									

	static double CalcDist2FromLine2Par( double a, double b,
													 double* x, double* y, int cnt );

	// Dist2 - means that this is distance squared ( dist^2 )  
	// to have distance result must be sqrt ( square rooted )														 
	static double CalcDist2FromLine( double a, double b, double c,
	                                 double x, double y );
	                                 
	// DISTANCE(POINT,PARABOLA) :
	// WARNING : currently it is chi2 !!! not real distance - maybe should be corrected  ...
	static double CalcDistFromParabola( double a, double b, double c, 
													double x, double y );

	static double CalcDist2FromLine2Par( double a, double b,
	                                     double x, double y );

	static double CalcChi2( double* x_values, double* y_values, int cnt,
	                        double a, double b, int exceptPos=-1 );
	                         
	// calculates Chi2 as sum of distances of points from a streight line 	                         
	// it is in the same power as CalcChi2 ( both distance^2 - without sqrt )
	static double CalcChi2_Dist( double* x_values, double* y_values, int cnt,
	                        double a, double b, double c,
	                        int exceptPos=-1 );
	                         

	static double CalcMaxChi2( double* x_values, double* y_values, int cnt,
   	                 double a, double b, int& pos  );

	static double FitLineChi2( double* x_values, double* y_values, int cnt,
										double& a, double& b, double& c );

	static BOOL_T FitLineHorizontal( double* x_values, double* y_values, int cnt, 
												double& c );

	static BOOL_T FitLine( double* x_values, double* y_values, int cnt,
								  double& a, double& b, double& c,
								  int exceptPos=-1 );

	static double getLineChi2( double x, double y, double a, double b, double c );

	static void RejectPoint( double* x_values, double* y_values, int cnt, int pos );

	static BOOL_T FindPointsOnLine( double* x_values, double* y_values, int cnt,
											  double max_chi2, 
											  double* line_x, double* line_y, int& cnt_on_line,
											  double& a, double& b, double& maxchi2_out );

	static BOOL_T FindPointsOnLine2( double* x_values, double* y_values, int cnt,
											   double max_chi2, 
											   double* line_x, double* line_y, int& cnt_on_line,
											   double& a, double& b, double& maxchi2_out );


	static BOOL_T FindPointsOnLine3( double* x_values, double* y_values, 
												int cnt,
	                                 double max_chi2,
                                    double* line_x, double* line_y, int& cnt_on_line,
                                    double& a, double& b, double& maxchi2_out,
                                    double* rejected_x, double* rejected_y, int& rejected_cnt );
                                    
	static BOOL_T FindPointsOnLine3New( sEventDesc* events, int cnt,
													double max_chi2_per_point,
													sEventDesc* events_on_line,int& cnt_on_line,
													double& a, double& b, double& c,
													double& maxchi2_out,
													sEventDesc* events_rejected, int& rejected_cnt,
													BOOL_T bCheckVelocity, double fVelocityError,
													double& minChi2PerEvent3Points, int ccd_idx,
													BOOL_T bAddFromSame=FALSE, eTrackCheckType_T type=eNormalTrack,
													double max_distance=-1,
													double MaxEstimatedPosDist=-1.00 );

	static BOOL_T TryToAddNewPoints( double* rejected_x, double* rejected_y, int rejected_cnt,
												double max_chi2_per_point,
												double* line_x, double* line_y, int& cnt_on_line,
												double& a, double& b,  double& maxchi2_out );

	static BOOL_T TryToAddNewPointsNew( sEventDesc* rejected, int rejected_cnt,
												double max_chi2_per_point,
												sEventDesc* line, int& cnt_on_line,
												double& a, double& b, double& c,
												double& maxchi2_out,
												int total_count, 
												BOOL_T (*fill_func)( void*, int, int, void* )=NULL,
												int cam_idx=0,BOOL_T bAddFromSame=FALSE,
												eTrackCheckType_T type=eNormalTrack,
												double max_distance=-1 );
												

	static int FindBest3Points( sEventDesc* events, int cnt,
	                            double& a, double& b, double& c,
	                            int* best3pos, double& minChi2_for3,
	                            BOOL_T bCheckVelocity, double fVelocityError,
	                            int ccd_idx, BOOL_T bAddFromSame=FALSE,
	                            eTrackCheckType_T type=eNormalTrack,
	                            double max_distance=-1,
	                            double MaxEstimatedPosDist=-1.00 );
	                                                          
	static BOOL_T CheckVelocity( sEventDesc* events, int cnt,
										  double fVelocityError , int maxTimeDiff=DEFAULT_MAX_DIFF_TIME,
										  int ccd_idx=0,eTrackCheckType_T type=eNormalTrack,
										  eHistoVariableType_T histo_type=eHistoVXRatioToOld,
										  double* rx_min=NULL, double* ry_min=NULL, 
										  eCheckVelocityReason_T check_velo_reason=eCheckVeloDefault,
										  double MaxEstimatedPosDist=-1.00 );
										  
	static BOOL_T CheckVelocity( sEventDesc* events, int cnt,
										  sEventDesc& newEvent, double fVelocityError,
										  int maxTimeDiff=DEFAULT_MAX_DIFF_TIME,
										  eTrackCheckType_T type=eNormalTrack,
										  eHistoVariableType_T histo_type=eHistoVXRatioToOld );	  

	static BOOL_T CheckVelocityCondition( double vx, double vy, 
													  double vx_new, double vy_new, 
													  double error,
													  int ccd_idx=DEFAULT_CCD_IDX_FOR_CHECK_COND,
													  eTrackCheckType_T type=eNormalTrack,
													  eHistoVariableType_T histo_type=eHistoVXRatioToOld );

	static BOOL_T CheckVelocityCondition( double vx, double vy, 
													  double vx_new, double vy_new, 
													  double error,
													  double& rx, double& ry,
													  int ccd_idx=DEFAULT_CCD_IDX_FOR_CHECK_COND,
													  eTrackCheckType_T type=eNormalTrack,
													  eHistoVariableType_T histo_type=eHistoVXRatioToOld,
													  eCheckVelocityReason_T check_velo_reason=eCheckVeloDefault );

	static BOOL_T CheckVelocityConditionNew( double v, double delta, double dt,
	                                         double v_new, double delta_new, double dt_new,
                                            double error,
                                            double& r,
                                            eCheckVelocityReason_T check_velo_reason,
                                            const char* szVeloDir );

	static BOOL_T CheckVelocityConditionQuick( double v, double v_err,
                 		                        double v_new, double v_new_err,
                                             double error, double& r,
                                             eCheckVelocityReason_T check_velo_reason,
                                             const char* szVeloDir );
                                                                                                                                                                    
	                                                                                                                                                                                                             													  

	static double CalcMaxDist( sEventDesc* events, int cnt );

	static double ReCalcMaxChi2( double maxDist, double maxchi2 );

	static BOOL_T DoAddFromSameFrame( eTrackCheckType_T track_type );

	static double GetGaussFast( double sigma, double mean );
	
	static void SetROOT_RandomSeed( int seed );

	// integration of gauss :
	static void* root_func;
	static double GaussIntegral( double x0, double y0, double x1, double y1 );
	static double CreateGaussFunc( double x, double y, double radius );

	// distance :
	static double calc_dist( const sEventDesc& point1 , const sEventDesc& point2 );
	static double calc_min_dist( sEventDesc* events, int cnt, const sEventDesc& point );
};


#endif

