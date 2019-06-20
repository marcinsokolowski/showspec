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
#include "myfits.h"
#include <string.h>
#include <stdio.h>
#include <cexcp.h>
#include "mathfunc.h"
// #include <datsrc.h>
// #include <grsrc.h>
#include <myutil.h>
#include <mylock.h>


void* CMyFit::root_func=NULL; 


#ifndef _NO_ROOT_
#include <TGraph.h>
#include <TF1.h>
#include <TF2.h>
#include <TRandom.h>

TF1* pGaussFunc=NULL;
CMyMutex gGlobalLock;
                                                                                
BOOL_T MyFits_FitGauss( double minValue, double& rms, double& center,
                      double &max, int binNo, double binWidth,
                      int* pCountTab, int& nstep, int first_bin )
{
	gGlobalLock.Lock();

	if( !pGaussFunc ){
		pGaussFunc = new TF1("fit_func","gaus",-100000,100000);
	}

	pGaussFunc->SetParameter( 2, rms );
	pGaussFunc->SetParameter( 1, center );
	pGaussFunc->SetParameter( 0, max );

	TGraph graph( binNo );

	double half = (binWidth/2.00);
	for(int i=0;i<binNo;i++){		
		float x_f = minValue+i*binWidth+half;
      float y_f = pCountTab[i];


		graph.SetPoint( i, x_f, y_f );
	}

	Double_t par[3];
	_TRACE_PRINTF_3("Fitting gauss....\n");

	// verbose :
	/*if( gPrintfLevel>=4 ){
		graph.Fit("gaus","V");
	}else{
		// quite :
		graph.Fit("gaus","Q");
	}*/


	int ret=0;  
	if( gPrintfLevel>=4 ){
		ret = graph.Fit( "fit_func", "VN0" ); // ms 20110330 - N added to avoid drawing 
	}else{
		ret = graph.Fit( "fit_func", "QN0" ); // ms 20110330 - N added to avoid drawing
	}

	pGaussFunc->GetParameters(par);

	/*rms = (graph.GetFunction("gaus"))->GetParameter(2);
	center = (graph.GetFunction("gaus"))->GetParameter(1);
	max = (graph.GetFunction("gaus"))->GetParameter(0);*/

	rms = par[2];
	center = par[1];
	max = par[0];
	
	_TRACE_PRINTF_2("Fit_gauss (mean,rms,norm) = ( %.4f , %.4f, %.4f ) ret=%d\n",center,rms,max,ret);

	gGlobalLock.UnLock();
	
	return (ret==0);
}

BOOL_T (*CMyFit::m_pFitGauss)( double minValue, double& rms, double& center,
	 			 double &max, int binNo, double binWidth,
	 			 int* pCountTab, int& nstep, int first_bin )=MyFits_FitGauss;

TRandom gRootRND;
CMyMutex gRootRNDLock;

double CMyFit::GetGaussFast( double sigma, double mean )
{
	gRootRNDLock.Lock();
	double ret = gRootRND.Gaus( mean, sigma );
	gRootRNDLock.UnLock();
	return ret;
}

void CMyFit::SetROOT_RandomSeed( int seed ){
	gRootRNDLock.Lock();
	printf("setting root seed = %d\n",seed);
	gRootRND.SetSeed( seed );
	gRootRNDLock.UnLock();
}

double CMyFit::GaussIntegral( double x0, double y0, double x1, double y1 )
{
	if( !root_func ){
		printf("Gauss function not initialized !!!\n");
		return -1;
	}
	TF2* func = (TF2*)root_func;
	double ret = func->Integral( x0, x1, y0, y1, 0.001 );
	return ret;
}

Double_t StarDistrGauss( Double_t* x, Double_t* y )
{
	Double_t valG = CMyMathFunc::GlobalGauss( x[0], x[1] );

	return valG;
}


double CMyFit::CreateGaussFunc( double x, double y, double radius ){
	if( root_func )
	{
		delete ((TF2*)root_func);
	}
	root_func = new TF2( "star_gauss", StarDistrGauss, x-radius, y-radius, x+radius, y+radius, 0 ); 	
}

// TF1* gParabolaFunc=NULL;
// CMyMutex gGlobalParabolaLock;

BOOL_T MyFits_FitParabola( double* x_values, double* y_values, int count,
                           double& a, double& b, double& c, double& chi2 )
{
   // TMinuit::Fitter should not be called in multiple threads !!!
	gGlobalLock.Lock();

   double min_x=1000000.00, max_x=-1000000;
   for(int i=0;i<count;i++){
     if( x_values[i] < min_x ){
       min_x = x_values[i];
     }
     if( x_values[i] > max_x ){
       max_x = x_values[i];
     }
   }

   min_x = min_x - 10.00;
   max_x = max_x + 10.00;

   if( gPrintfLevel >= 4 ){   
     printf("MyFits_FitParabola : fiting in range %.2f - %.2f\n",min_x,max_x);
   }

	TF1 gParabolaFunc("fit_parabola_func","[2]*x*x+[1]*x+[0]",min_x,max_x);

//	if( !gParabolaFunc ){
//		gParabolaFunc = new TF1("fit_func","gaus",-100000,100000);
//	}

	gParabolaFunc.SetParameter( 2, a );
	gParabolaFunc.SetParameter( 1, b );
	gParabolaFunc.SetParameter( 0, c );

	TGraph graph( count );

	for(int i=0;i<count;i++){		
		graph.SetPoint( i, x_values[i], y_values[i] );
	}

	Double_t par[3];
	_TRACE_PRINTF_3("MyFits_FitParabola : Fitting parabola....\n");

	int ret=0;
	if( gPrintfLevel>=4 ){
		ret = graph.Fit( "fit_parabola_func", "V" );
	}else{
		ret = graph.Fit( "fit_parabola_func", "Q" );
	}

	gParabolaFunc.GetParameters(par);

	a = par[2];
	b = par[1];
	c = par[0];
	
	chi2 = gParabolaFunc.GetChisquare();
	
	_TRACE_PRINTF_2("MyFits_FitParabola : fitted parabola %e*x^2+%e*x+%e=0\n",a,b,c);

	gGlobalLock.UnLock();
	
	return (ret==0);
}


#else
BOOL_T (*CMyFit::m_pFitGauss)( double minValue, double& rms, double& center,
	 			 double &max, int binNo, double binWidth,
	 			 int* pCountTab, int& nstep, int first_bin )=NULL;

BOOL_T MyFits_FitParabola( double* x_values, double* y_values, int count,
                           double& a, double& b, double& c, double& chi2 )
{
  printf("ERROR : without root parabola fitting in MyFits_FitParabola - NOT IMPLEMENTED !!!\n");
  return FALSE;
}                           

double CMyFit::CreateGaussFunc( double x, double y, double radius ){
	return 0;
}

double CMyFit::GetGaussFast( double sigma, double mean )
{
	return 0;
}

double CMyFit::GaussIntegral( double x0, double y0, double x1, double y1 )
{
	return CMyMathFunc::GlobalGauss( (x0+x1)/2.00, (y0+y1)/2.00 );
}

#endif

BOOL_T (*CMyFit::m_pFillFunc)( void* event_info, int type, int ccd_idx, void* event_info2 )=NULL;



// [MS] - change on 20041003 - before 0.1 , but tracks with 0.06 
// were not checked and velocity was not checked correctly :
// now set to 0.05
// 20041005 - set to 0.02
double gMinVeloToCheck=0.02;

// in case relative velocity error (dv/v)=(dx/x) > gMaxVeloRelativeError
// no velocity check is performed !
double gMaxVeloRelativeError=1.00;

// This is paramter allowing to specify maximum distance between frames 
// that is allowed for points in the 3 point seed :
// corresponds to paramter : CCD_TRACK_MAX_FRAME_DIST_FOR_3POINTS_SEED
int gMaxFrameDistFor3PointSeed=-1;

static void get_xy( sEventDesc* list, int cnt, double* x, double* y )
{
	for( register int i=0;i<cnt;i++){
		x[i] = list[i].x;
		y[i] = list[i].y;
	}
}

void sort_event_list( sEventDesc* list, int cnt )
{
	for( register int i=0;i<cnt-1;i++ ){
		register int min_value=list[i].timeUT;
		register int min_pos=i;
		for( register int j=i+1;j<cnt;j++){
			if( list[j].timeUT < min_value ){
				min_value = list[j].timeUT;
				min_pos = j;
			}
		}
		if( i!=min_pos ){
			sEventDesc tmp;
			tmp = list[i];
			list[i] = list[min_pos];
			list[min_pos] = tmp;
		}
	}	
}

static int find_frame( sEventDesc* list, int cnt, int frame )
{
	for( register int i=0;i<cnt;i++){
		if( list[i].frame == frame )
			return i;
	}
	return -1;
}




CMyFit::CMyFit()
{

}

CMyFit::~CMyFit()
{

}

double CMyFit::FitLineChi2( double* x_values, double* y_values, int cnt,
                            double& a, double& b, double& c )
{
	FitLine( x_values, y_values, cnt, a, b, c  );
	double chi2 = CalcChi2_Dist( x_values, y_values, cnt, a, b, c   );
	return chi2;
}

BOOL_T CMyFit::FitLineHorizontal( double* x_values, double* y_values, int cnt,
                                  double& c )
{
	return FALSE;
}

/* 2010-03-13 : changed interpretation of a,b from y = ax + b to 
                canonical ax + by + c = 0 to be able to handle vertical
                and horizontal lines !
*/            
BOOL_T CMyFit::FitLine( double* x_values, double* y_values, int cnt,
                        double& a, double& b, double& c,
                        int exceptPos/*=-1*/ )
{
	register double xy_sum = 0;
	register double x_sum  = 0;
	register double x2_sum = 0;
	register double y_sum  = 0; 	

	int usedCount=0;
	for(register int i=0;i<cnt;i++){
		if(exceptPos<0 || i!=exceptPos){
			xy_sum += x_values[i]*y_values[i];
			x_sum  += x_values[i];
			y_sum  += y_values[i];
			x2_sum += (x_values[i]*x_values[i]);
			usedCount++;
		}
	}

	
	double bottom = usedCount*x2_sum - x_sum*x_sum;
	
	if( bottom == 0 ){
	   a = 1;
	   b = 0;
	   c = - x_sum/usedCount;
	}else{
    	a = (usedCount*xy_sum-x_sum*y_sum)/bottom;
      b = -1;
    	c = (x2_sum*y_sum-x_sum*xy_sum)/bottom;
   }

	return TRUE;
}


double CMyFit::getLineChi2( double x, double y, double a, double b, double c )
{
   double line_y = 0.00;
   double chi2 = 0.00;
   if( b != 0 ){
	  line_y = -(a*x+c)/b;
     chi2 = (y-line_y)*(y-line_y);	   
   }else{
     double line_x = (-c/a);
     chi2 = (x-line_x)*(x-line_x);
   }
	
	return chi2;
}


double CMyFit::CalcChi2( double* x_values, double* y_values, int cnt,
                         double a, double b, int exceptPos/*=-1*/ )
{
	double sum = 0;
	for(register int i=0;i<cnt;i++){
		if(exceptPos<0 || exceptPos!=i){
			sum += getLineChi2( x_values[i], y_values[i], a, -1, b );
		}
	}
	return sum;
}

double CMyFit::CalcChi2_Dist( double* x_values, double* y_values, int cnt,
                         double a, double b, double c,
                         int exceptPos/*=-1*/ )
{
	double sum = 0;
	for(register int i=0;i<cnt;i++){
		if(exceptPos<0 || exceptPos!=i){
			sum += CalcDist2FromLine( a, b , c , x_values[i], y_values[i] );
		}
	}
	return sum;
}


double CMyFit::CalcMaxChi2( double* x_values, double* y_values, int cnt,
									 double a, double b, int& pos )
{
	double max_chi2 = -1000;
	
	pos = -1;
	for(register int i=0;i<cnt;i++){
		double chi2 = getLineChi2( x_values[i], y_values[i], a, -1, b );
		if(chi2 > max_chi2){
			max_chi2 = chi2;
			pos = i;
		}
	}

	return max_chi2;
}

void CMyFit::RejectPoint( double* x_values, double* y_values, int cnt, int pos )
{
	int i=pos+1;

	printf("rejected point : (%.5f,%.5f)\n",x_values[pos],y_values[pos]);
	while(i<cnt){
		x_values[i-1] = x_values[i];
		y_values[i-1] = y_values[i];
		i++;
	}
}

BOOL_T CMyFit::FindPointsOnLine( double* x_values, double* y_values, int cnt,
										   double max_chi2, 
										   double* line_x, double* line_y, int& cnt_on_line,
										 	double& a, double& b, double& maxchi2_out )
{
	double aa,bb,cc;
	int pos;


	double* wrk_x = new double[cnt];
	double* wrk_y = new double[cnt];

	memcpy( wrk_x, x_values, cnt*sizeof(double) );
	memcpy( wrk_y, y_values, cnt*sizeof(double) );
	
	
	FitLine( wrk_x, wrk_y, cnt, aa, bb, cc );
	double maxchi2 = CalcMaxChi2( wrk_x, wrk_y, cnt, aa, cc, pos );

	cnt_on_line = cnt;
	while(maxchi2>max_chi2 && pos>=0 && cnt_on_line>2){
		RejectPoint( wrk_x, wrk_y, cnt_on_line, pos );
		cnt_on_line--;											
		maxchi2 = CalcMaxChi2( wrk_x, wrk_y, cnt_on_line, aa, bb, pos );
	}

	memcpy( line_x, wrk_x, cnt_on_line*sizeof(double) );
	memcpy( line_y, wrk_y, cnt_on_line*sizeof(double) );
	
	delete [] wrk_x;
	delete [] wrk_y;

	maxchi2_out = maxchi2;
	return (maxchi2<=max_chi2 && cnt_on_line>2);
}


BOOL_T CMyFit::FindPointsOnLine2( double* x_values, double* y_values, int cnt,
										   double max_chi2, 
										   double* line_x, double* line_y, int& cnt_on_line,
										 	double& a, double& b, double& maxchi2_out )
{
	double aa,bb,cc;
	int pos;


	double* wrk_x = new double[cnt];
	double* wrk_y = new double[cnt];
	double* chi2  = new double[cnt];

	memcpy( wrk_x, x_values, cnt*sizeof(double) );
	memcpy( wrk_y, y_values, cnt*sizeof(double) );
	
	
	FitLine( wrk_x, wrk_y, cnt, aa, bb,cc );
	double maxchi2 = CalcMaxChi2( wrk_x, wrk_y, cnt, aa, cc, pos );
	double allChi2 = CalcChi2_Dist( wrk_x, wrk_y, cnt, aa, bb, cc );

	cnt_on_line = cnt;
	while(maxchi2>max_chi2 && pos>=0 && cnt_on_line>2){
		

		double dx1 = (x_values[1]-x_values[0]);
		double dy1 = (y_values[1]-y_values[0]);
		if(maxchi2<=500.00){
			// reject not chrionological :
			
			int bCont=1;
			while(bCont){
				bCont = 0;
				for(int k=2;k<cnt_on_line;k++){
					double dx = (wrk_x[k]-wrk_x[k-1]);		
					double dy = (wrk_y[k]-wrk_y[k-1]);

					if(dx*dx1<0 || dy*dy1<0){
						bCont = 1;						
						RejectPoint( wrk_x, wrk_y, cnt_on_line, k );
						cnt_on_line--;
						break;
					}
				}			
			}
		}

		double minChi2 = 10000000.00;		
		int minPos=-1;
		double bestA,bestB;
		for(register int i=0;i<cnt_on_line;i++){
			FitLine( wrk_x, wrk_y, cnt_on_line, aa, bb, cc, i );
			chi2[i] = CalcChi2_Dist( wrk_x, wrk_y, cnt_on_line, aa, bb, cc, i );	
			if(chi2[i]<minChi2){
				minChi2 = chi2[i];
				minPos = i;
				bestA = aa;
				bestB = bb;
			}
		}				


		// it was best to reject point at minPos :
		if(minPos>=0){
			RejectPoint( wrk_x, wrk_y, cnt_on_line, minPos );
			cnt_on_line--;			
			maxchi2 = CalcMaxChi2( wrk_x, wrk_y, cnt_on_line, bestA, bestB, pos );
			allChi2 = minChi2;
		}else{
			Assert(FALSE,"Could not fit line !!!");
		}
	}

	memcpy( line_x, wrk_x, cnt_on_line*sizeof(double) );
	memcpy( line_y, wrk_y, cnt_on_line*sizeof(double) );
	
	delete [] wrk_x;
	delete [] wrk_y;
	delete [] chi2;

	maxchi2_out = maxchi2;
	return (maxchi2<=max_chi2 && cnt_on_line>2);
}

double CMyFit::CalcMaxDist( sEventDesc* events, int cnt )
{
	double maxDist=-1;
	for(register int i=0;i<cnt;i++){
		for(register int j=(i+1);j<cnt;j++){
			double dist2 = ( (events[j].x-events[i].x)*(events[j].x-events[i].x)+
								  (events[j].y-events[i].y)*(events[j].y-events[i].y) );
			double dist = sqrt(dist2);
			if( dist>maxDist ){
				maxDist = dist;
			}
		}
	}
	return maxDist;
}

BOOL_T CMyFit::CheckVelocity( sEventDesc* events, int cnt,
                              double fVelocityError, int maxTimeDiff,
										int ccd_idx, eTrackCheckType_T type,
										eHistoVariableType_T histo_type,
										double* rx_min, double* ry_min,
										eCheckVelocityReason_T check_velo_reason,
										double MaxEstimatedPosDist /*=-1.00*/ )
{
	if(cnt<=2)
		return TRUE;

	// maybe to be sure add sort by timeUT here 
	// this function assumes table events is already sorted in this
	// maner ...
	int dt0 = (events[1].timeUT-events[0].timeUT);

	if( dt0!=0 ){
		double vx0 = (events[1].x-events[0].x)/dt0;
		double vy0 = (events[1].y-events[0].y)/dt0;
		double rx=0,ry=0;
			
		if( rx_min && ry_min ){
			(*rx_min) = 100000.000;
			(*ry_min) = 100000.000;
		}
	

		if( abs(dt0)>maxTimeDiff )
			return FALSE;

		for(register int i=1;i<(cnt-1);i++){	
			int dt = (events[i+1].timeUT-events[i].timeUT);

			if( dt!=0 ){
				// check only if events from different frame :
				double vx = (events[i+1].x-events[i].x)/dt;
				double vy = (events[i+1].y-events[i].y)/dt;			

				if( abs(dt)>maxTimeDiff )
	      		return FALSE;
		
				BOOL_T bCheck= CheckVelocityCondition( vx0, vy0, vx, vy, fVelocityError, rx , ry, ccd_idx, type, histo_type, check_velo_reason );
		
				if( rx_min && ry_min ){
					if( rx < (*rx_min) ){
						(*rx_min) = rx;
					}
					if( ry < (*ry_min) ){
						(*ry_min) = ry;
					}
				}

				if( !bCheck ){
					return FALSE;
				}
			}
		}	


		if( MaxEstimatedPosDist>0 && cnt>=3 ){
		  // calculate estimated position of event[2] starting from event[1]
		  double est_x = events[1].x + vx0*(events[2].timeUT-events[1].timeUT);
    	  double est_y = events[1].y + vy0*(events[2].timeUT-events[1].timeUT);

	     double dist_from_estimated = sqrt( (est_x-events[2].x)*(est_x-events[2].x) + (est_y-events[2].y)*(est_y-events[2].y));
    	  if( gDebugTracks ){
	        printf("DEBUG_TRACKS : CheckVelocity : estimated position (%.2f,%.2f) vs observed (%.2f,%.2f) , distance = %.2f compared to limit %.2f\n",
	            est_x,est_y,events[2].x,events[2].y,dist_from_estimated,MaxEstimatedPosDist);
    	  }
	     if( dist_from_estimated > MaxEstimatedPosDist ){	    
	       if( gDebugTracks ){
            printf("DEBUG_TRACKS : CheckVelocity returns FALSE due to estimated position disagreement\n");
          }
          return FALSE;
        }
    	}
	}
	
	return TRUE;
}

BOOL_T CMyFit::CheckVelocity( sEventDesc* events, int cnt,
                              sEventDesc& newEvent, 
										double fVelocityError, 
										int maxTimeDiff, eTrackCheckType_T type,
										eHistoVariableType_T histo_type )
{
	if( cnt<2 )
		return 1;
	time_t startTime=events[0].timeUT;
	time_t endTime=events[cnt-1].timeUT;
	sEventDesc* pStart=(&(events[0]));
	sEventDesc* pEnd=(&(events[cnt-1]));
	sEventDesc* pBefore=NULL;
	sEventDesc* pAfter=NULL;
	time_t beforeDT=1000000;
	time_t  afterDT=1000000;

	for(register int i=0;i<cnt;i++){
		if(events[i].timeUT<startTime){
			startTime = events[i].timeUT;
			pStart = &(events[i]);
		}
		if(events[i].timeUT>endTime){
			endTime = events[i].timeUT;
			pEnd = &(events[i]);
		}
		int dt = (newEvent.timeUT - events[i].timeUT);
		if( dt>0 && dt<beforeDT ){
			beforeDT = dt;
			pBefore = &(events[i]);
		}else{
			dt = -dt;
			if( dt>0 && dt<afterDT ){
				afterDT = dt;
				pAfter = &(events[i]);
			}
		}
	}
	if( pStart && pEnd && (pAfter || pBefore)){
		time_t dt = (pEnd->timeUT - pStart->timeUT);
		double vx = (pEnd->x - pStart->x)/dt;
		double vy = (pEnd->y - pStart->y)/dt;

		
		double vx_new,vy_new;
		if( pAfter ){
			if( abs(afterDT)>maxTimeDiff  ){
				return FALSE;
			}

			vx_new = (pAfter->x-newEvent.x)/afterDT;
			vy_new = (pAfter->y-newEvent.y)/afterDT;
			/*if( !CheckVelocityCondition( vx, vy, vx_new, vy_new, fVelocityError, 
												  DEFAULT_CCD_IDX_FOR_CHECK_COND, type, histo_type) ){
				return FALSE;
			}*/

			// NEW - 20041118 - velocity checked only if >1 pixels difference
			if( fabs(pAfter->y-newEvent.y)>1 && 
				 !CheckVelocityCondition( 1.00, vy, 1.00, vy_new, fVelocityError, 
												  DEFAULT_CCD_IDX_FOR_CHECK_COND, type, histo_type) ){
				return FALSE;
			}
			if( fabs(pAfter->x-newEvent.x)>1 && 
				 !CheckVelocityCondition( vx, 1.00, vx_new, 1.00, fVelocityError, 
												  DEFAULT_CCD_IDX_FOR_CHECK_COND, type, histo_type) ){
				return FALSE;
			}

		}
		if( pBefore ){
			vx_new = (newEvent.x-pBefore->x)/beforeDT;
			vy_new = (newEvent.y-pBefore->y)/beforeDT;
			if( abs(beforeDT)>maxTimeDiff  ){
				return FALSE;
			}

			/*if( !CheckVelocityCondition( vx, vy, vx_new, vy_new, fVelocityError, 
												  DEFAULT_CCD_IDX_FOR_CHECK_COND, type, histo_type ) ){ 
				return FALSE;
			}*/

			// NEW - 20041118 - velocity checked only if >1 pixels difference
			if( fabs(pBefore->y-newEvent.y)>1 && 
				 !CheckVelocityCondition( 1.00, vy, 1.00, vy_new, fVelocityError, 
												  DEFAULT_CCD_IDX_FOR_CHECK_COND, type, histo_type) ){
				return FALSE;
			}
			if( fabs(pBefore->x-newEvent.x)>1 && 
				 !CheckVelocityCondition( vx, 1.00, vx_new, 1.00, fVelocityError, 
												  DEFAULT_CCD_IDX_FOR_CHECK_COND, type, histo_type) ){
				return FALSE;
			}
		}
	}
	return TRUE;
}

BOOL_T CMyFit::CheckVelocityCondition( double vx, double vy,
                                       double vx_new, double vy_new,
                                       double error, int ccd_idx,
													eTrackCheckType_T type,
													eHistoVariableType_T histo_type )
{
/*	if( vx_new*vx<0 || vy_new*vy<0 || fabs( vx_new-vx )>error || fabs( vy_new-vy )>error ){
		return FALSE;
	}*/
	double vx_min=MIN(fabs(vx),fabs(vx_new));
	double vy_min=MIN(fabs(vy),fabs(vy_new));
	double vx_max=MAX(fabs(vx),fabs(vx_new));
	double vy_max=MAX(fabs(vy),fabs(vy_new));	
   double rx = ( vx_min/vx_max )*mysign(vx*vx_new);
   if(vx_max<gMinVeloToCheck)
     rx = 1.00;


   double ry = ( vy_min/vy_max )*mysign(vy*vy_new);
   if(vy_max<gMinVeloToCheck)
     ry = 1.00;

	if( m_pFillFunc ){
		(*m_pFillFunc)( &rx, eRXall, ccd_idx, NULL );
		(*m_pFillFunc)( &ry, eRYall, ccd_idx, NULL );
		(*m_pFillFunc)( &vx, eVXvsVX, ccd_idx, &vx_new );
		(*m_pFillFunc)( &vy, eVYvsVY, ccd_idx, &vy_new );
	}

	if(gDebugTracks>1){
	  printf("DEBUG_TRACK : CMyFit::CheckVelocityCondition  : (vx_min,vx_max)=( %.2f , %.2f ) , (vy_min,vy_max)=( %.2f , %.2f )\n",vx_min,vx_max,vy_min,vy_max);
     printf("DEBUG_TRACK : CMyFit::CheckVelocityCondition  : (rx,ry)=( %.2f , %.2f )\n",rx,ry);
           
     double rx2 = (vx_max-vx_min);
     double ry2 = (vy_max-vy_min);
     if( (vx_max+vx_min) != 0 ){
       rx2 = rx2/(vx_max+vx_min);
       printf("DEBUG_TRACK VELO_TEST : rx2 = %.2f\n",rx2);
     }
     if( (vy_max+vy_min) != 0 ){
       ry2 = ry2/(vy_max+vy_min);
       printf("DEBUG_TRACK VELO_TEST : ry2 = %.2f\n",ry2);
     }
   }

	// if max is slower then gMinVeloToCheck pixel / sec - do not perform this check :
	if( ( vx_max>=gMinVeloToCheck && ( vx_new*vx<0 || (vx_min/vx_max)<error ) ) || // bad vx 	
		 ( vy_max>=gMinVeloToCheck && ( vy_new*vy<0 || (vy_min/vy_max)<error ) ) ){ // or bad vy
		if(gDebugTracks>1)printf(" return FALSE\n");
      return FALSE;
	} 


	/*if( (vx_new*vx<0 && vx_max>=gMinVeloToCheck ) || 
		 (vy_new*vy<0 && vy_max>=gMinVeloToCheck ) || 
		 ( (vx_min/vx_max)<error && vx_max>=gMinVeloToCheck ) ||  // special for very slow 
		 ( (vy_min/vy_max)<error && vy_max>=gMinVeloToCheck ) ){  // special for very slow
		if(gDebugTracks>1)printf("return FALSE\n");
		return FALSE;
	}*/

	if(gDebugTracks>1)printf(" return TRUE\n");
	
	return TRUE;	
}


BOOL_T CMyFit::CheckVelocityCondition( double vx, double vy, 
													  double vx_new, double vy_new, 
													  double error,
													  double& rx, double& ry,
													  int ccd_idx, 
													  eTrackCheckType_T type,
													  eHistoVariableType_T histo_type,
													  eCheckVelocityReason_T check_velo_reason )
{
	double vx_min=MIN(fabs(vx),fabs(vx_new));
	double vy_min=MIN(fabs(vy),fabs(vy_new));
	double vx_max=MAX(fabs(vx),fabs(vx_new));
	double vy_max=MAX(fabs(vy),fabs(vy_new));	

	

	if( m_pFillFunc ){
		double rx = ( vx_min/vx_max )*mysign(vx*vx_new);
		if(vx_max<gMinVeloToCheck)
			rx = 1.00;


		double ry = ( vy_min/vy_max )*mysign(vy*vy_new);
		if(vy_max<gMinVeloToCheck)
			ry = 1.00;

		(*m_pFillFunc)( &rx, eRXall, ccd_idx, NULL );
		(*m_pFillFunc)( &ry, eRYall, ccd_idx, NULL );
		(*m_pFillFunc)( &vx, eVXvsVX, ccd_idx, &vx_new );
		(*m_pFillFunc)( &vy, eVYvsVY, ccd_idx, &vy_new );

		if(gDebugTracks>1)printf("DEBUG_TRACK : (rx,ry)=(%.2f,%.2f)\n",rx,ry);
	}

	// if max is slower then gMinVeloToCheck pixel / sec - do not perform this check :
	double rx2=-10.00,ry2=-10.00;
	if( vx_max!=0 ){
		rx = ( vx_min/vx_max )*mysign(vx*vx_new); // mysign NEW [20041006]
	}else{
		rx = 1.00;
	}
	if( vy_max!=0 ){
	   ry = ( vy_min/vy_max )*mysign(vy*vy_new); // mysign NEW [20041006]
	}else{
		ry = 1.00;
	}

	if( (vx_max-vx_min)!=0 ){
      rx2 = (vx_max-vx_min)/(vx_max+vx_min);	  
	}
	if( (vy_max-vy_min)!=0 ){
      ry2 = (vy_max-vy_min)/(vy_max+vy_min);	  
	}

	if(gDebugTracks>1){
	  printf("DEBUG_TRACK : CMyFit::CheckVelocityCondition_%d (vx_min,vx_max) = ( %.2f , %.2f ), (vy_min,vy_max)=( %.2f , %.2f ) , (rx,ry)=( %.6f , %.6f ) , (rx2,ry2)=( %.6f , %.6f ) , (vx,vx_new,vy,vy_new)=( %.6f , %.6f , %.6f , %.6f )=>",check_velo_reason,vx_min,vx_max,vy_min,vy_max,rx,ry,rx2,ry2,vx,vx_new,vy,vy_new);
   }

	if( ( vx_max>=gMinVeloToCheck && ( vx_new*vx<0 || rx<error ) ) || // bad vx 	
		 ( vy_max>=gMinVeloToCheck && ( vy_new*vy<0 || ry<error ) ) ){ // or bad vy
		if(gDebugTracks>1)printf(" return FALSE\n");
      return FALSE;
	} 

	if( gDebugTracks>1 ){
		printf("DEBUG_TRACK : (vx_min,vx_max) = (%.2f,%.2f), (vy_min,vy_max) = (%.2f,%.2f) ==> (rx,ry) = (%.2f,%.2f)\n",vx_min,vx_max,vy_min,vy_max,rx,ry);
	}

	/*if( (vx_new*vx<0 && vx_max>=gMinVeloToCheck ) || 
		 (vy_new*vy<0 && vy_max>=gMinVeloToCheck ) || 
		 ( (vx_min/vx_max)<error && vx_max>=gMinVeloToCheck ) ||  // special for very slow 
		 ( (vy_min/vy_max)<error && vy_max>=gMinVeloToCheck ) ){  // special for very slow
		if(gDebugTracks>1)printf("return FALSE\n");
		return FALSE;
	}*/
	if(gDebugTracks>1)printf(" return TRUE\n");
	
	return TRUE;	

}

BOOL_T CMyFit::CheckVelocityConditionNew( double v, double delta, double dt, 
													  double v_new, double delta_new, double dt_new,
													  double error,
													  double& r,
													  eCheckVelocityReason_T check_velo_reason,
													  const char* szVeloDir )
{
   double position_error=2.00; // 2 pixels - position error for flash / sat 
	double v_min=MIN(fabs(v),fabs(v_new));
	double v_max=MAX(fabs(v),fabs(v_new));

   delta = fabs(delta);
   delta_new = fabs(delta_new);	

	// if max is slower then gMinVeloToCheck pixel / sec - do not perform this check :
	double r2=-10.00;
	if( v_max!=0 ){
		r = ( v_min/v_max )*mysign(v*v_new); // mysign NEW [20041006]
	}else{
		r = 1.00;
	}

	if( (v_max-v_min)!=0 ){
      r2 = (v_max-v_min)/(v_max+v_min);	  
	}

	if(gDebugTracks>1){
	  printf("DEBUG_TRACK : CMyFit::CheckVelocityCondition_%d %s (v_min,v_max) = ( %.2f , %.2f ), r = %.6f , r2 = %.6f , (v,v_new)=( %.6f , %.6f )=>",check_velo_reason,szVeloDir,v_min,v_max,r,r2,v,v_new);
   }

   BOOL_T bDoCheckVelo=FALSE;
   // gMaxVeloRelativeError 

   double v_err = position_error/dt;
   double v_new_err = position_error/dt_new;
   if( fabs(v) > v_err || fabs(v_new) > v_new_err ){
      // at least one velocity is large enough for check 
      // otherwise velocity check is not performed 
      bDoCheckVelo = TRUE;
   }

   double v_diff_err = sqrt(v_err*v_err+v_new_err*v_new_err);
   if( fabs( v-v_new) < 3.00*v_diff_err ){
     // do not check if difference between velocities is smaller 
     // than 2sigma error , in fact 1 sigma would be enough
     // if not star+sat problem, when sat goes over a star
     // cluster includes the star and position of satellite equals 
     // star position causing large error !!!
     bDoCheckVelo = FALSE;
   }   

	if( bDoCheckVelo ){
	   if( v_new*v<0 || r<error ){
    		if(gDebugTracks>1)printf(" return FALSE\n");
         return FALSE;
      }
      
      // for slow tracks ( <1.00 px/sec === ~ 12px/frame - taka sobie kreska )
//      if( v_max<1.00 && fabs( v-v_new ) > 10.00*v_diff_err ){
        // if difference is much larger than error -> reject
//        if(gDebugTracks>1)printf(" return FALSE ( due to too large difference %.4f > 10.00*%.4f )\n",fabs( v-v_new),v_diff_err);
//        return FALSE;
//      }
	} 

	if( gDebugTracks>1 ){
		printf("DEBUG_TRACK : %s : (v_min,v_max) = (%.2f,%.2f) ==> rx = %.2f\n",szVeloDir,v_min,v_max,r);
	}

	if(gDebugTracks>1)printf(" return TRUE\n");
	
	return TRUE;	

}

BOOL_T CMyFit::CheckVelocityConditionQuick( double v, double v_err,
													  double v_new, double v_new_err,
													  double error, double& r,
													  eCheckVelocityReason_T check_velo_reason,
													  const char* szVeloDir )
{
	double v_min=MIN(fabs(v),fabs(v_new));
	double v_max=MAX(fabs(v),fabs(v_new));


	// if max is slower then gMinVeloToCheck pixel / sec - do not perform this check :
	double r2=-10.00;
	if( v_max!=0 ){
		r = ( v_min/v_max )*mysign(v*v_new); // mysign NEW [20041006]
	}else{
		r = 1.00;
	}

	if( (v_max-v_min)!=0 ){
      r2 = (v_max-v_min)/(v_max+v_min);	  
	}

	if(gDebugTracks>1){
	  printf("DEBUG_TRACK : CMyFit::CheckVelocityConditionQuick_%d %s (v_min,v_max) = ( %.2f , %.2f ), r = %.6f , r2 = %.6f , (v,v_new)=( %.6f , %.6f )=>",check_velo_reason,szVeloDir,v_min,v_max,r,r2,v,v_new);
   }

   BOOL_T bDoCheckVelo=FALSE;
   // gMaxVeloRelativeError 

   if( fabs(v) > v_err || fabs(v_new) > v_new_err ){
      // at least one velocity is large enough for check 
      // otherwise velocity check is not performed 
      bDoCheckVelo = TRUE;
   }

   double v_diff_err = sqrt(v_err*v_err+v_new_err*v_new_err);
   if( fabs( v-v_new) < 3.00*v_diff_err ){
     // do not check if difference between velocities is smaller 
     // than 2sigma error , in fact 1 sigma would be enough
     // if not star+sat problem, when sat goes over a star
     // cluster includes the star and position of satellite equals 
     // star position causing large error !!!
     bDoCheckVelo = FALSE;
   }   

	if( bDoCheckVelo ){
	   if( v_new*v<0 || r<error ){
    		if(gDebugTracks>1)printf(" return FALSE\n");
         return FALSE;
      }

      // for slow tracks ( <1.00 px/sec === ~ 12px/frame - taka sobie kreska )        
//      if( v_max<1.00 && fabs( v-v_new ) > 10.00*v_diff_err ){
        // if difference is much larger than error -> reject
//        if(gDebugTracks>1)printf(" return FALSE ( due to too large difference %.4f > 10.00*%.4f )\n",fabs( v-v_new),v_diff_err);
//        return FALSE;
//      }
	} 

	if( gDebugTracks>1 ){
		printf("DEBUG_TRACK : %s : (v_min,v_max) = (%.2f,%.2f) ==> rx = %.2f\n",szVeloDir,v_min,v_max,r);
	}

	if(gDebugTracks>1)printf(" return TRUE\n");
	
	return TRUE;	

}


//	FindBest3Points( wrk_x, wrk_y, frames, cnt, a, b, best3pos, minChi2_for3 );
int CMyFit::FindBest3Points( sEventDesc* events, int cnt,
									  double& a, double& b, double& c,
									  int* best3pos, double& minChi2_for3,
									  BOOL_T bCheckVelocity, double fVelocityError,
									  int ccd_idx, BOOL_T bAddFromSame,
									  eTrackCheckType_T type, 
									  double max_distance /*=-1*/,
									  double MaxEstimatedPosDist /*=-1.00*/ )
{
   BOOL_T bForceSameFrame = (type==eSingleCamTrack);   

	double pointsX[3],pointsY[3];
	int nTry=0;
	int nOK=0;
	double rx_min, ry_min;

	for(register int i=0;i<=cnt-3;i++){
		pointsX[0]=events[i].x;
		pointsY[0]=events[i].y;

		for(register int j=(i+1);j<=(cnt-2);j++){
			if( events[j].frame==events[i].frame ){
			   if( !bAddFromSame ){
    				// same frame 
	    			continue;
	    		}
			}else{
			  if( bForceSameFrame ){
			    continue;
			  }
			}

			if( max_distance > 0 ){
			   // if required skip points that are too far :
            if( calc_dist( events[j], events[i] ) > max_distance ){              
              continue;
            }			  
			}

			pointsX[1]=events[j].x;
	      pointsY[1]=events[j].y;
			for(register int k=j+1;k<=(cnt-1);k++){
				if( (events[j].frame==events[k].frame ||
					  events[j].frame==events[i].frame ) && !bAddFromSame ){
					// same frame
					continue;
				}
				
				if( bForceSameFrame ){
				  if( events[j].frame!=events[k].frame || events[j].frame!=events[i].frame ){
                continue;
				  }
				}

    			if( max_distance > 0 ){
	    		   // if required skip points that are too far :
              if( calc_dist( events[j], events[k] ) > max_distance || 
                  calc_dist( events[i], events[k] ) > max_distance ){    
                continue;
              }			  
    			}

				// have 3 points i,j,k now check them :
				double a3,b3,c3;
				pointsX[2]=events[k].x;
		      pointsY[2]=events[k].y;

            double chi2_to_show = 0.00;
            if( gDebugTracks ){
               double a_test,b_test,c_test;
    				FitLine( pointsX, pointsY, 3, a_test, b_test, c_test );
  	    			double chi2_3 = CalcChi2( pointsX, pointsY, 3, a_test, c_test );
  	    			double chi2_3_dist = CalcChi2_Dist( pointsX, pointsY, 3, a_test, b_test, c_test );
               printf("DEBUG_TRACK : FIND_3_POINTS CMyFit::FindBest3Points %.4f %.4f\n",(chi2_3/3),(chi2_3_dist/3));
               chi2_to_show = (chi2_3/3.00);
            }


				if( bCheckVelocity ){
					sEventDesc tocheck[3];
					tocheck[0]=events[i];
					tocheck[1]=events[j];
					tocheck[2]=events[k];
               BOOL_T bCheck = CheckVelocity( tocheck, 3,fVelocityError, 
                                             DEFAULT_MAX_DIFF_TIME, ccd_idx, 
                                             type, eHistoVXRatioTo3, &rx_min, 
                                             &ry_min, eCheckVeloDefault,
                                             MaxEstimatedPosDist );
	
					if( m_pFillFunc ){ 
						(*m_pFillFunc)( &rx_min, eHistoVXRatioTo3, ccd_idx, NULL );
		            (*m_pFillFunc)( &ry_min, eHistoVYRatioTo3, ccd_idx, NULL );	
					}
					
					if( !bCheck ){
						// 3 point do not satisfy velocity criteria 
						continue;
					}
            }
            
            if( MaxEstimatedPosDist > 0 ){
              
            }

            if( gMaxFrameDistFor3PointSeed > 0 ){
               // assumes that points are sorted by time ( which is true )
               // it is done before call of this function in ccd_analyse.cpp 
               // look for : sort_event_list               
               int max_frame_dist = (events[j].frame - events[i].frame);
               if( (events[k].frame-events[j].frame) > max_frame_dist ){
                 max_frame_dist = (events[k].frame-events[j].frame);
               }
                                             
               if( max_frame_dist > gMaxFrameDistFor3PointSeed ){                   
                   if( gDebugTracks ){
                     printf("DEBUG_TRACK : MaxFrameDistCHECK %d > %d -> 3 POINTS SKIPPED ( chi2 = %.4f )\n",max_frame_dist,gMaxFrameDistFor3PointSeed,chi2_to_show);
                   }
                   continue;
               }else{
                   if( gDebugTracks ){
                     printf("DEBUG_TRACK : MaxFrameDistCHECK %d <= %d -> 3 POINTS ACCEPTED FOR SEED ( chi2 = %.4f )\n",max_frame_dist,gMaxFrameDistFor3PointSeed,chi2_to_show);
                   }
               }
            }


				FitLine( pointsX, pointsY, 3, a3, b3, c3 );
				double chi2_3 = CalcChi2_Dist( pointsX, pointsY, 3, a3, b3, c3 );

				if( gShowChi2Points3 && chi2_3<1000.00 ){
					printf("chi2_3points = %.4f\n",chi2_3);fflush(0);
				}

				if(chi2_3<minChi2_for3 || nTry==0){					
					minChi2_for3 = chi2_3;
					best3pos[0]=i;
					best3pos[1]=j;
					best3pos[2]=k;		
					a = a3;
					b = b3;
					c = c3;
					nOK++;

					// just for testing reasons :
					/*sEventDesc tocheck[3];
               tocheck[0]=events[i];
               tocheck[1]=events[j];
               tocheck[2]=events[k];
               BOOL_T bVelo = CheckVelocity( tocheck, 3,fVelocityError );*/
				}
				nTry++;
			}
		}
	}
	return nOK;
}

double CMyFit::ReCalcMaxChi2( double maxDist, double maxchi2 )
{
	// return maxchi2;

	int r = (maxDist/100)+1;
	double ret = r*maxchi2;
	if(ret>100){
		ret = 100;
	}
	return ret;
}

BOOL_T CMyFit::FindPointsOnLine3New( sEventDesc* events, int cnt,
                             double max_chi2_per_point,
                             sEventDesc* events_on_line,int& cnt_on_line,
                             double& a, double& b, double& c,
                             double& maxchi2_out,
                             sEventDesc* events_rejected, int& rejected_cnt,
									  BOOL_T bCheckVelocity, double fVelocityError,
									  double& minChi2PerEvent3Points, int ccd_idx,
									  BOOL_T bAddFromSame, eTrackCheckType_T type,
									  double max_distance /*=-1 */,
									  double MaxEstimatedPosDist /*=-1.00*/ )
{
   BOOL_T bForceSameFrame = (type==eSingleCamTrack);
//   double localMaxChi2ForPointToMatchLine = max_chi2_per_point;
//   if( type==eSingleCamTrack  && MaxChi2ForPointToMatchLine>0 ){
//     localMaxChi2ForPointToMatchLine = MaxChi2ForPointToMatchLine;
//   }
  
	double aa,bb,cc;
   int pos;
	maxchi2_out = 0;
                                                                                
   if(cnt<3){
      return FALSE;
   }

	int cnt_original=cnt;
	rejected_cnt=0;                                                                                
   cnt_on_line=0;
	sEventDesc* wrk = new sEventDesc[cnt];
	memcpy( wrk, events, cnt*sizeof(sEventDesc) );	
	int best3[3];
   double minChi2_for3=10000000000.00;
   int nOK = FindBest3Points( wrk, cnt, a, b, c, best3, minChi2_for3,
										bCheckVelocity, fVelocityError, ccd_idx,
									   bAddFromSame, type, max_distance, 
									   MaxEstimatedPosDist );

	/*if(nOK>0){
		sEventDesc best[3];
		best[0] = wrk[best3[0]];
		best[1] = wrk[best3[1]];
		best[2] = wrk[best3[2]];
		double maxDist = CalcMaxDist( best, 3 );
		max_chi2_per_point = ReCalcMaxChi2( maxDist, max_chi2_per_point );
	}*/

	minChi2PerEvent3Points = (minChi2_for3/3);
	if(minChi2PerEvent3Points<=max_chi2_per_point){
		if( gShowChi2OfAdded ){
      	printf("Chi2_of_added %d %.8f\n",3,minChi2PerEvent3Points);
      }
      
		// now use 3 points as seed and add next :
		for(register int ii=0;ii<3;ii++){
			int bestpos=best3[ii];
			events_on_line[ii] = wrk[ bestpos ];
		}
		cnt_on_line = 3;

		// best3 is sorted by FindBest3Points
		for(int ii=2;ii>=0;ii--){
			// assuming sorted 
			int bestpos=best3[ii];
			for(int jj=bestpos;jj<(cnt-1);jj++){
				wrk[ jj ] = wrk[ jj+1 ];
			}
			cnt--;
		}

      if( gDebugTracks > 0 && bCheckVelocity ){
        // dump velocity information about good track candidate :
        sEventDesc tocheck[3];
        double rx_tmp,ry_tmp;
        memcpy(tocheck,events_on_line,3*sizeof(sEventDesc));
        BOOL_T bCheckTmp = CheckVelocity( tocheck,  3, fVelocityError, 
                                          DEFAULT_MAX_DIFF_TIME, ccd_idx, type,
                                          eHistoVXRatioTo3, &rx_tmp, &ry_tmp, 
                                          eCheckVeloGoodChi2 );
        printf("DEBUG_TRACKS : GOOD_CHI2_3POINTS (%.2f,%.2f) (%.2f,%.2f) (%.2f,%.2f) chi2 = %.2f\n",
                events_on_line[0].x,events_on_line[0].y,
                events_on_line[1].x,events_on_line[1].y,
                events_on_line[2].x,events_on_line[2].y,                                            
                minChi2PerEvent3Points );                
        printf("DEBUG_TRACKS : GOOD_CHI2_3POINTS_VELOCHECK CMyFit::FindPointsOnLine3New = %d\n",bCheckTmp);
      }


		// now try adding points until only very bad left :	
		int new_pos=3;
		double* line_x = new double[cnt_original];
		double* line_y = new double[cnt_original];

		// TODO : moze potrzeba dodac sprawdzanie czy dodawany punkt nie jest za daleko od toru 
		//        dla PLANE_TRACKow oczywiscie tylko !
		while(cnt>0){
			int checkpos=cnt-1;
			BOOL_T bAdded=FALSE;
			BOOL_T bVelOK=TRUE;
			events_on_line[new_pos] = wrk[checkpos];

			if( (bAddFromSame || find_frame( events_on_line, cnt_on_line, wrk[checkpos].frame )<0) && // in case events from same image are not allowed 
				 (!bCheckVelocity || CheckVelocity( events_on_line, cnt_on_line, wrk[checkpos], fVelocityError, DEFAULT_MAX_DIFF_TIME, type, eHistoVXRatioToOld ) ) &&
             (max_distance<=0 || calc_min_dist(events_on_line, cnt_on_line,wrk[checkpos])<max_distance )	&&
             (!bForceSameFrame || find_frame( events_on_line, cnt_on_line, wrk[checkpos].frame )>=0 ) // in case only events from same image are required ( plane tracks etc )
           ){

				get_xy( events_on_line, cnt_on_line+1, line_x, line_y );
				FitLine( line_x, line_y, new_pos+1, aa, bb, cc );
				double chi2 = CalcChi2_Dist( line_x, line_y, new_pos+1, aa, bb, cc );

				if( gShowChi2OfAdded ){
					printf("Chi2_of_added %d %.8f\n",(cnt_on_line+1),chi2);
				}

				if( gDebugTracks ){
				  printf("DEBUG_TRACK : MATCH_TO_3_POINTS CMyFit::FindPointsOnLine3New %.4f\n",chi2/(cnt_on_line+1));
				}

				// skip point
				if( (chi2/(cnt_on_line+1)) <= max_chi2_per_point ){
					// add to line if good 
					new_pos++;
					cnt_on_line++;					

					// update best line coeficents :
					a = aa;
   	         b = bb;
   	         c = cc;
					bAdded = TRUE;
				}
			}
			if( !bAdded ){
				events_rejected[rejected_cnt] = events_on_line[new_pos];
				rejected_cnt++;
			}
			cnt--;
		}

		delete [] line_x;
		delete [] line_y;
	}

	
	delete [] wrk;
	return (cnt_on_line>2);
}

//
// 1/ finding 3 best points - and adding next :
//
BOOL_T CMyFit::FindPointsOnLine3( double* x_values, double* y_values, 
										 	int cnt,
										   double max_chi2_per_point, 
										   double* line_x, double* line_y, int& cnt_on_line,
										 	double& a, double& b, double& maxchi2_out,
											double* rejected_x, double* rejected_y, int& rejected_cnt )
{
	double aa,bb,cc;
	int pos;

	if(cnt<3){
		return FALSE;
	}

	cnt_on_line=0;

	double* wrk_x = new double[cnt];
	double* wrk_y = new double[cnt];


	memcpy( wrk_x, x_values, cnt*sizeof(double) );
	memcpy( wrk_y, y_values, cnt*sizeof(double) );


	//printf("POTENTIAL TRACK POINTS :\n");
	//for(register int i=0;i<cnt;i++){
	//	printf("%d %d\n",(int)wrk_x[i],(int)wrk_y[i]);
	//}    	
	//printf("\n");

	int best3pos=-1;
   double minChi2_for3=10000000000.00;
	// FindBest3Points( wrk_x, wrk_y, frames, cnt, a, b, best3pos, minChi2_for3 );
		
	for(register int i=0;i<=cnt-3;i++){
		double a3,b3,c3;
		FitLine( wrk_x+i, wrk_y+i, 3, a3, b3, c3 );
		double chi2_3 = CalcChi2_Dist( wrk_x+i, wrk_y+i, 3, a3, b3, c3 );
		if(chi2_3<minChi2_for3 || i==0){
			minChi2_for3 = chi2_3;
			best3pos = i;

			// update best coeficients here :
			a = a3;
			b = b3;
		}
	}

	if(best3pos<0){
		delete [] wrk_x;
	   delete [] wrk_y;
		return FALSE;
	}

	Assert(best3pos>=0,"No 3 points found (all=%d) ???",cnt);
	if((minChi2_for3/3)<=max_chi2_per_point){
		// now use 3 points as seed and add next :
		memcpy( line_x, wrk_x+best3pos, 3*sizeof(double) );
		memcpy( line_y, wrk_y+best3pos, 3*sizeof(double) );
		cnt_on_line = 3;
	
		int numToCopy = (cnt-(best3pos+3));
		if(numToCopy>0){
			//memcpy( wrk_x+best3pos, wrk_x+best3pos+3, numToCopy*sizeof(double) );
			//memcpy( wrk_y+best3pos, wrk_y+best3pos+3, numToCopy*sizeof(double) );
			for(register int cc=-0;cc<numToCopy;cc++){
				wrk_x[best3pos+cc] = wrk_x[best3pos+3+cc];
				wrk_y[best3pos+cc] = wrk_y[best3pos+3+cc];
			}
		}
		cnt = cnt - 3;


		// now try adding points until only very bad left :	
		int new_pos=3;
		while(cnt>0){
			int checkpos=cnt-1;
			line_x[new_pos] = wrk_x[checkpos];
			line_y[new_pos] = wrk_y[checkpos];

			FitLine( line_x, line_y, new_pos+1, aa, bb, cc );
			double chi2 = CalcChi2_Dist( line_x, line_y, new_pos+1, aa, bb, cc );

			// skip point
			cnt--;
			if( (chi2/(cnt_on_line+1)) <= max_chi2_per_point){
				// add to line if good 
				new_pos++;
				cnt_on_line++;

				// update best line coeficents :
				a = aa;
            b = bb;
			}else{
				rejected_x[rejected_cnt] = line_x[new_pos];
				rejected_y[rejected_cnt] = line_y[new_pos];
				rejected_cnt++;
			}
		}
	}


	
	delete [] wrk_x;
	delete [] wrk_y;


	/*if(cnt_on_line>2){
		printf("POINTS ON TRACK :\n");
		for(register int i=0;i<cnt_on_line;i++){
			printf("%d %d\n",(int)line_x[i],(int)line_y[i]);
		}
		printf("line = %.5f * x + %.5f\n",a,b);
	}*/

	// maxchi2_out = maxchi2;
	return (cnt_on_line>2);
}

// Dist2 - means that this is distance squared ( dist^2 ) 
// to have distance result must be sqrt ( square rooted )
double CMyFit::CalcDist2FromLine( double a, double b, double c,
								 			 double x, double y )
{
	double up = CMyMathFunc::mysqr( a*x+b*y+c );
	double bottom = ( a*a+b*b ); 
	double ret = (up/bottom);
	
	return ret;
}

// WARNING : simple chi2 instead of real formula ( I have not calculated it yet ... )
// also in ^2 
double CMyFit::CalcDistFromParabola( double a, double b, double c, 
                                     double x, double y )
{
  double par_y = a*x*x + b*x + c;
  double diff = (y-par_y);
  double chi2 = diff*diff;
  
  return chi2;
}                                     

double CMyFit::CalcDist2FromLine2Par( double a, double b,
                                      double* x, double* y, int cnt )
{
	double sum=0;
	for(register int i=0;i<cnt;i++){
		sum += CalcDist2FromLine2Par( a, b, x[i], y[i] );
	}
	return sum;
}

double CMyFit::CalcDist2FromLine2Par( double a, double b,
                                      double x, double y )
{
	return CalcDist2FromLine( a, -1 , b , x, y );
}


BOOL_T CMyFit::TryToAddNewPointsNew( sEventDesc* rejected, int rejected_cnt,
											double max_chi2_per_point,
                                 sEventDesc* line, int& cnt_on_line,
                                 double& a, double& b, double& c,
                                 double& maxchi2_out,
											int total_count,
										   BOOL_T (*fill_func)( void*, int, int, void* ),
											int cam_idx, BOOL_T bAddFromSame,
											eTrackCheckType_T type,
											double max_distance /*=-1*/ )
{
   BOOL_T bForceSameFrame = (type==eSingleCamTrack);

	double original_a = a;
	double original_b = b;
	double original_c = c;

	int added_total=0;

	maxchi2_out = 0.00;

	if(rejected_cnt>0 && cnt_on_line>2){
		double* line_x = new double[ total_count ];	
		double* line_y = new double[ total_count ];

		// re-trying to fit rest :
		int added=1;
		int cnt_on_line_save = cnt_on_line;
		int iter=0;
		while(added && rejected_cnt && iter<1000){
			added=0;

			double min_chi2=100000.00;			
			int min_pos = -1;
			for(int i=0;i<rejected_cnt;i++){
			   int evt_same_frame_pos = find_frame( line, cnt_on_line, rejected[i].frame );			     
			   if( evt_same_frame_pos >= 0 && type==eNormalTrack ){
    			   double same_frame_dist = sqrt( CMyMathFunc::mysqr(line[evt_same_frame_pos].x-rejected[i].x)+CMyMathFunc::mysqr(line[evt_same_frame_pos].y-rejected[i].y) );
	    		   if( gDebugTracks ){
	  	      	     printf("DEBUG_TRACKS : TrType=%d, event %d-(%d,%d) cannot be added to track, due to another event from frame %d already found on line a = %.8f , b = %.8f , c = %.8f  in distance = %.2f\n",
	  	      	         type,
                        rejected[i].frame,(int)(rejected[i].x),(int)(rejected[i].y),
                        rejected[i].frame,
                        a,b,c,same_frame_dist);
                    if( same_frame_dist < 20.00 ){
                      printf("DEBUG TRACKS : TrType=%d , probably part of same object !!!\n",type);
                    }
               }        			               
			   }
				if( (bAddFromSame || evt_same_frame_pos<0) &&
				    (!bForceSameFrame || evt_same_frame_pos>=0 ) &&
				    (max_distance<=0 || calc_min_dist( line, cnt_on_line, rejected[i])<max_distance ) 
				 ){
					// chosing points from frames not already belonging to track :
					double point_chi2 = CalcDist2FromLine( a, b, c, rejected[i].x, rejected[i].y );
					if(point_chi2<min_chi2){
						min_pos = i;
						min_chi2 = point_chi2;
					}				
				}
			}

			if( fill_func ){
				// 2 - is eChi2ToOld - to add to histogram when
				// adding new points to exisiting track :
				if( type==eNormalTrack ){
					(*fill_func)( &min_chi2, eChi2ToOld, cam_idx, NULL );
				}else{
					if( DoAddFromSameFrame( type ) ){
						(*fill_func)( &min_chi2, eChi2ToOldSum, cam_idx, NULL );
					}
				}
			}

			if(min_chi2<max_chi2_per_point && min_pos>=0){
				// in case point satisfying chi2 criteria found - check
				// how new fit line will match :

				if( min_chi2 < max_chi2_per_point ){				
					double aa,bb,cc;

					line[cnt_on_line] = rejected[min_pos];
					get_xy( line, cnt_on_line+1, line_x, line_y ); 

					FitLine( line_x, line_y, cnt_on_line+1, aa, bb, cc );
	   	      double chi2 = CalcChi2_Dist( line_x, line_y, cnt_on_line+1, aa, bb, cc );
					// if( (chi2/(cnt_on_line+1)) <= max_chi2_per_point){
						cnt_on_line++;
						added++;
						added_total++;
						for(int j=min_pos;j<(rejected_cnt-1);j++){
							rejected[j] = rejected[j+1];
						}
						rejected_cnt--;
						a = aa;
						b = bb;
						c = cc;
						maxchi2_out = chi2/(cnt_on_line+1);
					// }	
				}
			}
			iter++;
		}
		
		delete [] line_x;
		delete [] line_y;
	}	

	return (added_total>0);	
}

BOOL_T CMyFit::TryToAddNewPoints( double* rejected_x, double* rejected_y, int rejected_cnt,
 											 double max_chi2_per_point,
											 double* line_x, double* line_y, int& cnt_on_line,
											 double& a, double& b,  double& maxchi2_out )
{
	double original_a = a;
	double original_b = b;

	int added_total=0;

	maxchi2_out = 0.00;

	if(rejected_cnt>0 && cnt_on_line>2){
		// re-trying to fit rest :
		int added=1;
		int cnt_on_line_save = cnt_on_line;
		int iter=0;
		while(added && rejected_cnt && iter<1000){
			added=0;

			double min_chi2=100000.00;			
			int min_pos = -1;
			for(int i=0;i<rejected_cnt;i++){
				double point_chi2 = CalcDist2FromLine( a, -1, b, rejected_x[i], rejected_y[i] );
				if(point_chi2<min_chi2){
					min_pos = i;
					min_chi2 = point_chi2;
				}				
			}

			if(min_chi2<max_chi2_per_point && min_pos>=0){
				// in case point satisfying chi2 criteria found - check
				// how new fit line will match :

				if( min_chi2 < max_chi2_per_point ){				
					double aa,bb,cc;

					line_x[cnt_on_line] = rejected_x[min_pos];
	            line_y[cnt_on_line] = rejected_y[min_pos];

					FitLine( line_x, line_y, cnt_on_line+1, aa, bb, cc );
	   	      double chi2 = CalcChi2_Dist( line_x, line_y, cnt_on_line+1, aa, bb, cc );				
					// if( (chi2/(cnt_on_line+1)) <= max_chi2_per_point){
						cnt_on_line++;
						added++;
						added_total++;
						for(int j=min_pos;j<(rejected_cnt-1);j++){
							rejected_x[j] = rejected_x[j+1];
							rejected_y[j] = rejected_y[j+1];
						}
						rejected_cnt--;
						a = aa;
						b = bb;
					// }	
				}
			}
			iter++;
		}
	}	

	return (added_total>0);	
}



BOOL_T CMyFit::FitParabola( double* x_values, double* y_values, int count,
                              double& a, double& b, double& c, double& chi2 )
{
  return MyFits_FitParabola( x_values, y_values, count, a, b ,c, chi2 );
}                              

BOOL_T CMyFit::FitGauss( double minValue, double& rms, double& center, 
							 double &max, int binNo, double binWidth, 
							 int* pCountTab, int& nstep, int first_bin)
{
	if( CMyFit::m_pFitGauss ){
		return (*CMyFit::m_pFitGauss)( minValue, rms, center, max, binNo, binWidth, pCountTab, nstep,  first_bin);
	}

	if(binNo>100){		
		 return FALSE;
	}
	double T[100], Y[100], X[3], DELTAY[100], a[300], scrat[9], cx[9], r;
	int n=binNo, nr=3, nred=3, list[3];
	
	X[0]=max; X[1]=center; X[2]=rms;	// Initial fit values
	
	for(register int i=0;i<9;i++){
		scrat[i] = 0;
		cx[i] = 0;
	}
	for(register int i=0;i<100;i++){
		T[i]=0;
		Y[i]=0;
		DELTAY[i]=0;
	}
	for(register int i=0;i<300;i++){
		a[i] = 0;			
	}

	double delta=binWidth*0.5;
	for(register int ii=first_bin; ii<(first_bin+binNo); ii++)	// Generating histogram points
	{
		int i = (ii-first_bin);

		T[i] = minValue + delta;
		Y[i] = pCountTab[i];
		DELTAY[i]=1;
		delta += binWidth;
	}

	int nstep2 = nstep;	
//	lsqnon_(lsqgss_, T, Y, DELTAY, &n, &nr, &nred, list, X, cx, &r, a, scrat, &nstep2);
    printf("ERROR IN CODE : lsqnon_ function not compiled !\n");
	
	max=X[0]; 
	center=X[1]; 
	rms= fabs(X[2]);	// to be sure - there is sigma^2 in expresion so 	
							// it is possible to get sigma<0 from fit
							// which causes VERY BAD THINGS ( clusters like WHOLE FRAME )

	nstep = (int)nstep2;
	return (nstep2>0);
}


void show_points( const char* cmt, sEventDesc* list, int cnt )
{
	printf("%s",cmt);
	for(int ii=0;ii<cnt;ii++){
		double rx=0,ry=0;
		if(ii>0){
			double dx = (list[ii].x-list[ii-1].x);
			double dy = (list[ii].y-list[ii-1].y);
			double dt = (list[ii].timeUT-list[ii-1].timeUT);
			
			if( fabs(dt)>0 ){
				rx = (dx/dt);
				ry = (dy/dt);
			}
		}

 		printf("%d-(%.2f,%.2f)-(%.2f,%.2f),",list[ii].frame,list[ii].x,list[ii].y,rx,ry);
   }
   printf("\n");
}


BOOL_T CMyFit::DoAddFromSameFrame( eTrackCheckType_T track_type )
{
   return ( track_type==ePlaneTrack || track_type==eTrackOnSumedFrame || track_type==eSingleCamTrack );
}

double CMyFit::calc_dist( const sEventDesc& point1 , const sEventDesc& point2 )
{
  double dist = sqrt( (point1.x-point2.x)*(point1.x-point2.x) + (point1.y-point2.y)*(point1.y-point2.y));
  return dist;
}

double CMyFit::calc_min_dist( sEventDesc* events, int cnt, const sEventDesc& point )
{
  double min_dist = 10000000.0;
  
  for(int i=0;i<cnt;i++){
    double dist = calc_dist( events[i], point );
    
    if( dist < min_dist ){
      min_dist = dist;
    }
  }
  
  return min_dist;
}

