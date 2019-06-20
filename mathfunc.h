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
#ifndef _MATH_FUNC_H__
#define _MATH_FUNC_H__


#include <math.h>
#include "mathdefs.h"

int mysign( double val );
int mysign_non_zero( double val );

// trygonomatric on DEGREES :
#define RADEG       (180.0/PI_VALUE)
#define DEGRAD      (PI_VALUE/180.0)

double rev( double x );
double sind( double x );
double cosd( double x);
double tand( double x);
double asind(double x);
double acosd(double x);
double atand(double x);
double atan2d(double y,double x);


class CMyMathFunc 
{
public:
	CMyMathFunc();
	
	// static double Erf( double x );
	static long double Erfc( double x );
	static long double ErfcPositive( double x );
	
	static inline double mysqr( double x ){ return (x*x); }

	static double my_atan( double y, double x );

	static inline double GetAngleFromSin( double sin_alfa_1, double x1, double y1 ){
	   double alfa_1_tmp = asin( sin_alfa_1 );
	   double alfa_1 = alfa_1_tmp;

   	if(x1>=0 && y1>=0)
      	alfa_1 = alfa_1_tmp;
	   if(x1>=0 && y1<0)
   	   alfa_1 = (2*PI_VALUE) - alfa_1;
	   if(x1<0 && y1>=0)
   	   alfa_1 = PI_VALUE - alfa_1;
	   if(x1<0 && y1<0)
   	   alfa_1 = PI_VALUE + alfa_1;


		return alfa_1;
	}


	static double round( double x, int places );


	// 0 - not found 
	// 1 - found
	static int find_zero_place( double (*func)( double x ), double& x_zero,
								double x0, double x1, double delta=0.01  );

	static int calc_sqr_eq( double a, double b, double c, 
									double& delta, double& x1, double& x2 );

	static int calc_rot_z( double x, double y, double z, double angle_in_rad,
								  double& x_prim, double& y_prim, double& z_prim );
	static int calc_rot_y( double x, double y, double z, double angle_in_rad,
								  double& x_prim, double& y_prim, double& z_prim );
	static int calc_rot_y_declin( double x, double y, double z, double angle_in_rad,
								  double& x_prim, double& y_prim, double& z_prim );

	static int shift_vec( double x, double y, double z,
								 double vec_x, double vec_y, double vec_z,
								 double& x_prim, double& y_prim, double& z_prim );						

   /* OPIS :
   Na podstawie 4 najblizszych komorek wyliczana jest wartosc w 
   zadanej pozycji wazona odleglosciami do srodkow tych 4 
   najblizszych komorek
   	offset - 2 te same funkcje roznia sie tylko pierwszym paramterem float vs double   
   */
	static double offset(double* offset,int nx,int ny,
	              double xc,double yc,double xmin,double xmax,
	              double ymin,double ymax);
	                                                    
	static double offset(float* offset,int nx,int ny,
	              double xc,double yc,double xmin,double xmax,
	              double ymin,double ymax);

	// nazwa jest mylaca, tak naprawde ta funkcja znajduje na podstawie 
	// poprawek dla listy gwiazd katalogowych xl,yl,ml,n poprawki
   // w komorkach ( caly chip dzielony na nc x nc komorek )
   // i wyliczana jest mediana z poprawek w kazdej komorce
   // oraz srednia poprawka na calym chipie
   // Ta funkcja dzieli chip na nc x nc ( default : 32 x 32 ) , komorek , 
   // potem w kazdej znajduje gwiazy podane w INPUT, nastepnie szuka w 
   // sasiednich komorkach, az do promienia rmax gwiazd, chyba ze osiagnie 
   // nmin ( minimalna wymagana liczbe gwiazd ) , wtedy przestaje poszukiwac, 
   // przestaje takze jesli przekroczy promien poszukiwan rmax.    
   // OUTPUT : srednia poprawka w kazdej z nc x nc komorek, srednia poprawka
   // na calym chipie, rozrzut poprawek na calym chipie
	static int smooth( double* xl,double* yl,double* ml, int n,
	            double xmin,double xmax,double ymin,double ymax,
	            float** off,int nc,int nmin,int rmax,
	            double* sm,double* ss, int bShowMap=0 );

	static void dsortindx(int n,double* arrin,int* indx);    

	static double bilinear_interpol( double x, double y,
	                        double h1, double h2, double h3, double h4 );
	                              

	// calculates line passing through 2 points, in ax+by+c=0 form
	// and parametric : [ (x0+tx*p),(y0+ty*p) ]
	static int calc_line( double x0, double y0, double x1, double y1,
							double& a, double& b, double& c,
							double& tx,double& ty );

	static int calc_line_crossing( double a1,double b1, double c1,
											 double a2,double b2, double c2,
											 double& x0, double& y0 );

	static double line_value( double a,double b, double c, double x );
											 
	static int is_inside( double x0, double y0, 
								 double a1,double b1, double c1,
	                      double a2,double b2, double c2,											 	
	                      double a3,double b3, double c3,
	                      double a4,double b4, double c4 );

	static int is_above( double a,double b, double c, double x0, double y0 );
	static int is_below( double a,double b, double c, double x0, double y0 );
	
	static double mysqrt( double x );

	// calculatin gauss :
	static double m_X;
	static double m_Y;
	static double m_PSF;
	static double m_Norm;
	
	static double GlobalGauss( double x, double y );
	static double GaussIntegral( double x0, double y0, double x1, double y1 );

	static double ParLenFunc( double p );
	static double ParLen( double a, double b, double c,
								 double x0, double x1 );

	static int FindParStartEnd( double a, double b, double c,
                                  double x0_last, double s_in,
                                  double vx, 
                                  double& x0_out, double& x1_out, double& s_out,
                                  double err=0.1 );

};


#endif
