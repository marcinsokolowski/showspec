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
#include "mathfunc.h"
// #include "mathdefs.h"
#include <stdio.h>
#include <stdlib.h>
#include "myfits.h"

static double gMult=100000000000.000;
static long double gDwa = 2.00;
static long double gJeden = 1.00;


double CMyMathFunc::m_X=0;
double CMyMathFunc::m_Y=0;
double CMyMathFunc::m_PSF=2.00;
double CMyMathFunc::m_Norm=1.00;


// ------------------------------------------------------------
// function which reduces an angle to between 0 and 360 degrees
double rev( double x )
{
   return  x - floor(x/360.0)*360.0;
}


double sind( double x )
{
	return sin((x)*DEGRAD);
}

double cosd( double x){
	return cos((x)*DEGRAD);
}

double tand( double x){
	return tan((x)*DEGRAD);
}

double asind(double x){
	return (RADEG*asin(x));
}

double acosd(double x){
	return (RADEG*acos(x));
}

double atand(double x){
	return (RADEG*atan(x));
}

double atan2d(double y,double x){
	return (RADEG*atan2((y),(x)));
}


int mysign( double val )
{
	if(val>0)
		return 1;
	if(val<0)
		return -1;
	return 0;
}

int mysign_non_zero( double val )
{
	if(val>=0)
      return 1;
	else
		return -1;
}


CMyMathFunc::CMyMathFunc()
{
	
}


/*double CMyMathFunc::Erfc( double x )
{
	double ret = 0;
	ret = 1.00 - Erf(x);
	return ret;
}*/



long double CMyMathFunc::Erfc( double x )
{
	if(x>=0){
		return ErfcPositive( x );
	}else{
		long double val = (-ErfcPositive( -x ));	
		long double dwa = 2.0000000000000000000000000000000000000000000;
		long double ret = (dwa+val);
		// printf("ret=%.50Lf val=%.50Lf\n",ret,val);
		return ret;
	}
}

long double CMyMathFunc::ErfcPositive( double x )
{
	long double dx = 0.02;
	long double dx2 = (dx*0.5);
	long double x0 = x + dx2;
	long double sum = 0;
	
	while(fabs(x0)<10.00){
		long double val = exp( -(x0*x0) );			
		sum += val;
		// printf("val(%.4Lf) = %.30Lf sum=%.30Lf\n",x0,val,sum);	
		x0 += dx;
	}
	sum = sum*dx;

	long double sqrtPI = sqrt(PI_VALUE);
	long double OneOversqrtPI = gJeden/sqrtPI;
	sum = (sum*gDwa)*OneOversqrtPI;
	// printf("Erfc(%f) = %.50Lf\n",x,sum);
	return sum;
}

/*double CMyMathFunc::mysqr( double x )
{
	return (x*x);
}*/


/*double CMyMathFunc::GetAngleFromSin( double sin_alfa_1, double x1, double y1 )
{
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
}*/

double CMyMathFunc::my_atan( double y, double x )
{
	double a = atan2( fabs(y), fabs(x) );
	if( x<0 && y>=0 )
		a = PI_VALUE - a;
	if( x<0 && y<0 )
		a = PI_VALUE + a;
	if( x>0 && y<0 )
		a = TWO_PI_VALUE - a; 

	return a;
}



double CMyMathFunc::round( double x, int places )
{
	double mult = pow(10,places);
	int x_int = x*mult;
	double ret = double(x_int)/mult;

	return ret;
}


int CMyMathFunc::find_zero_place( double (*func)( double x ), double& x_zero,
											 double x0, double x1, double delta )
{
	double len = (x1-x0);
	double f_x0 = (*func)( x0 );
	double f_x1 = (*func)( x1 );

	if(f_x0*f_x1<0){
		while(fabs(x1-x0)>delta){
			double x_new = x0+(x1-x0)/2;
			f_x0 = (*func)( x0 );
			f_x1 = (*func)( x1 );
			double f_new = (*func)( x_new );
			if( f_x0*f_new < 0 ){
				// x1 := x_new
				x1 = x_new;
			}else{
				if( f_new*f_x1 < 0 ){
					// x0 := x_new
					x0 = x_new;
				}else{
					return 0;
				}
			}
			len = (x1-x0);
		}
		x_zero = (x1+x0)*0.5;
		return 1;
	}else{
		return 0;
	}

	return 0;
}

int CMyMathFunc::calc_sqr_eq( double a, double b, double c,
                 	            double& delta, double& x1, double& x2 )
{	
	int ret=0;
	delta = b*b - 4*a*c;
	if( delta>0 ){
		printf("here1\n");
		ret=2;
		x1 = (-b+sqrt(delta))/(2*a);
		x2 = (-b-sqrt(delta))/(2*a);
	}else{
		if( delta<-1){
			return 0;
		}
		if( delta > -0.00000005 ){
			printf("here2\n");
			ret=1;
			x1 = (-b)/(2*a);
			x2 = x1;
		}
	}
	return ret;
}

int CMyMathFunc::calc_rot_z( double x, double y, double z, double angle_in_rad,
                 	           double& x_prim, double& y_prim, double& z_prim )
{
	double cos_angle = cos(angle_in_rad);
	double sin_angle = sin(angle_in_rad);
	x_prim = cos_angle*x - sin_angle*y;
	y_prim = sin_angle*x + cos_angle*y;
	z_prim = z;

	return 1;
}

int CMyMathFunc::calc_rot_y( double x, double y, double z, double angle_in_rad,
                 	           double& x_prim, double& y_prim, double& z_prim )
{
	double cos_angle = cos(angle_in_rad);
	double sin_angle = sin(angle_in_rad);
	x_prim = cos_angle*x + sin_angle*z;
	y_prim = -sin_angle*x + cos_angle*z;
	z_prim = z;

	return 1;
}

int CMyMathFunc::calc_rot_y_declin( double x, double y, double z, double angle_in_rad,
		                 	            double& x_prim, double& y_prim, double& z_prim )
{
	double cos_angle = cos(angle_in_rad);
	double sin_angle = sin(angle_in_rad);
	x_prim = cos_angle*x - sin_angle*z;
	y_prim = sin_angle*x + cos_angle*z;
	z_prim = z;

	return 1;
}

int CMyMathFunc::shift_vec( double x, double y, double z,
                         double vec_x, double vec_y, double vec_z,
                         double& x_prim, double& y_prim, double& z_prim )
{
	x_prim = x + vec_x;
	y_prim = y + vec_y;
	z_prim = z + vec_z;

	return 1;
}


/***********************************************************
* get offset from the offset_grid

OPIS :
Na podstawie 4 najblizszych komorek wyliczana jest wartosc w 
zadanej pozycji wazona odleglosciami do srodkow tych 4 
najblizszych komorek
***********************************************************/
double CMyMathFunc::offset(double* offset,int nx,int ny,
			  double xc,double yc,double xmin,double xmax,
			  double ymin,double ymax)
{
  int ix0,iy0,ix1,iy1;
  double dx,dy,fx,fy,f,f0,f1;

  dx = (xmax -xmin)/nx;
  dy = (ymax -ymin)/ny;
  fx = (xc-xmin)/dx-0.5;
  fy = (yc-ymin)/dy-0.5;
  if(fx < 0) ix0 = 0;
  else if(fx < nx-2 ) ix0 = fx;
  else ix0 = nx-2;
  if(fy < 0) iy0 = 0;
  else if(fy < ny-2 ) iy0 = fy;
  else iy0 = ny-2;

  ix1 = ix0+1;
  iy1 = iy0+1;
  fx -= ix0;
  fy -= iy0;
  f0 = (1-fx)*offset[ix0+nx*iy0]+fx*offset[ix1+nx*iy0];
  f1 = (1-fx)*offset[ix0+nx*iy1]+fx*offset[ix1+nx*iy1];
  f = (1-fy)*f0+fy*f1;
  return(f);
}

double CMyMathFunc::offset(float* offset,int nx,int ny,
			  double xc,double yc,double xmin,double xmax,
			  double ymin,double ymax)
{
  int ix0,iy0,ix1,iy1;
  double dx,dy,fx,fy,f,f0,f1;

  dx = (xmax -xmin)/nx;
  dy = (ymax -ymin)/ny;
  fx = (xc-xmin)/dx-0.5;
  fy = (yc-ymin)/dy-0.5;
  if(fx < 0) ix0 = 0;
  else if(fx < nx-2 ) ix0 = fx;
  else ix0 = nx-2;
  if(fy < 0) iy0 = 0;
  else if(fy < ny-2 ) iy0 = fy;
  else iy0 = ny-2;

  ix1 = ix0+1;
  iy1 = iy0+1;
  fx -= ix0;
  fy -= iy0;
  f0 = (1-fx)*offset[ix0+nx*iy0]+fx*offset[ix1+nx*iy0];
  f1 = (1-fx)*offset[ix0+nx*iy1]+fx*offset[ix1+nx*iy1];
  f = (1-fy)*f0+fy*f1;
  return(f);
}



/**************************************************************************
* smooth (with median) function
*   ml(xl[i],yl[i]), i=0...n-1;
*   into the (float)matrix off[nc,nc]
*   minimum nuber of medianed points: nmin
*   maximum radius for collecting points: rmax grid_points

OPIS :
Ta funkcja dzieli chip na nc x nc ( default : 32 x 32 ) , komorek , 
potem w kazdej znajduje gwiazy podane w INPUT, nastepnie szuka w 
sasiednich komorkach, az do promienia rmax gwiazd, chyba ze osiagnie 
nmin ( minimalna wymagana liczbe gwiazd ) , wtedy przestaje poszukiwac, 
przestaje takze jesli przekroczy promien poszukiwan rmax. 
**************************************************************************/
int CMyMathFunc::smooth(double* xl,double* yl,double* ml,
	   int n,
	   double xmin,double xmax,double ymin,double ymax,
	   float** off,int nc,int nmin,int rmax,
	   double* sm,double* ss, int bShowMap )
{
  int i,j,k,ii,jj,m,ix,iy,in,iss;
  int *ic=NULL,**id=NULL,mk,r,*ind;
  double *md, smm, sss, d;

  // image is divided in 32x32 squares (each ~ 2030/32 ~ 67 pixeli ) 
  // preparation of cells, id is table 2D to quickly find given star
  // it is 2D index in fact , ic table is counter for each cell 
  // counting number of stars in cell
  *off = (float *)malloc(nc*nc*sizeof(float));
  id = (int **) calloc(nc*nc,sizeof(int *)); // calloc fills with 0 
  ic = (int *) calloc(nc*nc,sizeof(int));    // calloc fills with 0
  for(i=0; i<n; i++){
    if(ml[i] > 30)continue;
    ix = (xl[i]-xmin)/(xmax-xmin)*nc;
    iy = (yl[i]-ymin)/(ymax-ymin)*nc;
    if(ix < 0) ix = 0;
    if(ix > nc-1) ix = nc-1;
    if(iy < 0) iy = 0;
    if(iy > nc-1) iy = nc-1;
    in = iy*nc+ix;
    id[in] = (int*)realloc(id[in], (ic[in]+1)*sizeof(int*));
    *(id[in]+ic[in]) = i;
    ic[in]++;
  }
  iss = 0;
  smm = sss = 0.;

  // Dla kazdej komorki wybieramy teraz co najmniej nmin gwiazd w jej 
  // otoczeniu ( max rmax) , i zapisujemy sobie je w tablicy 
  // ich poprawki beda uzyte do najlepszej sredniej poprawki w danej komorce 
  for(j=0; j<nc; j++){
    for(i=0; i<nc; i++){
      mk = 0;
      md = NULL;
      r=0;
      while ( mk < nmin){
        for(jj=j-r; jj<=j+r; jj++){
          if(jj < 0 || jj >= nc)continue;
          for(ii=i-r; ii<=i+r; ii++){
            if(ii < 0 || ii >= nc)continue;
//            d = sqrt ( (double)((ii-i)*(ii-i)+(jj-j)*(jj-j)));
//            if( r> 0 && d > r)continue;

				// check only borders, skip inside ( checked earlier )
            if(jj>j-r && jj <j+r && ii>i-r && ii<i+r) continue;

				// add all stars in given radius to list of stars for normalization
				// of current cell (one of 32x32 cells to which CCD was divided )
            in = ii + jj*nc;
            for(k=0; k<ic[in]; k++){
              if(mk%100 == 0)md = (double*)realloc(md, (mk+100)*sizeof(double));
              md[mk++] =  ml[*(id[in]+k)];
            }
          }
        }
        r++;
        if(r > rmax)break;
      }
      if(mk == 0){
        (*off)[i+j*nc] = 0.;
        iss++;
      }else{
        ind = (int *)malloc(mk*sizeof(int));
        dsortindx(mk,md,ind);
        smm += ((*off)[i+j*nc] = md[ind[mk/2]]);
        sss += (md[ind[(int)(0.83*mk)]] - md[ind[(int)(0.17*mk)]])/2;
        iss++;
        free(ind);
      }
      if(md)free(md);
    }
  }  

  if( bShowMap ){
  	 printf("------------------------------ STARS IN CELLS ------------------------------\n");
    for(j=(nc-1); j>=0; j--){
  	 	for(i=0; i<nc; i++){
  	 		printf("%.4d ",ic[i+j*nc]);
		}
		printf("\n");
	 }
	 printf("----------------------------------------------------------------------------\n");
  }

  if(iss) sss /= iss, smm /= iss;
  *ss = sss;
  *sm = smm;

  for(i=0; i<nc*nc; i++)if(id[i])free(id[i]);
  free(id);
  free(ic);
  for(i=0; i<nc*nc; i++) (*off)[i] -= smm;
  return(0);
}

void CMyMathFunc::dsortindx(int n,double* arrin,int* indx)
{
  int l,ir,i,j,indxt;
  double q;
 
  if(n<=0)return;
  for ( i=0; i<n; i++ ) indx[i]=i;
  if ( n <= 1 ) return;
  l=n/2;
  ir=n-1;
  while ( 10 ) {
    if ( l > 0 ) {
      l--;
      indxt=indx[l];
      q=arrin[indxt];
    } else {
      indxt=indx[ir];
      q=arrin[indxt];
      indx[ir]=indx[0];
      ir--;
      if ( ir == 0 ) {
        indx[0]=indxt;
        return;
      }
    }
    i=l;
    j=l+l+1;
    while ( j <= ir ) {
      if ( j < ir )
        if ( arrin[indx[j]] < arrin[indx[j+1]] ) j++;
      if ( q < arrin[indx[j]] ) {
        indx[i]=indx[j];
        i=j;
        j+=j+1;
      } else {
        j=ir+1;
      }
    }
    indx[i]=indxt;
  }
}

double CMyMathFunc::bilinear_interpol( double x, double y, 
										double h1, double h2, double h3, double h4 )
{
	double h0 = h1 + (h2-h1)*x + (h3-h1)*y + (h1-h2-h3+h4)*x*y;
	return h0;	
}

int CMyMathFunc::calc_line( double x0, double y0, double x1, double y1,
                     double& a, double& b, double& c,
                     double& tx,double& ty )
{
	a = (y0-y1);
	b = (x1-x0);
	c=(y1*x0 - y0*x1);
	
//	double r=sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) );
	tx = (x1-x0);
	ty = (y1-y0);

	double test1 = a*x0+b*y0+c;
	double test2 = a*x1+b*y1+c;

	if( fabs(test1) > 0.00001 || fabs(test2)>0.00001 ){
		printf("ERROR in CMyMathFunc::calc_line (%f,%f)-(%f,%f) : a=%f,b=%f,c=%f\n",x0,y0,x1,y1,a,b,c);
		return 0;
	}

	return 1;
}


int CMyMathFunc::calc_line_crossing( double a1,double b1, double c1,
                 		                double a2,double b2, double c2,
                     	             double& x0, double& y0 )
{
	x0 = ( b2*c1 - b1*c2 ) / ( a2*b1 - a1*b2 );
	y0 = - ( a1*x0 + c1 ) / b1;
	return 1;
}

double CMyMathFunc::line_value( double a,double b, double c, double x )
{
	double ret = -(a*x+c)/b;
	return ret;
}

int CMyMathFunc::is_above( double a,double b, double c, double x0, double y0 )
{
	if( b != 0 ){
		double y = -(a*x0+c)/b;
		if( y0 >= y ){
			return 1;
		}
	}else{
		double x = -(c/a);
		if( x0 <= x )
			return 1;
	}
	return 0;
}

int CMyMathFunc::is_below( double a,double b, double c, double x0, double y0 )
{
	if( b != 0 ){
		double y = -(a*x0+c)/b;
		if( y0 <= y ){
			return 1;
		}
	}else{
		double x = -(c/a);
		if( x0 <= x )
			return 1;
	}
	return 0;
}

int CMyMathFunc::is_inside( double x0, double y0,
									 double a1,double b1, double c1,
   	                      double a2,double b2, double c2,
      	                   double a3,double b3, double c3,
         	                double a4,double b4, double c4 )
{
	if( is_above( a1, b1, c1, x0, y0 ) ){
		if( is_above( a2, b2, c2, x0, y0 ) ){
			if( is_below( a3, b3, c3, x0, y0 ) ){
				if( is_below( a4, b4, c4, x0, y0 ) ){
					return 1;
				}
			}
		}
	}
	return 0;
}

double CMyMathFunc::mysqrt( double x )
{
	return ::sqrt( x );
}

double CMyMathFunc::GlobalGauss( double x, double y )
{
	double n = (m_Norm/(2*3.1415*m_PSF*m_PSF));
	double valG = n*exp( -( (x-m_X)*(x-m_X) + (y-m_Y)*(y-m_Y))/(2.00*m_PSF*m_PSF) );

	return valG;
}

double CMyMathFunc::GaussIntegral( double x0, double y0, double x1, double y1 )
{
	return CMyFit::GaussIntegral( x0, y0, x1, y1 );
}

double CMyMathFunc::ParLenFunc( double p )
{
	double ret = p*sqrt(1.0+p*p)+asinh(p);
	return ret;
}


double CMyMathFunc::ParLen( double a, double b, double c,
                             double x0, double x1 )
{
	double s=0.00;
	if( a != 0 ){
		s = (1.00/(4.00*a))*(ParLenFunc(2*a*x1+b)-ParLenFunc(2*a*x0+b));
	}
	
	return s;
}                             
            
// RETURNS :              
// -1 : when cannot find position where parabola lenght equals len
//        because it is outside CCD 
int CMyMathFunc::FindParStartEnd( double a, double b, double c,
                                  double x_last, double s_in,
                                  double vx, 
                                  double& x0_out, double& x1_out, double& s_out,
                                  double err )
{
  int max_loop=1000;

  if( vx > 0 ){
    // moving right : x0 - constant , changing x1 :
    double x1_start = x_last; // lower limit for x1
    double x1_end = x_last + 2*s_in; 	     // upper limit for  x1
    
    // finding "zero place" with bisection method 
    int i=0;
    double x1_test = (x1_start+x1_end)/2.00;    
    double len = ParLen( a, b, c, x_last, x1_test );
    while ( abs(len-s_in) > err && i<max_loop ){ // was abs(x1_start-x1_end) > err
      x1_test = (x1_start+x1_end)/2.00;
      len = ParLen( a, b, c, x_last, x1_test ); // calculate lenght of parabola 
      
      if( len < s_in ){
        // too short -> x1 > x1_test :
        x1_start = x1_test;
      }else{
        // too long -> x1 < x1_test :
        x1_end = x1_test;
      }
      i++;
    }    
    
    x0_out = x_last;
    x1_out = x1_test;
    s_out = ParLen( a, b, c, x0_out, x1_out );
  }else{
    // moving left : x1 - constant , changing x0 :
    double x0_start = x_last - 2*s_in; // lower limit for x0
    double x0_end = x_last;         // upper limit for x0
    
    // finding "zero place" with bisection method 
    int i=0;
    double x0_test = (x0_start+x0_end)/2.00;
    double len = ParLen( a, b, c, x0_test, x_last );
    while ( abs(len-s_in) > err && i<max_loop ){ // was abs(x0_start-x0_end) > err
      x0_test = (x0_start+x0_end)/2.00;
      len = ParLen( a, b, c, x0_test, x_last );      
      
      if( len < s_in ){
        // too short -> x1 > x1_test :
        x0_end = x0_test;
      }else{
        // too long -> x1 < x1_test :
        x0_start = x0_test;
      }
      i++;
    }    

    x0_out = x0_test;    
    x1_out = x_last;
    s_out = ParLen( a, b, c, x0_out, x1_out );
  }


  return (x0_out>0 && x0_out<2048 && x1_out>0 && x1_out<2048);
}
