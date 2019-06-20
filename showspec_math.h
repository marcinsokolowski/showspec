#ifndef _SHOWSPEC_CONST_H__
#define _SHOWSPEC_CONST_H__

const double k_b=1.3806488E-23; // J/K 
const double   c=299792458.00;  // m/s 
const double rad2deg_const = (180.00/M_PI); // RADIAN -> DEGREES 
const double deg2rad_const = (M_PI/180.00); // DEGREES -> RADIAN 

inline double sqr(double val)
{
   return (val*val);
}

inline double fabs_inline(double val)
{
   if( val < 0 ){
      return (-val);
   }
   
   return val;
}

#endif
