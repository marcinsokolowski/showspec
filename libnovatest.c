#include <stdio.h>
#include <stdlib.h>
#include "gal2hor.h"

int main( int argc, char* argv[] )
{
  double glon_deg = atof(argv[1]);
  double glat_deg = atof(argv[2]);
  double azim,alt,ra,dec;
  
  time_t ux;
  time(&ux);
  int year,mon,day,h,m;
  double s;
  get_ymd_hms_ut(ux,year,mon,day,h,m,s);
  
  double geo_long_deg=-6.7342055;
  double geo_lat_deg =37.10406944;
  gal2hor( glon_deg, glat_deg, geo_long_deg, geo_lat_deg, ux, azim, alt, ra, dec );
  
  printf("Galactic (GLON,GLAT) = (%.4f,%.4f) [deg]\n",glon_deg,glat_deg);
  printf("At unix time         = %d = %04d%02d%02d_%02d%02d%02d\n",(int)ux,year,mon,day,h,m,(int)s);
  printf("Horizontal coord     = (%.4f,%.4f) [deg]\n",azim, alt);
  printf("Equatorial coord     = (%.4f,%.4f) [deg] = (%.6f [h],%.6f [deg])\n",ra,dec,ra/15.00,dec);
}