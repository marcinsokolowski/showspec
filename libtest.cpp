#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "libnova_interface.h"
#include "bg_globals.h"

int main(int argc,char* argv[])
{
	double uxtime = get_dttm();
	if( argc > 1 && strcmp(argv[1],"-") ){
		uxtime = atof(argv[1]);		
	}
	if( argc > 2 && strcmp(argv[2],"-") ){
		geo_long = atof(argv[2]);
	}
	
	double jd;
	double sid_time_h = get_local_sidereal_time( (double)uxtime, jd );

	printf("Sidereal time at GeoLong=%.8f [deg] = %.8f [h] = %.8f [deg]\n",geo_long,sid_time_h,sid_time_h*15.00);
}
