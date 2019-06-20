#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "libnova_interface.h"
#include "bg_globals.h"

void usage()
{
   printf("radec2azh UXTIME [SITE default Muresk]\n");
   printf("Sites : mwa(or mro) , ebo (EBO, ebo201309, ebo201312), wond\n");
   exit(-1);
}

int main(int argc,char* argv[])
{
	double uxtime = get_dttm();
	if( argc<=1 || strncmp(argv[1],"-h",2)==0 ){
	   usage();
	}
	
	if( argc > 1 && strcmp(argv[1],"-") ){
		uxtime = atof(argv[1]);		
	}
	if( argc > 2 && strcmp(argv[2],"-") ){
      set_geo_location( argv[2] );
	}

	double ra = 156.87722403644997;  
   double dec = -26.703319000000004;
   if( argc > 3 && strcmp(argv[3],"-") ){	   
      ra = atof(argv[3]);
   }
   if( argc > 4 && strcmp(argv[4],"-") ){	   
      dec = atof(argv[4]);
   }
	
	double jd;
//	void radec2azh( double ra, double dec, time_t unix_time, double& out_azim, double& out_alt );
   double azim,alt;
   radec2azh( ra,dec, uxtime, geo_long, geo_lat, azim, alt);
   
   printf("(AZIM,ALT) = ( %.8f , %.8f )\n",azim,alt);

}
