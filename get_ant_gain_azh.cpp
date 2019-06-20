#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <healpix_base.h>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include "bg_globals.h"
#include "bg_components.h"
#include "bg_units.h"
#include "antenna_pattern.h"

using namespace std;

// options :
string gAntennaPatternFile=DEFAULT_ANT_PATTERN_FILE;

eAntennaPatternType gAntennaPatternType=eAntPattDipolOverGndScreen;
CAntPatternParser* gAntennaPattern=NULL;

void usage()
{
  printf("get_ant_gain FREQ PHI[deg] THETA[deg]\n");
  
  exit(0);
}

void parse_cmdline(int argc, char *argv[]) {
    char optstring[] = "h";
    int opt;

    while ((opt = getopt(argc, argv, optstring)) != -1) {
        switch (opt) {
            case 'h':
                usage();
                break;

            default:
                fprintf(stderr,"Unknown option %c\n",opt);
                usage();
        }
    }
}


double azim2phi( double azim_deg )
{
   double alpha_deg = 20.00; // 20-30 deg rotation of South direction from the X-axis 
   
   double phi_simul_deg = azim_deg + (90.00 - alpha_deg);
   if( phi_simul_deg >= 360 ){
      phi_simul_deg = phi_simul_deg - 360;
   }
   
   return phi_simul_deg;                     
}

void print_parameters()
{
  printf("#################################################\n");
  printf("PARAMTERS :\n");
  printf("#################################################\n");
  printf("#################################################\n");
}

int main(int argc,char* argv[])
{
  if( argc < 2 ){
     usage();
  }

  double freq_mhz = atof(argv[1]);
  double azim = atof(argv[2]);
  double zenith_dist = atof(argv[3]);
  
  double phi_deg = azim2phi(azim);
  
  printf("(Azim,Alt) = (%.2f,%.2f) [deg] -> (Phi,Theta) = (%.2f,%.2f) [deg]\n",azim,(90-zenith_dist),phi_deg,zenith_dist);

  // parse command line :  
  parse_cmdline(argc,argv);
  print_parameters();

  if( strlen(gAntennaPatternFile.c_str()) ){
     gAntennaPattern = new CAntPatternParser();
     if( gAntennaPattern->Read( gAntennaPatternFile.c_str() ) <= 0 ){
        printf("ERROR : could not antenna pattern file %s -> cannot continue\n",gAntennaPatternFile.c_str());
        exit(-1);
     }else{  
        double gain = gAntennaPattern->GetGain(freq_mhz,phi_deg,zenith_dist);   
        printf("Gain = %.8f\n",gain);
     }
  }
}  

