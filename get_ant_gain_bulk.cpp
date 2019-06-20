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
    char optstring[] = "hg";
    int opt;

    while ((opt = getopt(argc, argv, optstring)) != -1) {
        switch (opt) {
            case 'h':
                usage();
                break;

            case 'g':
                CAntPatternParser::m_bEnableInterpolation = 1;
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
  printf("Interpolation enabled = %d\n",CAntPatternParser::m_bEnableInterpolation);
  printf("#################################################\n");
}

int main(int argc,char* argv[])
{
  string infile = argv[1];
  string outfile = argv[2];
  double freq_mhz = 140.0;

  // parse command line :  
  parse_cmdline(argc-2,argv+2);
  print_parameters();
    

  char buff[1024];  
  FILE* inf = fopen(infile.c_str(),"r");
  FILE* outf = fopen(outfile.c_str(),"w");

  gAntennaPattern = new CAntPatternParser();
  if( gAntennaPattern->Read( gAntennaPatternFile.c_str() ) <= 0 ){
     printf("ERROR : could not antenna pattern file %s -> cannot continue\n",gAntennaPatternFile.c_str());
     exit(-1);
  }                  
  
  while (1) {
     if(fgets(buff,1024,inf)==0)
        break;
     if(buff[0]=='#')
        continue;
     if(strstr(buff,"nan"))
        continue;

     double zenith_dist,power,azim;
     int ncols = sscanf(buff,"%lf\t%lf\t%lf\n",&zenith_dist,&power,&azim);
  
     double phi_deg = azim2phi(azim);  
     printf("(Azim,Alt) = (%.2f,%.2f) [deg] -> (Phi,Theta) = (%.2f,%.2f) [deg] , power = %e\n",azim,(90-zenith_dist),phi_deg,zenith_dist,power);

     double gain = gAntennaPattern->GetGain(freq_mhz,phi_deg,fabs(zenith_dist));   
     printf("OUTLINE : %.8f %.8f %.8f\n",zenith_dist,gain,azim);
     
     fprintf(outf,"%.8f %.8f %.8f\n",zenith_dist,gain,azim);
  }

  fclose(inf);
  fclose(outf);
}  

