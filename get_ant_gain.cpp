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

// frequency range 
int gDumpFreqRange=0;
double gStartFreqMHZ=0;
double gEndFreqMHZ=0;
double gStepFreqMHZ=0;

void usage()
{
  printf("get_ant_gain FREQUENCY[MHz] PHI[deg] THETA[deg] -f ANT_PATTERN_FILE.dat -b FREQ_START_MHZ -e FREQ_END_MHZ -s FREQ_STEP_MHZ\n");
  printf("-f ANT_PATTERN_FILE.dat - FEKO pattern file [default %s]\n",gAntennaPatternFile.c_str());
  
  exit(0);
}

void parse_cmdline(int argc, char *argv[]) {
    char optstring[] = "hf:b:e:s:";
    int opt;

    while ((opt = getopt(argc, argv, optstring)) != -1) {
        switch (opt) {
            case 'f':
                if( optarg ){
                    gAntennaPatternFile = optarg;
                }
                break;

            case 'b':
                if( optarg ){
                    gStartFreqMHZ = atof(optarg);
                    gDumpFreqRange=1;
                }
                break;

            case 'e':
                if( optarg ){
                    gEndFreqMHZ = atof(optarg);
                    gDumpFreqRange=1;
                }
                break;

            case 's':
                if( optarg ){
                    gStepFreqMHZ = atof(optarg);
                    gDumpFreqRange=1;
                }
                break;

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

double freq_mhz = 0;
double phi_deg = 0;
double zenith_dist = 0;

void print_parameters()
{
  printf("#################################################\n");
  printf("PARAMTERS :\n");
  printf("#################################################\n");
  printf("(Phi,Theta)  = (%.2f,%.2f) [deg]\n",phi_deg,zenith_dist);
  printf("Pattern file = %s\n",gAntennaPatternFile.c_str());
  printf("#################################################\n");
}

int main(int argc,char* argv[])
{
  if( argc < 2 || strncmp(argv[1],"-h",2)==0 ){
     usage();
  }

  // parse command line :  
  parse_cmdline(argc-3,argv+3);
  print_parameters();

  freq_mhz = atof(argv[1]);
  phi_deg = atof(argv[2]);
  zenith_dist = atof(argv[3]);                   
             

  if( strlen(gAntennaPatternFile.c_str()) ){
     gAntennaPattern = new CAntPatternParser();
     printf("Reading antenna pattern file |%s| ...\n",gAntennaPatternFile.c_str());fflush(stdout);
     if( gAntennaPattern->Read( gAntennaPatternFile.c_str() ) <= 0 ){
        printf("ERROR : could not antenna pattern file %s -> cannot continue\n",gAntennaPatternFile.c_str());
        exit(-1);
     }else{  
        if( gDumpFreqRange <= 0 ){
           double gain = gAntennaPattern->GetGain(freq_mhz,phi_deg,zenith_dist);   
           printf("Gain = %.8f\n",gain);
        }else{
           for(double freq_mhz=gStartFreqMHZ;freq_mhz<=gEndFreqMHZ;freq_mhz+=gStepFreqMHZ ){
              double gain = gAntennaPattern->GetGain(freq_mhz,phi_deg,zenith_dist);
              printf("%.2f %.8f\n",freq_mhz,gain);
           }
        }
     }
  }
}  

