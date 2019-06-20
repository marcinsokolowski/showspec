// PROGAM WRITTERN BY MARCIN SOKOLOWSKI 
//  It integrates whole sky according to given antenna pattern :
// WARNING :
//   - if antenna pattern is taken from FEKO output CSV files and G_tot is there as last column, it contains efficiency and may be used as it is 
//       BY default it is normalized to integral of G_tot over whole 4PI -> meaning that efficiency is lost and should be accounted for later
//       with option -q norm_4pi=1 - divided by 4PI intstead of integral -> antenna efficiency is already accounted for 
//   TODO : compare to methods !

// version 1.00 - before big optimization changes 
// version 2.00 - after speed changes (sun,antenna pattern caching)
#define SHOWSPEC_VERSION "Version2.00"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include "MollweideSkyMap.h"
#include "gal2hor.h"
#include "bg_globals.h"
#include "libnova_interface.h"
#include "ce300e.h"
#include "bg_components.h"
#include "bg_units.h"
#include "antenna_pattern.h"
#include "weather_station.h"
#include "skyinfo.h"
#include "bg_fits.h"

#include "random.h"
#include "myfile.h"
#include "myparser.h"
#include "mystring.h"
#include "mystrtable.h"
#include "basestructs.h"

// HEALPIX includes :
#include <healpix_base.h>
#include <healpix_map.h>
#include <arr.h>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>

// signal :
#include <signal.h>

#include "skymap_cache.h"
#include "showspec_math.h"


using namespace std;

CSkyMapCache gSkyMapCache;
                  
// options :
string infile_list="spec_files_list";
int gDebug=0;
double start_ux_time=get_dttm();
int    time_interval=0; // was 1800 1/2 hour
double  dt=1800; // 1/2 hour 
string gUxListFile; // list of ux times to generate, according to real data 
//double geo_long=127.885437;
//double geo_lat=-31.826232;
double gFreqStart=0;     // MHz 
double gFreqEnd=1000.00; // MHz 

// Default location : Muresk : 
// double geo_long=116.68609722;
// double geo_lat=-31.74625000;
// string gSiteName="Muresk";
// Caiguna : (125.47768611,-32.34921667)

string gOutPutDir="./";
int gDoSaveAllMaps=0;
int gSaveAntPatternDB=1;
mystring gInDir="$(BIGHORNS)/software/analysis/gsm/10Mhz_to_500Mhz_step10Mhz/";
// string gInDir="/opt/caastro/bighorns/software/analysis/gsm/10Mhz_to_500Mhz_step10Mhz/";
string gOutSpecFileName;
int gOnlyAboveHorizont=0; // set -1 to turn off horizon cut 
double gThetaMax=-10000.00;
int gSaveSkyMapFits=0;
string gMapFileBaseName="";
double gSaveSkyMapAtFreq=-1000;
int gSaveFitsFile=0;
int gZenithProjRadius=400;
string gZenithProjBaseName="zenith_projection";
string gFitsFile;

// autocontinue option :
int gAutoContinue=0;
double gLastUxTimeProcessed=-1000.00;
double gLastFreqProcessed=-1000.00;

// sensor weather file : UXTIME TEMP HUM DEWPOINT 
string gAntEffFile;
int gAddAntEffCorrection=0;
vector<cValue> gAntEffTable;
CWeatherStation gWeather;
double gT_amb=290.00;

// included foregrounds :
int gAddGalaxy=1;
int gAddSun=1;
int gAddSources=0;
double gSourceRadiusArcSec=233.00; 
double gExcludeSourcesRadiusArcSec=-1; // size of HealPIX pixel ~232.6315 arcsec 
double gExternalNoiseT=0;

// SUN :
//   double sun_T_b = 600000.00; // from Claude Mercier and Gilbert Chambe 2009 paper 
//   double sun_radius_arcsec = 1818.00/2.00;
// Lantos & Avignon - Sun ~ 60000Jy 
// Aaron Chippendale :
double sun_T_b = 1000000.00; // Lantos et all 1975 etc - difficult to dind one good value 
double sun_radius_arcsec = (31.00/2.00)*60.00; // Diameter ~ 31 arcmin 
double sun_radius_deg = sun_radius_arcsec / 3600.00;         

string gObjectListFile;
vector<cObjectInfo> gObjectList;


// HEALPIX global variables :
// fill coordinates just once and then re-use (OPTIMIZATION) :
long int nside = Healpix_Base::npix2nside(12*512*512);
long int norder = Healpix_Base::nside2order( nside );   
Healpix_Base* pHealpix=NULL;
MollweideSkyMap gMap(100);
long double dOMEGA = 0;

// Ionospheric effects :
enum eIonoAbsModelType {eNoIonoAbsorption=0,eNoIonoHarish=1,eNoIonoAbhiDatta=2,eIonoLossSimple=3};
int gTestTskyWithTau=0;
int gIncludeIonosphere=0;
eIonoAbsModelType gIncludeIonoAbsorption=eNoIonoAbsorption;
int gIncludeIonoElectronTemp=0;
double gPlasmaFrequency = 5.00; // MHz 
double gIonoAtt0 = 1.00;
int gIonoLossWithPathLength = 1;
double gTypicalLossAt100MHz = 0.01; // [dB] at nighttime

// IONOSPHERIC REFRACTION in F-layer
enum eIonoRefractionModelType {eNoIonoRefraction,eIonoRefractionHarishAbhiData=1,eIonoRefractionSimple=3,eIonoRefractionSimpleWithThetaDep=4};
eIonoRefractionModelType gIncludeIonoRefraction=eNoIonoRefraction;
double gTEC_total = 10e16; // x 10^16 [e/m^2] - TEC as measured by GPS satellites 
double gFlayerHeight = 300e3; // [m]
double gFlayerD = 200e3; // [m]
double gTypicalRefractionAngleArcMin = 1; // arcmin at 100 MHz 
double gTypicalRefractionZenithAngle = 45.00;

// IONOSPHERIC LOSSES - as per Harish's paper :
// D-layer parameters 
double gHeightD = 75e3; // [m] = 75km 60-90km
double gDeltaHeightD = 30e3; // [m] 30km
double gRadiusEarth = 6371e3; // [m] = 6300 km
double gCollisionFreq = 10 * 1e6; // 10 MHz
// double gPlasmaFreq    = 9.00*sqrt(2.5*1e8); // in Hz
double gN_e = 0.44e8; // 1/m^3
double gPlasmaFreq    = 6.00*1e4; // this is what matches Harishe's plots - why ???
double gLightSpeed = 300e6; // m/s
double gElectronTemp = 800; // K 
double gDLayerTEC_fraction = 0.0008; // 8e-4 - D-layer contribtion to TEC (from A.Datta paper)

// TROPOSPHERIC REFRACTION:
int gIncludeTropoRefraction=0;
double gAtmPressure=1000.00; // millibar = hPa 

CAntPattern gAntennaPattern;
int gNorm4PI=0;

// fast but not so precise version :
int gRunFastVersion=0;

// intergation info :
double gIntTime=(200*0.2730667); // = 54.6133400seconds 

// program execution and timing :
time_t program_start = get_dttm();

// some local functions :
double parse_ni( const char* infile );
double integrate_galactic_spectrum_FAST( const char* infile , double unix_time, ofstream* ofile, double& out_pattern_sum, double& out_sun_sum, int first_freq  );
void AutoDetectExtensions( vector<string>& list_of_files );
double iono_loss( double freq, double theta_deg );
double iono_loss_simple( double freq, double theta_deg ); // 0.01 @ 100 MHz 
double iono_loss_adatta( double freq, double theta_deg );
double refraction_angle( double freq, double theta_deg );
double tropo_refraction_angle( double theta_deg );
double tau_g( double theta_deg );
void test_refraction_angle() ;

void usage()
{
  printf("showspec IN_FILE_LIST -d -s START_UX -i TIME_INTERVAL_SEC -r TIME_STEP -l GEO_LONGITUDE -a GEO_LATITUDE -o OUTPUT_DIR -m -p IN_DIR -f OUT_SPEC_FILENAME -b ANTENNA_PATTERN_FILE -t PATTERN_TYPE -c HORIZON_ALT_CUT -e 20120101.sens -g -j ANTENNA_ROTATION_ANGLE_in_deg -k IONO_ATT0[default 1db] -q PARAM=VALUE\n");
  printf("showspec program integrates sky model with antenna pattern to obtain spectrum. %s\n",SHOWSPEC_VERSION);
  printf("By default show spectrum for one given time\n");
  printf("\n\noptions :\n");
  printf("-d : enable debug\n");
  printf("-e SENSOR_FILE : taking into account antenna efficiency [default BICON table, otherwise option -z]\n");
  printf("-s START_UX : start unix time\n");
  printf("-i TIME_INTERVAL_SEC : time interval\n");
  printf("-r TIME_STEP : time resolution\n");
  printf("-l GEO_LONGITUDE [deg], default %.2f [deg]\n",geo_long);
  printf("-a GEO_LATITUDE  [deg], default %.2f [deg]\n",geo_lat);
  printf("-o OUTPUT_DIR\n");
  printf("-m : save all spectrum maps\n");
  printf("-i IN_DIR : path to directory with input files\n");
  printf("-f OUT_SPEC_FILENAME : output .spec filename [default UT.spec]\n");
  printf("-b ANTENNA_PATTERN_FILE : FEKO output file, by default default dipol pattern is used, (to use biconincal simulation put %s)\n",DEFAULT_ANT_PATTERN_FILE);
  printf("-t PATTERN_TYPE : antenna pattern type 0-dipol over gnd-screen, 1-isotropic [default 0], 2 - Gaussian Beam exp(-(theta/theta0)^2), use -q ant_beam_width=50 to specify Gaussian Beam Width theta0 [default = %.2f]\n",gAntennaPattern.m_AntennaBeamWidth);
  printf("-c HORIZON_ALT_CUT : set minimum elevation above the horizon which counts [default 0 deg]\n");
  printf("-n : do not add contribution from the galaxy\n");  
  printf("-z ant_eff_file.txt : file with antenna efficiency in format : Freq[MHz] Eff [by default hardcoded version for biconical simulation is used]\n");
  printf("-g : 2D interpolation of antenna pattern is used (TGraph2D::Interpolate)\n");
  printf("-j ANTENNA_ROTATION_ANGLE_in_deg : default for Muresk %.2f [deg]\n",CAntPattern::m_RotFromAxis);
  printf("\t\tIt is an angle between Antenna-E direction (PHI=0 in FEKO) and GEOGRAPHICAL-E direction (AZ=270 in libnova / AZ=90 - Randall)\n");
  printf("\t\tFEKO antenna directions         : E PHI=0 , N PHI=90, W PHI=180, S PHI=270\n");
  printf("\t\tASTRONOMY (libnova - used here) : E AZ=270, N AZ=180, W AZ=90  , S AZ=0\n"); 
  printf("-k : include ionosphere effects [disabled by default]\n");
  printf("-p IN_DIR : path where .out files and ni_list file of Angelica sky model are [default %s]\n",gInDir.c_str());
  printf("-q PARAM=VALUE : setting extra paramters like : sun=0 or sun=1 to disable/enable adding solar flux, freq_start, freq_end\n");
  printf("\t\tOther -q parameters : \n");
  printf("\t\t\t\tdebug_level=0             - debug level [default = %d]\n",gDebug);
  printf("\t\t\t\tfreq_start=10            - to specify starting frequency [default = %.2f MHz]\n",gFreqStart);
  printf("\t\t\t\tfreq_end=10              - to specify starting frequency [default = %.2f MHz]\n",gFreqEnd);
  printf("\t\t\t\tsave_map_fits=1          - to save skymaps to fits files\n");
  printf("\t\t\t\tsave_map_at_freq=%.2f    - saves skymap to fits file only for selected frequency [default all]\n",gSaveSkyMapAtFreq);
  printf("\t\t\t\tsave_pattern_map_db=%d   - if save antenna pattern map in db scale (default)\n",gSaveAntPatternDB);
  printf("\t\t\t\tmap_file_base=angelica   - basename of map file [default %s\\n",gMapFileBaseName.c_str());
  printf("\t\t\t\tzenith_basename=%s       - basename of zenith projection fits file [default %s]\n",gZenithProjBaseName.c_str(),gZenithProjBaseName.c_str());
  printf("\t\t\t\tinterpol=1               - enables interpolation of antenna pattern (makes it smoothed)\n");
  printf("\t\t\t\tant_eff=1                - enables usage of antenna efficiency\n");
  printf("\t\t\t\tobject_list=molonglo.txt - load list of objects from TXT file in format : NAME RA_DEG_DEC DEC_DEG_DEC FLUX[Jy] ERR_FLUX[Jy] [optional RADIUS_ARCSEC default %.2f]\n",gSourceRadiusArcSec);
  printf("\t\t\t\tsun=0                    - enables/disabled Sun model adding\n");
  printf("\t\t\t\tgalaxy=0                 - enables/disabled Galaxy model adding\n");
  printf("\t\t\t\tsources=0                - enables/disabled point sources adding (object_list=???)\n");
  printf("\t\t\t\texclude_sources_radius=0 - enables/disabled exclusion of point sources (object_list=???)\n");
  printf("\t\t\t\tsources_radius=0         - sources matching radius (defult = %.2f arcsec)\n",gSourceRadiusArcSec);
  printf("\t\t\t\tsite=muresk              - known sites muresk (DEFAULT), mwa, ebo, ebo201309, ebo201312 (Eyre Bird Observatory), muresk (DEFAULT), wond20140405, wond20140406, wond\n");
  printf("\t\t\t\tt_amb=290                - ambient temperature by default %.2f K is used, or use option -e WEATHER.sens to provide file\n",gT_amb);
  printf("\t\t\t\tnorm_4pi=1               - normalizes (1/4PI) instead of integral from G_tot over 4PI - according to Adrian it should include Ant_eff\n");
  printf("\t\t\t\tfits=1                   - enables saving data to FITS files\n");
  printf("\t\t\t\tcache_on=1               - enabled caching of skymodel files (.out files), by default [CacheON=%d]\n",gSkyMapCache.GetCacheONOFF());
  printf("\t\t\t\tant_cache_on=1           - enabled caching of antenna pattern, by default [CacheON=%d]\n",gAntennaPattern.m_bCacheON);
  printf("\t\t\t\tmax_feko_theta           - maximum value of theta angle simulated in FEKO [default = %.2f deg]\n",CAntPatternParser::m_SimulMaxTheta);
  printf("\t\t\t\tauto_detect_theta        - auto-detect theta angle simulated in FEKO [be default %d , and maxTheta = %.2f deg]\n",CAntPatternParser::m_AutoDetectTheta,CAntPatternParser::m_SimulMaxTheta);
  printf("\t\t\t\tant_beam_width           - antenna beam width, if formula (option -t is used), default value = %.2f [deg]\n",DEFAULT_ANT_BEAM_WIDTH);
  printf("\t\t\t\tsuppress_above_zenangle  - if suppress antenna pattern (formula or FEKO) above given zenithal angle [default = %.2f deg]\n",CAntPattern::gSuppressAntPatternAboveZenAngle);
  printf("\t\t\t\tfast=1                   - fast version, uses much less points (not full 3mln pixel Angelica resolution)\n");
  printf("\t\t\t\tauto_continue=1          - automatically continue if files .spec,.plot are found in local directory ( it will repeat last line in spec/plot files )\n");
  printf("\t\t\t\tbinary_input=1           - read Angelica model from binary .dat files instead of .out files, should be also auto-detected\n");
  printf("\t\t\t\tux_list_file=ux.txt      - list of unix times to generate sky model for, should be according to ux-times of real integratiions, by default not set\n");
  printf("\t\t\t\tant_rotation_deg=%.2f    - rotation of antenna from E-W axis, default is for Muresk %.2f deg (same as option -j)\n",CAntPattern::m_RotFromAxis,CAntPattern::m_RotFromAxis);
  printf("\t\t\t\tinclude_ionosphere=%d    - to include attenuation effects of the ionosphere (see Harish's paper MNRAS 2013), for related parameters (n_e) - see list below :\n",gIncludeIonosphere);
  printf("\t\t\t\tiono_refraction=%d       - include ionospheric refraction in the F-layer (Harish & A.Datta), by default %s\n",gIncludeIonoRefraction,(gIncludeIonoRefraction>0 ? "enabled" : "disabled"));
  printf("\t\t\t\tiono_absorption=%d       - include ionospheric absorption in the D-layer (Harish & A.Datta), by default %s, 1-HARISH's, 2-ABHI DATTA, 3-simple 0.01dBx(100/f)^2xLength(theta_deg)-optional see iono_loss_with_path_len\n",gIncludeIonoAbsorption,(gIncludeIonoAbsorption>0 ? "enabled" : "disabled"));
  printf("\t\t\t\tiono_electemp=%d         - include emission from electrons in the D-layer (Harish & A.Datta) with T_e=%.2f [K], by default %s\n",gIncludeIonoElectronTemp,gElectronTemp,(gIncludeIonoElectronTemp>0 ? "enabled" : "disabled"));
  printf("\t\t\t\telectron_temp=%.2f       - temperature of electrons in the D-layer default = %.2f [K]\n",gElectronTemp,gElectronTemp);
  printf("\t\t\t\t\ttec=%e                 - TEC value (total as measured by GPS satellites), default = %e, in units of 10^16 [e/m^2], so just put 10, 20 or 5 and no 10e16\n",gTEC_total,gTEC_total);
  printf("\t\t\t\t\tn_e=2.5e8                - concentration [1/m^3] of electrons\n");
  printf("\t\t\t\t\tinterpol=%d            - enable ROOT based interpolation of antenna pattern (default = %d)\n",CAntPatternParser::m_bEnableInterpolation,CAntPatternParser::m_bEnableInterpolation);
  printf("\t\t\t\t\tzenith_proj_radius=%d  - zenith projection FITS file size (2*r+1) [default = %d] -> fits file size = %d x %d\n",gZenithProjRadius,gZenithProjRadius,(2*gZenithProjRadius+1),(2*gZenithProjRadius+1));
  printf("\t\t\t\t\ttest_tsky_with_tau=%d  - test integration of T_sky(theta,phi)*P(theta,phi)*tau(theta) for the 2015 ionospheric paper, by default DISABLED [%d]\n",gTestTskyWithTau,gTestTskyWithTau);
  printf("\t\t\t\t\tiono_loss_with_path_len=%d - ionospheric loss scaled by length of path in the D-layer [H_D=%.2f km, dh_D=%.2f km, Re=%.2f km] , by default %d\n",gIonoLossWithPathLength,gHeightD,gDeltaHeightD,gRadiusEarth,gIonoLossWithPathLength);
  printf("\t\t\t\t\ttypical_loss_at_100mhz=%.2f - typical loss at 100 MHz for simple absorption [default %.2f dB]\n",gTypicalLossAt100MHz,gTypicalLossAt100MHz);
  printf("\t\t\t\t\ttypical_refraction_angle_arcmin=%.2f - typical refraction angle at  100 MHz for simple refraction model [default %.2f arcmin]\n",gTypicalRefractionAngleArcMin,gTypicalRefractionAngleArcMin);
  printf("\t\t\t\t\ttypical_refraction_zenith_angle_deg=%.2f - zenithal angle to provide typical refraction scale at  100 MHz for simple refraction model [default %.2f arcmin]\n",gTypicalRefractionZenithAngle,gTypicalRefractionZenithAngle);
  printf("\t\tTROPOSHERIC PARAMETERS:\n");
  printf("\t\t\t\tinclude_tropo_refr=%d    - whether enable tropospheric refraction [default %d]\n",gIncludeTropoRefraction,gIncludeTropoRefraction);
  printf("\t\t\t\tatm_pressure=%.2f        - atmospheric pressure in millibars = hPa [default %.2f millibars]\n",gAtmPressure,gAtmPressure);
  
#ifdef _SUPERFAST_OPTIMIZED_
  printf("\n\n");
  printf("WARNING : some options are not available in version compiled with -D_SUPERFAST_OPTIMIZED_ - check makefile to enable them, disabled options : galaxy=0 , -k (ionospheric effects)\n");
  printf("Unavailable options :\n");
  printf("\t\tsave_map_at_freq=? - saving sky image to fits file\n");
  printf("Please, recompile the program with -D_SUPERFAST_OPTIMIZED_ commented out in the Makefile\n");
  printf("\n\n");
#endif

  exit(0);
}

void parse_cmdline(int argc, char *argv[]) {
    char optstring[] = "gwmndhs:i:r:l:a:o:p:f:b:t:x:c:e:z:j:k:q:";
    int opt;

    while ((opt = getopt(argc, argv, optstring)) != -1) {
        switch (opt) {
            case 'n':
               gAddGalaxy=0;
               break;
            case 'g':
               CAntPatternParser::m_bEnableInterpolation = 1;
#ifdef _SUPERFAST_OPTIMIZED_
                            printf("ERROR : antenna pattern interpolation could not be enabled, due to pre-compiler directive (recompile program without -d _SUPERFAST_OPTIMIZED_ or disable by : -q interpol=0 )\n");
                            exit(-1);
#endif
#ifdef _NO_ROOT_
                            printf("ERROR : antenna pattern interpolation could not be enabled, due to pre-compiler directive (recompile program without -d _NO_ROOT_ or disable by : -q interpol=0 )\n");
                            exit(-1);
#endif                         
               
               break;
            case 'b':
                if( optarg ){
                   CAntPattern::m_AntennaPatternFile = optarg;
                }
                break;

            case 'c':
                if( optarg ){
                   gOnlyAboveHorizont = atof( optarg );
                }
                break;

            case 'e':
                if( optarg ){
                   gWeather.m_szWeatherFile = optarg;
                }
                break;

            case 't':
                if( optarg ){
                   gAntennaPattern.m_AntennaPatternType = (eAntennaPatternType)(atol(optarg));
                   CAntPattern::gAntennaPatternFormula = 1;
                }
                break;
            case 'd':
                gDebug++;
                break;
            case 'h':
                usage();
                break;
            case 'k':
                gIncludeIonosphere = 1;
                #ifdef _SUPERFAST_OPTIMIZED_
                   printf("WARNING : the program compilation options (_SUPERFAST_OPTIMIZED_) require optimisation -> ionospheric effects are disabled !!!\n");
                #endif
                
                if( optarg ){
                   gIonoAtt0 = atof(optarg);
                }
                break;
            case 's':
                start_ux_time = atof(optarg);
                break;
            case 'i':
                time_interval = atoi(optarg);
                break;
            case 'r':
                dt = atof(optarg);
                break;
            case 'l':
                geo_long = atof(optarg);
                break;
            case 'a':
                geo_lat = atof(optarg);
                break;
            case 'o':
                gOutPutDir = optarg;
                break;
            case 'f':
                gOutSpecFileName = optarg;
                break;
            case 'p':
                gInDir = optarg;
                break;
            case 'm':
                gDoSaveAllMaps=1;
                break;
            case 'x':
                if( optarg ){
                   gThetaMax = atof(optarg);
                }
                break;
            case 'w':
                gAddAntEffCorrection = 1;
                break;
            case 'z':
                if( optarg ){
                   gAntEffFile = optarg;
                }
                break;

            case 'j':
                if( optarg ){
                   CAntPattern::m_RotFromAxis = atof(optarg);
                }
                break;

            case 'q':
                if( optarg ){
                   MyParser pars = optarg;
                   CEnvVar envvar;
                   if( pars.GetVarAndValue(envvar) ){
                      if( strcmp(envvar.szName.c_str(),"debug_level")==0 ){
                         gDebug = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"sun")==0 ){
                         gAddSun = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"galaxy")==0 ){
                         gAddGalaxy = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"sources")==0 ){
                         gAddSources = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"exclude_sources_radius")==0 ){
                         gExcludeSourcesRadiusArcSec = atof(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"sources_radius")==0 ){
                         gSourceRadiusArcSec = atof(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"freq_start")==0 ){
                         gFreqStart = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"freq_end")==0 ){
                         gFreqEnd = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"save_map_fits")==0 ){
                         gSaveSkyMapFits = atol(envvar.szValue.c_str());
#ifdef _SUPERFAST_OPTIMIZED_
                         printf("\n\n\nWARNING : this option is not available when program is compiled with _SUPERFAST_OPTIMIZED_ directive !!!\n");
#endif                         
                      }
                      if( strcmp(envvar.szName.c_str(),"save_map_at_freq")==0 ){
                         gSaveSkyMapAtFreq = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"save_pattern_map_db")==0 ){
                         gSaveAntPatternDB = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"zenith_basename")==0 ){
                         gZenithProjBaseName = envvar.szValue.c_str();
                      }
                      if( strcmp(envvar.szName.c_str(),"interpol")==0 ){
                         CAntPatternParser::m_bEnableInterpolation = atol(envvar.szValue.c_str());
                         if( CAntPatternParser::m_bEnableInterpolation > 0 ){
#ifdef _SUPERFAST_OPTIMIZED_
                            printf("ERROR : antenna pattern interpolation could not be enabled, due to pre-compiler directive (recompile program without -d _SUPERFAST_OPTIMIZED_ or disable by : -q interpol=0 )\n");
                            exit(-1);
#endif
#ifdef _NO_ROOT_
                            printf("ERROR : antenna pattern interpolation could not be enabled, due to pre-compiler directive (recompile program without -d _NO_ROOT_ or disable by : -q interpol=0 )\n");
                            exit(-1);
#endif                         
                         }
                      }
                      if( strcmp(envvar.szName.c_str(),"ant_eff")==0 ){
                         gAddAntEffCorrection = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"object_list")==0 ){
                         gObjectListFile = envvar.szValue.c_str();
                      }
                      if( strcmp(envvar.szName.c_str(),"map_file_base")==0 ){
                         gMapFileBaseName = envvar.szValue.c_str();
                      }
                      if( strcmp(envvar.szName.c_str(),"t_amb")==0 ){
                         gT_amb = atof( envvar.szValue.c_str() );
                      }
                      if( strcmp(envvar.szName.c_str(),"norm_4pi")==0 ){
                         gNorm4PI=1;
                      }
                      if( strcmp(envvar.szName.c_str(),"fits")==0 ){
                         gSaveFitsFile = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"max_feko_theta")==0 ){
                         CAntPatternParser::m_SimulMaxTheta = atof(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"auto_detect_theta")==0 ){
                         CAntPatternParser::m_AutoDetectTheta = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"ant_beam_width")==0 ){
                         gAntennaPattern.m_AntennaBeamWidth = atof(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"suppress_above_zenangle")==0 ){
                         CAntPattern::gSuppressAntPatternAboveZenAngle = atof(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"ant_rotation_deg")==0 ){
                         CAntPattern::m_RotFromAxis = atof(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"cache_on")==0 ){
                         gSkyMapCache.SetCache( atol(envvar.szValue.c_str()) );
                      }
                      if( strcmp(envvar.szName.c_str(),"ant_cache_on")==0 ){
                         gAntennaPattern.m_bCacheON = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"fast")==0 ){
                         gRunFastVersion = atol(envvar.szValue.c_str());
                         CAntPattern::m_bCacheON = 1;
                      }
                      if( strcmp(envvar.szName.c_str(),"auto_continue")==0 ){
                         gAutoContinue = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"binary_input")==0 ){
                         CSkyMapCache::m_bBinaryFile = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"ux_list_file")==0 ){
                         gUxListFile = envvar.szValue.c_str();
                      }
                      if( strcmp(envvar.szName.c_str(),"n_e")==0 ){
                         gN_e = atof( envvar.szValue.c_str() );
                         gPlasmaFreq = 9.00 * sqrt(gN_e);
                      }
                      if( strcmp(envvar.szName.c_str(),"include_ionosphere")==0 ){                         
                         gIncludeIonoAbsorption = ((eIonoAbsModelType)atol(envvar.szValue.c_str()));
                         if( gIncludeIonoAbsorption > 0 ){
                            gIncludeIonosphere = 1;
                            gIncludeIonoRefraction = eIonoRefractionSimpleWithThetaDep;
                         }else{
                            gIncludeIonosphere = 0;
                            gIncludeIonoRefraction  = eNoIonoRefraction;
                         }
                      }
                      if( strcmp(envvar.szName.c_str(),"iono_refraction")==0 ){
                         gIncludeIonoRefraction = (eIonoRefractionModelType)atol(envvar.szValue.c_str());
                         if( gIncludeIonoRefraction > 0 ){
                            if( atol(envvar.szValue.c_str())==1 || atol(envvar.szValue.c_str())==2 ){
                               gIncludeIonoRefraction = eIonoRefractionHarishAbhiData;
                            }else{
                               gIncludeIonoRefraction = ((eIonoRefractionModelType)atol(envvar.szValue.c_str()));
                            }
                            gIncludeIonosphere = 1;
                         }
                      }
                      if( strcmp(envvar.szName.c_str(),"iono_absorption")==0 ){
                         gIncludeIonoAbsorption = ((eIonoAbsModelType)atol(envvar.szValue.c_str()));
                         if( gIncludeIonoAbsorption > 0 ){
                            gIncludeIonosphere = 1;
                         }
                      }
                      if( strcmp(envvar.szName.c_str(),"iono_loss_with_path_len")==0 ){
                         gIonoLossWithPathLength = atol(envvar.szValue.c_str());                         
                      }
                      if( strcmp(envvar.szName.c_str(),"iono_electemp")==0 ){
                         gIncludeIonoElectronTemp = atol(envvar.szValue.c_str());
                         // if( gIncludeIonoAbsorption > 0 ){
                         //    gIncludeIonosphere = 1;
                         // }
                      }
                      if( strcmp(envvar.szName.c_str(),"electron_temp")==0 ){
                         gElectronTemp = atof(envvar.szValue.c_str());
                         // if( gIncludeIonoAbsorption > 0 ){
                         //    gIncludeIonosphere = 1;
                         // }
                      }
                      if( strcmp(envvar.szName.c_str(),"typical_loss_at_100mhz")==0 ){
                         gTypicalLossAt100MHz = atof(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"typical_refraction_angle_arcmin")==0 ){
                         gTypicalRefractionAngleArcMin = atof(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"typical_refraction_zenith_angle_deg")==0 ){
                         gTypicalRefractionZenithAngle = atof(envvar.szValue.c_str());
                      }

                      if( strcmp(envvar.szName.c_str(),"test_tsky_with_tau")==0 ){
                         gTestTskyWithTau = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"tec")==0 ){
                         gTEC_total = atof(envvar.szValue.c_str())*1e16;
                      }
                      if( strcmp(envvar.szName.c_str(),"zenith_proj_radius")==0 ){
                         gZenithProjRadius = atol(envvar.szValue.c_str());
                      }
                      // TROPOSPHERE :
                      if( strcmp(envvar.szName.c_str(),"include_tropo_refr")==0 ){
                         gIncludeTropoRefraction = atol(envvar.szValue.c_str());
                      }
                      if( strcmp(envvar.szName.c_str(),"atm_pressure")==0 ){
                         gAtmPressure = atol(envvar.szValue.c_str());
                      }
                        
                      if( strcmp(envvar.szName.c_str(),"site")==0 ){
                         if( strcmp(envvar.szValue.c_str(),"mwa") == 0 || strcmp(envvar.szValue.c_str(),"mro") == 0 ){
                            geo_long = 116.670815;
                            geo_lat  = -26.703319;
                            gSiteName = "MWA";
                         }
                         if( strcmp(envvar.szValue.c_str(),"muresk") == 0 || strcmp(envvar.szValue.c_str(),"MURESK") == 0 ){
                            geo_long=116.68609722;
                            geo_lat=-31.74625000;
                            gSiteName = "MURESK";
                         }
                         if( strcmp(envvar.szValue.c_str(),"ebo") == 0 || strcmp(envvar.szValue.c_str(),"EBO") == 0 || strcasecmp(envvar.szValue.c_str(),"ebo201309") == 0 ){
                            // 2013-09
                            geo_long = 126.30268889; // was 126.3013;
                            geo_lat  = -32.25256389; // was -32.246391;
                            gSiteName = "EBO201309";
                         }
                         if( strcasecmp(envvar.szValue.c_str(),"ebo201312") == 0  ){
                            // 2013-12 :
                            geo_long = 126.29777778; // was 126.3013;
                            geo_lat  = -32.25250000; // was -32.246391;
                            gSiteName = "EBO201312";
                         }
                         if( strcasecmp(envvar.szValue.c_str(),"wond20140406") == 0  || strcasecmp(envvar.szValue.c_str(),"wond")==0 ){
                            // 2014-04-06 :
                            geo_long = 118.43999167; // - 27deg 51'10.31''
                            geo_lat  = -27.85286389; // 118deg 26' 23.97''
                            gSiteName = "WONDINONG_20140406";
                         }
                         if( strcasecmp(envvar.szValue.c_str(),"wond20140405") == 0  ){
                            // 2014-04-05 :
                            // TODO: change according to Randall's info !
                            geo_long = 118.43999167; // - 27deg 51'10.31''
                            geo_lat  = -27.85286389; // 118deg 26' 23.97''
                            gSiteName = "WONDINONG_20140406";
                         }
                      }
                   }
                }
                break;


            default:
                fprintf(stderr,"Unknown option %c\n",opt);
                usage();
        }
    }
    
    gInDir.env2str();
    
    if( strlen(gAntEffFile.c_str()) > 0 ){
       printf("Reading efficiency file %s ...",gAntEffFile.c_str());
       int n_eff_count = read_file( gAntEffFile.c_str() , gAntEffTable );
       printf("Read %d points\n",n_eff_count);
    }
}


cObjectInfo* FindPointSource( double ra, double dec, double radius_arcsec_param )
{
   for(int i=0;i<gObjectList.size();i++){
      cObjectInfo& object = gObjectList[i];
      double radius_arcsec = radius_arcsec_param;      
      if( object.radius_arcsec > 0 ){
         radius_arcsec = object.radius_arcsec;
      }
      
      if( get_ang_distance_deg( object.ra, object.dec, ra, dec ) < (radius_arcsec/3600.00) ){
         return &object;
      }
   }
   
   return NULL;
}


int IsPointSource( double ra, double dec, double radius_arcsec_param )
{
   for(int i=0;i<gObjectList.size();i++){
      cObjectInfo& object = gObjectList[i];
      double radius_arcsec = radius_arcsec_param;      
      if( object.radius_arcsec > 0 ){
         radius_arcsec = object.radius_arcsec;
      }
      
      if( get_ang_distance_deg( object.ra, object.dec, ra, dec ) < (radius_arcsec/3600.00) ){
         return 1;
      }
   }
   
   return 0;
}

int read_sources( const char* infile )
{
   MyIFile ifile( infile );
   gObjectList.clear();

   const char* pLine=NULL;
   while( pLine = ifile.GetLine( TRUE ) ){
      if( mystring::get_first_non_white( pLine )=='#' )
         continue; 
                                      
      MyParser pars=pLine;
      CMyStrTable items;  
      pars.GetItems( items );                                                        

      cObjectInfo source;
      source.radius_arcsec = -1;
      source.period = 0.00;
      source.pulse_length = 0.00;
      
      source.name = items[0].c_str();
      source.ra = atof(items[1].c_str())*15.00; // hours -> deg 
      source.dec = atof(items[2].c_str());      
      source.flux = atof(items[3].c_str());
      source.err_flux = atof(items[4].c_str());

      if( items.size() >= 6 ){
         source.radius_arcsec = atof(items[5].c_str());
      }else{
         source.radius_arcsec = -1;
      }

      if( items.size() >= 7 ){
         source.period = atof(items[6].c_str());
      }
      if( items.size() >= 8 ){
         source.pulse_length = atof(items[7].c_str());
      }
      
      gObjectList.push_back(source);
   }
   
   printf("Read %d sources from file %s\n",gObjectList.size(),infile);

   return gObjectList.size();
}

inline double get_sun_brightness( double ra, double dec, double azim, double alt, // coordinates of integrated sky point 
                                  double sun_ra, double sun_dec                   // Sun coordinates for current moment 
                                )
{
   if( alt > 0 ){
//      double dist_to_sun = get_ang_distance_deg( ra, dec, sun_ra, sun_dec );   

      // local calculation of Sun distance to make it faster :
      double cos_x = sin(dec*deg2rad_const)*sin(sun_dec*deg2rad_const) + cos(dec*deg2rad_const)*cos(sun_dec*deg2rad_const)*cos( (ra-sun_ra)*deg2rad_const );
      double x = acos( cos_x );
      double dist_to_sun = x * rad2deg_const;
            
      
      if( dist_to_sun < sun_radius_deg ){
         return sun_T_b;
      }
   }
                                                
   return 0.00;                                                      
}


double get_sun_brightness( double ra, double dec, time_t unix_time )
{
   double sun_ra, sun_dec, sun_az, sun_alt;
   get_sun_all( unix_time, geo_long, geo_lat, sun_ra, sun_dec, sun_az, sun_alt );

   double dist_to_sun = get_ang_distance_deg( ra, dec, sun_ra, sun_dec );
   
   if( dist_to_sun < sun_radius_deg ){
      double part_az, part_alt;
      radec2azh( ra, dec, unix_time, geo_long, geo_lat, part_az, part_alt );
      
      if( part_alt > 0 ){
         return sun_T_b;
      }
   }
   
   return 0.00;
}

arr<double>* fill_healpix_array( vector<cSkyInfo>& sky_intensities , Healpix_Map<double>& healpix_skymap  )
{
   int sky_intensities_count = sky_intensities.size();
   arr<double>* healpix_arr = new arr<double>(sky_intensities_count);

  // fill healpix array:
  int i=0;
  for(vector<cSkyInfo>::iterator it=sky_intensities.begin();it!=sky_intensities.end();it++){
     cSkyInfo& skyinfo = (*it);         
     (*healpix_arr)[i] = skyinfo.intensity;
     i++;
  }
  

//     Healpix_Map<double> healpix_skymap(norder,RING);
  healpix_skymap.Set(*healpix_arr,RING);

  printf("fill_healpix_array - filled healpix array with %d values\n",sky_intensities_count);  
  return healpix_arr;
}


double integrate_galactic_spectrum( const char* infile , double unix_time, ofstream* ofile, double& out_pattern_sum, double& out_sun_sum, int first_freq  )
{
  // uses or not cache - depending on cache_on=1 or =0 parameter (by default disabled)
  vector<cSkyInfo>* pSkyIntensities = gSkyMapCache.get_skyintensities( infile );  
  if( !pSkyIntensities ){
     printf("ERROR : could not get skyintensities from in-file = %s\n",infile);
     exit(-1);     
  }
  vector<cSkyInfo>& sky_intensities = (*pSkyIntensities);

  Healpix_Map<double> healpix_skymap(norder,RING);
  arr<double>* sky_arr = fill_healpix_array( sky_intensities , healpix_skymap );

//  vector<cSkyInfo>& sky_intensities = CSkyMapCache::m_SkyIntensities;
//  CSkyMapCache::read_infile( infile );

  int npix=0;
  string line;

  MollweideSkyMap map2(100);   
 
  double JD_time = uxd2jd( unix_time );
  double freq_mhz = parse_ni( infile );
//  double ant_gain_num = get_ce300e_gain(freq_mhz);
  double ant_gain_num = 1.00;

// MODELLING 2012-06-27
//  double gm_filter = get_gm_transmiss(freq_mhz);
  double gm_filter = 1.00;

  double t_amb = gT_amb;
  if( gWeather.GetCount() > 0 ){
     t_amb = gWeather.GetAmbTempK(unix_time);
  }
       
  printf("First freq flag = %d\n",first_freq);
  printf("Npix=%d -> nside=%d -> norder=%d -> npix=%d , dOMEGA=%.10Lf [srad] (OK?)\n",(12*512*512),(int)nside,(int)norder,(int)pHealpix->Npix(),dOMEGA);
  printf("Freq = %.2f [MHz], ant_gain = %.4f, gm_filter = %.4f , t_amb = %.2f [K]\n",freq_mhz,ant_gain_num,gm_filter,t_amb);
        
  long double sum=0.0000;  
  char outline[1024],outline_short[256];
  long double antenna_pattern_sum = 0.00;
  double antenna_pattern_sum_below_hor = 0.00;
  long double sum_below_hor=0.00;
  long double sun_sum=0.0000;
  double sum_tau_g = 0.00, sum_tau_g_bar=0.00;

/*  int freq_index = -1;
  CFreqPattern* pPatternFreq = NULL;
  CFreqPattern* pPatternFreqNext = NULL;
  if( gAntennaPattern.GetAntParsedFile() ){
     freq_index = (gAntennaPattern.GetAntParsedFile())->FindFreq(freq_mhz);
     pPatternFreq = (gAntennaPattern.GetAntParsedFile())->GetFreqPattern( freq_index );
     pPatternFreqNext = (gAntennaPattern.GetAntParsedFile())->GetFreqPattern( freq_index+1 );

     printf("integrate_galactic_spectrum : Freq index for %.2f MHz = %d -> pointers are 0x%x and 0x%x\n",freq_mhz,freq_index,pPatternFreq,pPatternFreqNext);
     
     if( !pPatternFreq || !pPatternFreqNext ){
        printf("ERROR in code : could not find antenna pattern for frequency index %d or %d (freq_mhz=%.2f MHz)\n",freq_index,(freq_index+1),freq_mhz);
        exit(-1);
     }
  }*/

  int bAddSun = gAddSun;
  double sun_ra, sun_dec, sun_az, sun_alt;
  if( bAddSun ){
     // check if above horizon :
     get_sun_all( unix_time, geo_long, geo_lat, sun_ra, sun_dec, sun_az, sun_alt );
     
     if( sun_alt <= -1.00 ){
         // ignore Sun if deep below horizon :
         bAddSun = 0;
     }   
  }

  double min_iono_att_db = 1000000.00;
  double max_iono_att_db = -1000000.00;

  int sky_intensities_count = sky_intensities.size();
#ifndef _SUPERFAST_OPTIMIZED_
  arr<double>* healpix_arr = NULL;
  arr<double>* healpix_antpat_arr = NULL;
  arr<double>* healpix_deltaTa_arr = NULL;  
  CBgFits* p_zenith_projection = NULL;
  CBgFits* p_zenith_projection_ant = NULL;
  CBgFits* p_zenith_projection_mult = NULL;

  if( gSaveSkyMapFits > 0 && (gSaveSkyMapAtFreq<0 || fabs(gSaveSkyMapAtFreq-freq_mhz)<1 ) ){
     healpix_arr = new arr<double>(sky_intensities_count);
     healpix_antpat_arr = new arr<double>(sky_intensities_count);
     healpix_deltaTa_arr = new arr<double>(sky_intensities_count);
          
     p_zenith_projection = new CBgFits(2*gZenithProjRadius+1,2*gZenithProjRadius+1);
     p_zenith_projection_ant = new CBgFits(2*gZenithProjRadius+1,2*gZenithProjRadius+1);
     p_zenith_projection_mult = new CBgFits(2*gZenithProjRadius+1,2*gZenithProjRadius+1);
  }
#endif

  if( first_freq > 0 ){
     printf("INFO : first frequency for given time (JD_time=%.8f) -> initializing (AZIM,ALT) data in sky_intensities array ...",JD_time);fflush(stdout);
     time_t start = get_dttm();
     for(vector<cSkyInfo>::iterator it=sky_intensities.begin();it!=sky_intensities.end();it++){
        cSkyInfo& skyinfo = (*it);
        
       // at first frequency for given UNIX_TIME fill coorindates info as sky changes in time :
       // using libnova - means in convention of AZIM : S=0, W=90, N=180, E=270 deg 
       //       gal2hor( skyinfo.glon, skyinfo.glat, geo_long, geo_lat, unix_time, skyinfo.azim, skyinfo.alt, skyinfo.ra, skyinfo.dec );
       gal2hor_fromJD( skyinfo.glon, skyinfo.glat, geo_long, geo_lat, JD_time, skyinfo.azim, skyinfo.alt, skyinfo.ra, skyinfo.dec );
       
/* I am testing alternative way to add point sources in the very end as a separate loop

#ifndef _SUPERFAST_OPTIMIZED_
       if( gObjectList.size() > 0  ){
          double radius_arcsec = gSourceRadiusArcSec;
          if( gExcludeSourcesRadiusArcSec > 0 ){
             radius_arcsec = gExcludeSourcesRadiusArcSec;
          }
          skyinfo.pPointSource = FindPointSource( skyinfo.ra, skyinfo.dec, radius_arcsec );
       }
#endif */
     }
     printf("took %d sec\n",(get_dttm()-start));fflush(stdout);
  } 
  
  // get antenna pattern array (map) :
  time_t start = get_dttm();
  double* antenna_pattern_array = NULL;  
  if( gIncludeIonoRefraction <= 0 ){
     printf("Getting antenna pattern map ...");fflush(stdout);
     antenna_pattern_array = gAntennaPattern.antenna_pattern( infile,  freq_mhz, sky_intensities );
     printf("took %d sec\n",(get_dttm()-start));fflush(stdout);
  }else{
     printf("WARNING : ionospheric refraction is turned on (-q refraction=1) - antenna pattern is calculated for every single point later in the code\n");
  }
  
  if( gDebug >= 3 && antenna_pattern_array ){
     printf("Antenna gain at %.2f MHz debug :",freq_mhz);
     for(int ll=0;ll<sky_intensities.size();ll+=10000){
        printf(" %.20f",antenna_pattern_array[ll]);       
     }
     printf("\n");
  }

  start = get_dttm();
  printf("Starting loop over %d elements ...",sky_intensities_count);fflush(stdout);
//  for(int k=0;k<sky_intensities_count;k++){ // OPTIMIZE : use iterator !
  for(vector<cSkyInfo>::iterator skyinfo=sky_intensities.begin();skyinfo!=sky_intensities.end();skyinfo++){
//    cSkyInfo& skyinfo = (*it); 
//    double intensity=skyinfo->intensity;

    double phi_deg = skyinfo->coord.phi * rad2deg_const;     // OPTIMIZE value times constant rad2deg 
    double theta_deg = skyinfo->coord.theta * rad2deg_const; // OPTIMIZE value times constant rad2deg 
    
    
    pointing ptg( skyinfo->coord.theta , skyinfo->coord.phi );
    double intensity_test = healpix_skymap.interpolated_value( ptg );           
    if( fabs(intensity_test-skyinfo->intensity) > 0.0001 ){
       printf("ERROR at (phi,theta) = (%.4f , %.4f ) [deg] : interpolated value different from standrd by more than 0.0001 (%.8f != %.8f)\n",phi_deg,theta_deg,skyinfo->intensity,intensity_test);
    }    
    if( gDebug > 0 ){
       printf("TESTING intensity (%.4f , %.4f ) [deg] = %.8f vs %.8f\n",phi_deg,theta_deg,skyinfo->intensity,intensity_test);
    }
    skyinfo->intensity = intensity_test; // overwriting with test value

//    double glon = skyinfo->glon;
//    double glat = skyinfo->glat;
//    double xp = skyinfo->xp;
//    double yp = skyinfo->yp;

    if( theta_deg >= 180.00 || npix>=3145728 ){	
       printf("WARNING : last line skipped !\n");
       phi_deg = 0.00;
       break;
    }

    // calculate horizontal coordinates from galactic, assuming that antena is beam is pointing to ZENITH : 
//    double azim,alt,ra,dec,zen_dist;    
//    azim = skyinfo->azim;
//    alt  = skyinfo->alt;
//    ra   = skyinfo->ra;
//    dec  = skyinfo->dec;
    double zen_dist = (90.00-skyinfo->alt);

#ifndef _SUPERFAST_OPTIMIZED_
    if( fabs(skyinfo->glon)<0.1 && fabs(skyinfo->glat)<0.1 ){
       printf("DEBUG : GAlaxy center (GLON,GLAT)=(%.2f,%.2f) -> (RA,DEC)=(%.2f,%.2f) -> (AZIM,ALT)=(%.2f,%.2f)\n",skyinfo->glon, skyinfo->glat, skyinfo->ra, skyinfo->dec, skyinfo->azim, skyinfo->alt);
    }
#endif

    double delta_sky = skyinfo->intensity;    
    if( gOnlyAboveHorizont >= 0 ){
       if( skyinfo->alt < gOnlyAboveHorizont /* && alt>0 TEST with signal from underground added - according to antenna-over-gnd-screen-pattern */ ){
          // skip below horizont if required 
//          delta_sky = 0.00;
          delta_sky = t_amb; // 290K below minium elevation 
       }
       
       if( CAntPatternParser::m_SimulMaxTheta < 180 ){
         if( zen_dist > CAntPatternParser::m_SimulMaxTheta ){
            delta_sky = 0;
         }
       }
    }
    
#ifndef _SUPERFAST_OPTIMIZED_
    if( gAddGalaxy <= 0 ){
       int bAdd=0;
       if( gAddSources ){
          if( skyinfo->pPointSource ){
             // if it is point source and point sources are required set add flag :
             bAdd = 1;     
          }
       }   
       
       if( !bAdd ){   
          delta_sky = 0.00;     
       }
    }else{
       if( gExcludeSourcesRadiusArcSec > 0 ){
          // excluding point sources if needed :
          if( skyinfo->pPointSource ){
             delta_sky = 0.00;
          }
       }
    }

    if( gDebug > 0 ){
       if( (npix % 1000) == 0 ){
          printf("npix = %d / %d = %.2f %%...\n",npix,sky_intensities.size(),(double(npix)/double(sky_intensities.size()))*100.00);fflush(stdout);
       }
    }
#endif

    // OPTIMIZED Solar code (should be fast and OK) :
    double T_sun = 0.00;
    if( bAddSun >0 && skyinfo->alt>0 ){ // if Sun is enabled in the code and pixel Alt>0 :
       double fabs_dec = fabs_inline( skyinfo->dec );
       if( fabs_dec < 26.00 ){ // if close enough from the celestial Equator 
          // (RA1-RA2)=10deg -> angular distance = 8.90773 deg at DEC=27 deg 
          // (RA1-RA2)= 5deg -> angular distance = 4.45474 deg at DEC=27 deg          
          // this preliminary check is to avoid calculation of angular distance using sin/cos (CPU expensive) calls :
          double ra_diff = fabs_inline(skyinfo->ra - sun_ra);
          double dec_diff = fabs_inline(skyinfo->dec - sun_dec);
          
          if( ra_diff<5.00 && dec_diff<5.00 ){
             // only in case we are close to Sun call get_sun_brightness (expensive sin/cos calculations)
             T_sun = get_sun_brightness( skyinfo->ra, skyinfo->dec, skyinfo->azim, skyinfo->alt, sun_ra, sun_dec );
             delta_sky += T_sun;
             // sun_sum += T_sun;
          }
       }
    }
    
    delta_sky += gExternalNoiseT;

    // NEW 2013-09-17 - not needed. As Randall explained to me, the sky can have lower temperature than ambient. At GHz-s it might be even ~20K that's why 
    // VLA is a Cassagrain not to have spil-over from the sky, but GMRT is not and can pick-up ground temperature 
    // in our case we can have T_ant < T_amb (at higher frequencies, but there is this ripple which probably spoils my picture)
    // if <T_ambient -> set T_ambient :
//    if( delta_sky < gT_amb ){
//       delta_sky = gT_amb;
//    }

   // INTEGRAL ~= SUM :
   double pattern = 0;
   if( antenna_pattern_array && CAntPattern::gAntennaPatternFormula<=0 ){
      pattern = antenna_pattern_array[npix];       
   }else{
      // 2015-03-17 seems that it was bug phi_deg = glon !!! not azim !!! or whatever should be here 
      // a OK but for the Gaussian Beam (symetric it does not matter ) 
      pattern = gAntennaPattern.antenna_pattern_formula( freq_mhz, phi_deg, zen_dist );
   }
   double pattern_for_pattern_sum = pattern;
   

#ifndef _SUPERFAST_OPTIMIZED_
    if( gThetaMax > 0 ){
       if( theta_deg > gThetaMax ){
          delta_sky = 0.00;
       }
    }

    if( gIncludeIonosphere ){
/*       // Formula from : https://www.cv.nrao.edu/~demerson/ionosphere/atten/atten.htm
       // it is valid for "typical" nighttime, but it is said that that in nighttime attenuations are typically a factor of 10-30 lower 
       // it means sky model should match better for nighttimes and worse for day time, I should double check it again - but look 
       // at low frequencies 
       double h_km = 100.00;
       double R_earth = 6367.00;
       double takeoff_angle = skyinfo->alt*(M_PI/180.00);
//       double att_dt = 2200.00 / (freq_mhz*sqrt(sin(takeoff_angle)*sin(takeoff_angle) + 2.00*(h_km/R_earth)));
//       double att_num = exp(log(10)*(att_dt/10.0));
       
       // Several models for :
       // 2200/f_mhz^2 ( https://www.cv.nrao.edu/~demerson/ionosphere/atten/atten.htm )
       // I assume freq_1db = sqrt(2200) :
       double freq_1db = sqrt(2200.00);
       // EDGES assumation freq_1db = 18 MHz : 
//       double freq_1db = 18.00;
       // John Kennewell presentation from Global Eor Workshop ( http://www.caastro.org/files/30/1357872276/terkildsen_raiono.pdf ) 
       // but this looks rather like order of magnitude 
       // double freq_1db = 10.00;
       double freq_scaled = (freq_mhz/freq_1db);
              
       double att_db = gIonoAtt0 / (freq_scaled*freq_scaled*sqrt(sin(takeoff_angle)*sin(takeoff_angle) + 2.00*(h_km/R_earth)));
       double att_num = exp(log(10)*(att_db/10.0));  

       if( att_db < min_iono_att_db ){
          min_iono_att_db = att_db;
       }
       if( att_db > max_iono_att_db ){
          max_iono_att_db = att_db;
       }
          
       delta_sky = delta_sky / att_num;*/
    
       if( gIncludeIonoAbsorption > 0 ){
          double att_num = 1.00;
          
          switch( gIncludeIonoAbsorption ){
             case eNoIonoAbhiDatta :
                // as in Abhi Datta paper 2014 : 
                att_num = iono_loss_adatta( freq_mhz, zen_dist );
                break;
             
             case eNoIonoHarish :
                // As In Harish's paper 2013 MNRAS :
                att_num = iono_loss( freq_mhz, zen_dist );
                break;
                
             case eIonoLossSimple :
                att_num = iono_loss_simple( freq_mhz, zen_dist );
                break;
             
             defaul :
                att_num = iono_loss_simple( freq_mhz, zen_dist );
                break;
          }
          
 //       double att_num = exp(log(10)*(att_db/10.0));
          if( gDebug > 10 ){
             printf("DEBUG : iono-loss = %.4f -> delta_sky := %.2f * %.2f + [%.2f * %.2f]*%d\n",att_num,att_num,delta_sky,(1-att_num),gElectronTemp,gIncludeIonoElectronTemp);
          }                                 

          delta_sky = delta_sky * att_num;          
          if( gIncludeIonoElectronTemp > 0 ){
             delta_sky = delta_sky + (1-att_num)*gElectronTemp;
          }
       }
       
       if( gIncludeIonoRefraction > 0 ){
          double delta_theta_deg = refraction_angle( freq_mhz, zen_dist );
          double delta_theta_arcsec = delta_theta_deg*3600.00;

          zen_dist = zen_dist - delta_theta_deg;
          if( zen_dist<0 ){
             zen_dist = 0;
             if( gDebug>0 ){
                printf("WARNING : invalid zen_dist<0 after refraction set zen_dist=0 - perhaps it can be accepted\n");             
             }
          }
          if( zen_dist > 180 ){
             zen_dist = 180;
             if( gDebug>0 ){
                printf("WARNING : invalid zen_dist>180 after refraction set zen_dist=180 - perhaps it can be accepted\n");
             }
          }
          pattern = gAntennaPattern.antenna_gain_from_azh( freq_mhz, skyinfo->azim, zen_dist );
          if( gDebug > 10 ){
             printf("DEBUG : refraction delta_theta = %.2f [secsec] : theta -> theta - delta_theta = %.8f - %.8f = %.8f [deg] -> pattern_new = %.8f (vs old = %.8f)\n",delta_theta_arcsec,zen_dist,delta_theta_deg,(zen_dist - delta_theta_deg),pattern,gAntennaPattern.antenna_gain_from_azh( freq_mhz, skyinfo->azim, zen_dist));
          }                                   
       }       
    }
#else
    if( gIncludeIonosphere ){
       printf("ERROR : program is compiled with precompiler directive _SUPERFAST_OPTIMIZED_ which disables ionospheric corrections, but parameters require it - please verify :\n");
       printf("either recompile the program or disable ionospheric corrections by setting -q include_ionosphere=0\n");
       exit(-1);
    }    
#endif   

    if( gIncludeTropoRefraction > 0 ){
       double delta_theta_deg = tropo_refraction_angle( zen_dist );
       double delta_theta_arcsec = delta_theta_deg*3600.00;

       zen_dist = zen_dist - delta_theta_deg;
       if( zen_dist<0 ){
          zen_dist = 0;
          if( gDebug>0 ){
             printf("WARNING : invalid zen_dist<0 after refraction set zen_dist=0 - perhaps it can be accepted\n");             
          }
       }
       if( zen_dist > 180 ){
          zen_dist = 180;
          if( gDebug>0 ){
             printf("WARNING : invalid zen_dist>180 after refraction set zen_dist=180 - perhaps it can be accepted\n");
          }
       }
       pattern = gAntennaPattern.antenna_gain_from_azh( freq_mhz, skyinfo->azim, zen_dist );
       if( gDebug > 10 ){
          printf("DEBUG : TROPOSPHERIC refraction delta_theta = %.2f [secsec] : theta -> theta - delta_theta = %.8f - %.8f = %.8f [deg] -> pattern_new = %.8f (vs old = %.8f)\n",delta_theta_arcsec,zen_dist,delta_theta_deg,(zen_dist - delta_theta_deg),pattern,gAntennaPattern.antenna_gain_from_azh( freq_mhz, skyinfo->azim, zen_dist));
       }                                         
    } 

    
    // TEST :
//    if( skyinfo->alt < 0 ){
//       pattern = 0.001;
//    }
    
    double delta = pattern*ant_gain_num*gm_filter*delta_sky;

    // Calculation of tau_g for the 2015 paper :
    double tau_g_val = tau_g(zen_dist);
    sum_tau_g += tau_g_val*pattern;
    sum_tau_g_bar += tau_g_val*pattern*delta_sky;
    if( gTestTskyWithTau > 0 ){
       delta = delta * tau_g_val;
    }
//    if( k == 10000*4 ){
//       printf("odo pattern = %.20f\n",pattern);
//    }
    sum += delta;    
    antenna_pattern_sum += pattern_for_pattern_sum;            
    if( skyinfo->alt < 0 ){
       antenna_pattern_sum_below_hor += pattern_for_pattern_sum;       
       sum_below_hor += delta;
    }
    if( T_sun > 0 ){
       sun_sum += pattern*ant_gain_num*gm_filter*T_sun; // to optimize Sun/No Sun into one 
    }
    
    
#ifndef _SUPERFAST_OPTIMIZED_
    // filling antenna pattern map :
    if( healpix_antpat_arr && healpix_deltaTa_arr && healpix_arr ){
       if( gSaveAntPatternDB > 0 ){
          (*healpix_antpat_arr)[npix] = 10.00*log10(pattern);
       }else{
          (*healpix_antpat_arr)[npix] = pattern;
       }        
       (*healpix_deltaTa_arr)[npix] = delta;
       (*healpix_arr)[npix] = skyinfo->intensity;
    }
    
    if( p_zenith_projection ){
       // skyinfo->azim, zen_dist
//       double x = sin(deg2rad(zen_dist))*cos(deg2rad(skyinfo->azim));
//       double y = sin(deg2rad(zen_dist))*sin(deg2rad(skyinfo->azim));

       // LIBNOVA : az starts from S=0,W=90,N=180,E=270
       // X=E, Y=N, -X=W, -Y=S 
       // az -> 360 - az 
       // az' -> az' - 90;
       // az = (360-az) - 90 = 270 - az;       
       double x = sin(deg2rad(zen_dist))*cos(deg2rad(270-skyinfo->azim));
       double y = sin(deg2rad(zen_dist))*sin(deg2rad(270-skyinfo->azim));

       // to make it proper astronomy E - points left -> mirror ?
       x = -x; // is that ok ???
       
       int x_int = (int)((1+x)*gZenithProjRadius);
       int y_int = (int)((1+y)*gZenithProjRadius);
       
//       if( p_zenith_projection->getXY(x_int,y_int) < 0 ){
       if( x_int>=0 && y_int>=0 && x_int<p_zenith_projection->GetXSize() && y_int<p_zenith_projection->GetYSize() ){
          p_zenith_projection->setXY(x_int,y_int,skyinfo->intensity);
          if( gSaveAntPatternDB > 0 ){
             p_zenith_projection_ant->setXY(x_int,y_int,10.00*log10(pattern));
          }else{
             p_zenith_projection_ant->setXY(x_int,y_int,pattern);
          }
          p_zenith_projection_mult->setXY(x_int,y_int,delta);
       }
//       }
    }

//    double theta_deg = 90.00 - (coord.theta/M_PI)*180.00;
//    printf("%d : %.4f (%.8f,%.8f)\n",npix,intensity,phi_deg,theta_deg);

    if( gDebug>0 || ofile ){
       double x_azh,y_azh;
       pointing coord_azh( deg2rad(zen_dist) , deg2rad(skyinfo->azim) );
       map2.project(coord_azh,x_azh,y_azh);
    
       sprintf(outline_short,"%.4f %.4f %.4f %.8Lf %.4f %.4f %.4f %.4f %.8f %.8f\n",skyinfo->xp,skyinfo->yp,skyinfo->intensity,sum,skyinfo->azim,zen_dist,x_azh,y_azh,delta_sky,delta);
       if( gDebug ){
         if( gDebug>=2 || (npix%10000)==0 || npix>=3145710 ){
            sprintf(outline,"%d : %.4f %.8f %.8f %.4f %.4f %.4f %.4f %.4f %.4f %.8f %.8f %.8f %.8Lf\n",npix,skyinfo->intensity,phi_deg,theta_deg,skyinfo->xp,skyinfo->yp,skyinfo->azim,zen_dist,x_azh,y_azh,delta_sky,pattern,delta,sum);
            printf("%s",outline);
         }
       }
    
// yp ??? does not compile - don't know what it is ...
//       if( ofile ){
//         if( fabs(yp) < 1.4141 ){
//           (*ofile) << outline_short;
//         }
//       }    
    }
#endif
    
//    cout << line << endl;
    npix++;
  }
    
  double tau_g_bar = double(sum_tau_g_bar)/double(sum);
  if( gDebug > 0 ){ printf("TEST = %.8f / %.8Lf = %.8f , sum_tau_g_bar = %.8f / %.8Lf = %.8f\n",sum_tau_g,antenna_pattern_sum,double(sum_tau_g)/double(antenna_pattern_sum),sum_tau_g_bar,sum,tau_g_bar); }
  
  sun_sum = sun_sum *dOMEGA;
  sum = sum* dOMEGA;        
  sum_tau_g = double(sum_tau_g)/double(antenna_pattern_sum);
  antenna_pattern_sum = antenna_pattern_sum*dOMEGA;
  antenna_pattern_sum_below_hor = antenna_pattern_sum_below_hor*dOMEGA;
  sum_below_hor = sum_below_hor*dOMEGA;

  // setting values of output parameters :
  out_pattern_sum = antenna_pattern_sum;  
  out_sun_sum = sun_sum;
  
  printf("%.2f MHz : npix=%d , sum = %.8Lf, sum_below_hor = %.8Lf , pattern_sum = %.8Lf , sum_norm = (sum/pattern_sum) = %.8Lf, gain_num=%.2f , pattern_sum_below_hor = %.8f, tau_g = %.8f , tau_g_bar = %.8f \n",freq_mhz,npix,sum,sum_below_hor,antenna_pattern_sum,(sum/antenna_pattern_sum),ant_gain_num,antenna_pattern_sum_below_hor,sum_tau_g,tau_g_bar);
  printf("Pure loop took %d [sec] , uxtime = %d\n",(get_dttm()-start),(int)get_dttm());
  if( gIncludeIonosphere ){
     printf("Min/max ionospheric attenuation = %.8f / %.8f\n",min_iono_att_db,max_iono_att_db);
  }

//  if( sun_sum > 0 ){
  printf("Sun contribution = %.8Lf [K] = %.20Lf/%.20Lf\n",(sun_sum/antenna_pattern_sum),sun_sum,antenna_pattern_sum);
//  }

#ifndef _SUPERFAST_OPTIMIZED_
//  if( gSaveSkyMapFits > 0 && (gSaveSkyMapAtFreq<0 || fabs(gSaveSkyMapAtFreq-freq_mhz)<1 ) ){
  if( healpix_arr && healpix_antpat_arr && healpix_deltaTa_arr ){
     // saving map to file :
     char szSkyMapFile[1024];
     if( strlen(gOutPutDir.c_str()) ){
        sprintf(szSkyMapFile,"%s/%sangelica_%07.2fMHz.fits",gOutPutDir.c_str(),gMapFileBaseName.c_str(),freq_mhz);
     }else{
        sprintf(szSkyMapFile,"%sangelica_%07.2fMHz.fits",gMapFileBaseName.c_str(),freq_mhz);
     }

     printf("Writing skymap to FITS-healpix file %s ...",szSkyMapFile);fflush(stdout);
//   Healpix_Map<double> healpix_skymap( (Healpix_Map<double>)*pHealpix );
     Healpix_Map<double> healpix_skymap(norder,RING);
     healpix_skymap.Set(*healpix_arr,RING);
     string szMapFile=szSkyMapFile;
     fitshandle out_fits;
     out_fits.create(szMapFile);
     write_Healpix_map_to_fits(out_fits,healpix_skymap,PLANCK_FLOAT32);
//  prepare_Healpix_fitsmap(out_fits,
     out_fits.close();
     printf("DONE OK !\n");


     double healpix_ant_sum = 0.00;
     for(int i=0;i<sky_intensities.size();i++){  
        // WARNING : these lines caused program to crash if put after saving antenna pattern to the fits file.
        // I don;t understand why, but apparently HEALPIX is rather poorly written  ... or my code is 
        double pattern_num = (*healpix_antpat_arr)[i];
        if( gSaveAntPatternDB > 0 ){
           pattern_num = db2num( (*healpix_antpat_arr)[i] );
        }
        
        healpix_ant_sum += pattern_num;
     }                            
     printf("TEST : healpix ant_sum = %.2f\n",healpix_ant_sum);
     if( strlen(gOutPutDir.c_str()) ){
        sprintf(szSkyMapFile,"%s/%santpat_%07.2fMHz.fits",gOutPutDir.c_str(),gMapFileBaseName.c_str(),freq_mhz);
     }else{
        sprintf(szSkyMapFile,"%santpat_%07.2fMHz.fits",gMapFileBaseName.c_str(),freq_mhz);
     }
     printf("Writing antenna pattern to FITS-healpix file %s ...",szSkyMapFile);fflush(stdout);
//   Healpix_Map<double> healpix_skymap( (Healpix_Map<double>)*pHealpix );
     Healpix_Map<double> healpix_antpat(norder,RING);
     healpix_antpat.Set(*healpix_antpat_arr,RING);
     szMapFile=szSkyMapFile;
     fitshandle out_fits2;
     out_fits2.create(szMapFile);
     write_Healpix_map_to_fits(out_fits2,healpix_antpat,PLANCK_FLOAT64);
//  prepare_Healpix_fitsmap(out_fits,
     out_fits2.close();
     printf("DONE OK !\n");


     if( strlen(gOutPutDir.c_str()) ){
        sprintf(szSkyMapFile,"%s/%sdTa_%07.2fMHz.fits",gOutPutDir.c_str(),gMapFileBaseName.c_str(),freq_mhz);
     }else{
        sprintf(szSkyMapFile,"%sdTa_%07.2fMHz.fits",gMapFileBaseName.c_str(),freq_mhz);
     }
     printf("Writing delta_Ta to FITS-healpix file %s ...",szSkyMapFile);fflush(stdout);
// first need to normalize it as it is done normally (in calling function):

// CHECK CONSISTANCE of projected (healpix image) integration again final result :
     double test_sum_healpix = 0;
     for(int i=0;i<sky_intensities.size();i++){
//        double norm_val = (*healpix_deltaTa_arr)[i] / healpix_ant_sum;
//        (*healpix_deltaTa_arr)[i] = norm_val;
        (*healpix_deltaTa_arr)[i] = (*healpix_deltaTa_arr)[i] / healpix_ant_sum;
        test_sum_healpix += (*healpix_deltaTa_arr)[i];
     }
     
//   Healpix_Map<double> healpix_skymap( (Healpix_Map<double>)*pHealpix );
     Healpix_Map<double> healpix_deltaTa(norder,RING);
     healpix_deltaTa.Set(*healpix_deltaTa_arr,RING);
     szMapFile=szSkyMapFile;
     fitshandle out_fits3;
     out_fits3.create(szMapFile);
     write_Healpix_map_to_fits(out_fits3,healpix_deltaTa,PLANCK_FLOAT64);
//  prepare_Healpix_fitsmap(out_fits,
     out_fits3.close();
     printf("DONE OK !\n");
     printf("TEST_SUM_HEALPIX = %.2f [K]\n",test_sum_healpix);

     if( p_zenith_projection ){
        mystring szOutFits;

        // need to normalise to pattern sum :
        double test_ant_sum=0;
        for(int y=0;y<p_zenith_projection->GetYSize();y++){
           for(int x=0;x<p_zenith_projection->GetXSize();x++){
              double val = p_zenith_projection_ant->getXY(x,y);
              if( gSaveAntPatternDB > 0 ){
                 val = db2num( val );
              }                                 
              test_ant_sum += val;
           }
        }
        double test_sum=0.00,test_check_sum=0.00;
        for(int y=0;y<p_zenith_projection->GetYSize();y++){
           for(int x=0;x<p_zenith_projection->GetXSize();x++){           
              double val = p_zenith_projection_mult->getXY(x,y);
              test_check_sum += val;
              double new_val = val / test_ant_sum;
              p_zenith_projection_mult->setXY( x, y, new_val );
              test_sum += new_val;
           }
        }        
        printf("TEST_SUM = %.2f [K] (vs %.2f [K]) (in zenithal projection, the difference from final result comes from projection errors etc)\n",test_sum,(test_check_sum/test_ant_sum));
        
        szOutFits << gOutPutDir.c_str() << "/" << gZenithProjBaseName.c_str() << ".fits";
        printf("Saving SKYMODE zenith projection to fits file %s ...\n",szOutFits.c_str());fflush(stdout); 
        p_zenith_projection->WriteFits(szOutFits.c_str());
        
        szOutFits="";
        szOutFits << gOutPutDir.c_str() << "/" << gZenithProjBaseName.c_str() << "_ant.fits";
        printf("Saving ANT zenith projection to fits file %s ...\n",szOutFits.c_str());fflush(stdout);
        p_zenith_projection_ant->WriteFits(szOutFits.c_str());

        szOutFits="";
        szOutFits << gOutPutDir.c_str() << "/" << gZenithProjBaseName.c_str() << "_mult.fits";
        printf("Saving MULT zenith projection to fits file %s ...\n",szOutFits.c_str());fflush(stdout);
        p_zenith_projection_mult->WriteFits(szOutFits.c_str());
     }   
                 
     if( healpix_arr )
        delete healpix_arr;
     if( healpix_antpat_arr )
        delete healpix_antpat_arr;
     if( healpix_deltaTa_arr )
        delete healpix_deltaTa_arr;               
     if( p_zenith_projection )
        delete p_zenith_projection;
  }else{
     printf("INFO : saving of skymap to fits file is not required\n");
  }  
#else 
  printf("WARNING : saving of skymap is not compiled at all\n");
#endif

  return sum;
}

double calc_spec( const char* infile, double unix_time, int save_map, int first_freq );

void print_parameters()
{
  printf("showspec program integrating sky model with antenna pattern to obtain spectrum. Version : %s Started at unixtime = %d \n",SHOWSPEC_VERSION,get_dttm());
  printf("#################################################\n");
  printf("PARAMTERS :\n");
  printf("#################################################\n");
  printf("list of infiles  = %s\n",infile_list.c_str());
  printf("Object list file = %s\n",gObjectListFile.c_str());
  printf("Start UT         = %s\n",get_gmtime_string_bg(start_ux_time).c_str());
  printf("Time interval    = %d [sec]\n",time_interval);
  printf("Time step        = %.8f [sec]\n",dt);
  printf("Save all maps    = %d\n",gDoSaveAllMaps);
  printf("Input directory  = %s\n",gInDir.c_str());
  printf("Binary input     = %d\n",CSkyMapCache::m_bBinaryFile);
  printf("Antenna pattern type :");
  if( strlen(CAntPattern::m_AntennaPatternFile.c_str()) ){
     printf("From file %s\n",CAntPattern::m_AntennaPatternFile.c_str());
  }
  if( strlen(CAntPattern::m_AntennaPatternFile.c_str())==0 || CAntPattern::gAntennaPatternFormula>0 ){
     string szFormulaName="DIPOL OVER GND-SCREEN";
     if( gAntennaPattern.m_AntennaPatternType == eAntPattGaussian ){
        szFormulaName="GAUSSIAN exp(-(theta/theta0)^2)";
     }
     printf("According to formula %d (%s)\n",gAntennaPattern.m_AntennaPatternType,szFormulaName.c_str());
     printf("\t\tAntenna beam width = %.2f [deg]\n",gAntennaPattern.m_AntennaBeamWidth);
  }
  printf("Suppress above zenithal angle = %.2f [deg]\n",CAntPattern::gSuppressAntPatternAboveZenAngle);
  printf("Antenna cache ON = %d\n",gAntennaPattern.m_bCacheON);
  printf("Theta max        = %.2f [deg]\n",gThetaMax);
  printf("weather file     = |%s|\n",gWeather.m_szWeatherFile.c_str());
  printf("T_amb            = %.2f [K]\n",gT_amb);
  printf("Site location    = (%.4f,%.4f) [deg] [%s]\n",geo_long,geo_lat,gSiteName.c_str());
  printf("gAddGalaxy       = %d\n",gAddGalaxy);
  printf("Enable antenna pattern interpolation = %d\n",CAntPatternParser::m_bEnableInterpolation);
  printf("Antenna pattern file = %s\n",CAntPattern::m_AntennaPatternFile.c_str());
  printf("Only above horizon = %d\n",gOnlyAboveHorizont);
  printf("Debug           = %d\n",gDebug);
  printf("Include ionosphere = %d\n",gIncludeIonosphere );
  printf("\t\tF-layer refraction = %d\n",gIncludeIonoRefraction);
  printf("\t\t\tTypical refraction at 100 MHz = %.2f [arcmin]\n",gTypicalRefractionAngleArcMin);
  printf("\t\t\tTypical refraction at 100 MHz at zenith angle = %.2f [deg]\n",gTypicalRefractionZenithAngle);
  printf("\t\tD-layer absorption = %d\n",gIncludeIonoAbsorption);  
  printf("\t\tD-layer electron emission = %d\n",gIncludeIonoElectronTemp);  
  printf("\t\tTEC (total)        = %e [e/m^2]\n",gTEC_total);  
  printf("\t\tn_e                = %e [1/m^3]\n",gN_e);
  printf("\t\tplasma freq        = %.4f [MHz]\n",gPlasmaFreq/1e6);   
  printf("\t\tT_electron (D-layer) = %.2f [K]\n",gElectronTemp);
  printf("\t\tD-layer height     = %.2f [m]\n",gHeightD);
  printf("\t\tD-layer width      = %.2f [m]\n",gDeltaHeightD);
  printf("\t\tD-layer Earth radius = %.2f [m]\n",gRadiusEarth);
  printf("\t\tTypical loss at 100 MHz = %.2f [dB]\n",gTypicalLossAt100MHz);
  printf("Include troposhere (refraction) = %d\n",gIncludeTropoRefraction);
  printf("\t\tPressure = %.2f [milibars = hPa]\n",gAtmPressure);     
  printf("Output dir      = %s\n",gOutPutDir.c_str());
  printf("Output spec file = %s\n",gOutSpecFileName.c_str());
  printf("Save output to fits = %d\n",gSaveFitsFile);
  printf("Zenith projection resolution = %d x %d\n",(2*gZenithProjRadius+1),(2*gZenithProjRadius+1));
  printf("Antenna efficiency included = %d\n",gAddAntEffCorrection);
  printf("Antenna efficiency file     = %s\n",gAntEffFile.c_str());
  printf("Antenna rotation angle      = %.2f [deg]\n",CAntPattern::m_RotFromAxis);
  printf("Norm to 4PI                 = %d\n",gNorm4PI);  
  printf("Include Sun model = %d\n",gAddSun);
  printf("Frequency range = %.2f - %.2f [MHz]\n",gFreqStart,gFreqEnd);
  printf("Save sky maps to fits = %d , fits prefix = |%s|  ",gSaveSkyMapFits,gMapFileBaseName.c_str());  
  if( gSaveSkyMapAtFreq > 0 ){
     printf("( at frequency = %.2f [MHz] )\n",gSaveSkyMapAtFreq);  
  }else{
     printf("( at all frequencies )\n");
  }
  printf("\t\tgSaveAntPatternDB = %d\n",gSaveAntPatternDB);
  printf("Caching of skymodel enabled = %d\n",gSkyMapCache.GetCacheONOFF());
  printf("Maximum value of theta (simulated in FEKO) = %.2f [deg]\n",CAntPatternParser::m_SimulMaxTheta);
  printf("Auto-detetction of theta is                = %d [0-disabled,1-enabled]\n",CAntPatternParser::m_AutoDetectTheta);
  printf("Automatic continuation enabled = %d\n",gAutoContinue);
  if( strlen(gUxListFile.c_str()) > 0 ){
     printf("Unix list file = %s\n",gUxListFile.c_str());
  }else{
     printf("INFO : list file of unix times not provided\n");
  }
  printf("TEST of integral of T_sky(theta,phi)*P(theta,phi)*Tau(theta) = %s\n",(gTestTskyWithTau ? "ENABLED" : "DISABLED" ));
  printf("#################################################\n");
}

void catch_ctrlc (int sig)
{
  printf("Quitting\n");
  //  signal (sig, catch_ctrlc);
  printf("Program execution took : %d sec\n",(get_dttm()-program_start));
  exit(0);  
}
      

int main(int argc,char* argv[])
{
  if(argc >= 2){
    infile_list = argv[1];    
  }else{
    usage();
  }

  // parse command line :  
  parse_cmdline(argc,argv);
  print_parameters();
  
  // setting signal handlers :
  signal(SIGINT, catch_ctrlc);
  signal(SIGTERM, catch_ctrlc);

  if( gIncludeIonoRefraction != eNoIonoRefraction ){
     test_refraction_angle();
  }      

  if( strlen(gObjectListFile.c_str()) ){
     read_sources(gObjectListFile.c_str());
  }

  // initialization of global variables :
  nside = Healpix_Base::npix2nside(12*512*512);
  norder = Healpix_Base::nside2order( nside );
  pHealpix = new Healpix_Base(norder,RING);// RING OR NEST - which one ??? chyba raczej RING 
             // dump gives reasonable results : means at phi=0.00/360.00 [deg] there high values >10,000 and at 90,270 [deg] there
             // are smallest values , which agrees with galactic coordinates Long=0,Lat~0 ( Center of the Galaxy )
  dOMEGA = ((long double)(4*M_PI))/((long double)pHealpix->Npix());;
  

  // reading weather file if required 
  if( strlen(gWeather.m_szWeatherFile.c_str()) ){
     int weather_points = gWeather.ReadTempFile();
     if( weather_points <= 0 ){
        printf("ERROR : could not read weather file %s, option -e required -> cannot continue !\n",gWeather.m_szWeatherFile.c_str());
        exit(-1);
     }else{
        printf("OK : read %d points from weather file %s\n",weather_points,gWeather.m_szWeatherFile.c_str());
     }
  }

  if( MyFile::DoesFileExist(infile_list.c_str()) ){
     printf("INFO : using local ni-list file |%s|\n",infile_list.c_str());
  }else{
     if( strlen(gInDir.c_str()) > 0 ){
        char szTmp[1024];
        sprintf(szTmp,"%s/%s",gInDir.c_str(),infile_list.c_str());
        infile_list = szTmp;
     }
  }

  if( strlen(gOutPutDir.c_str()) > 0 ){
    string szMkDir="mkdir -p ";
    szMkDir += gOutPutDir;
    printf("command : %s\n",szMkDir.c_str());
    system(szMkDir.c_str());
  }

  vector<string> infiles_list;
  if( bg_read_list( infile_list.c_str(), infiles_list ) <= 0 ){
    printf("ERROR : list file %s not found, cannot continue\n",infile_list.c_str());
    exit(-1);
  }
  AutoDetectExtensions( infiles_list );
  
  vector<double> ni_list;
  for(int i=0;i<infiles_list.size();i++){
     double freq_mhz = parse_ni( infiles_list[i].c_str() );
     
     if( gFreqStart<=freq_mhz && freq_mhz<=gFreqEnd ){
        ni_list.push_back(freq_mhz);
     }
  }
  
  // initialization of antenna pattern structure :
  gAntennaPattern.Init( ni_list );  

  int n_steps = (int)(time_interval / dt);
  vector<cValue> unix_times_to_process;
  if( strlen(gUxListFile.c_str()) > 0 ){
     if( read_file(gUxListFile.c_str(),unix_times_to_process)  > 0 ){
        n_steps = unix_times_to_process.size();
     }
  }else{
     double t = start_ux_time;
     while( t <= (start_ux_time+time_interval) ){
        cValue tmp;
        tmp.x = t;

        unix_times_to_process.push_back(tmp);
        t = t + dt; 
     }
  }

  // 
  CBgFits* pOutFits = NULL;

  // LOOP OVER TIME :   
  double first_freq=-1,delta_freq_mhz=-1;
  int i=0;
//  double t = start_ux_time;
//  while( t <= (start_ux_time+time_interval) ){
  for(int integration=0;integration<unix_times_to_process.size();integration++){
    double t = unix_times_to_process[integration].x;
  
    printf("################################################# t = %.8f #################################################\n",t);
    if( t < start_ux_time ){
       printf("SKIPPED : unix time %.8f < start_time = %.8f -> ignored\n",t,start_ux_time);
       continue;
    }

    // LOOP OVER FREQUENCIES :
    time_t start = get_dttm(); 
    CBgArray spectrum;
    int bFirstFreq=1;
    for(int f=0;f<infiles_list.size();f++){      
      const char* ni_file_tmp = infiles_list[f].c_str();
      char ni_file[1024];
      strcpy(ni_file,ni_file_tmp);
      if( strlen(gInDir.c_str()) > 0 ){
         sprintf(ni_file,"%s/%s",gInDir.c_str(),ni_file_tmp);
      }
      double freq_mhz = parse_ni( ni_file );

      // automatically checking frequency resolution       
      if( f==0 ){
         first_freq = freq_mhz;
      }else{
         if( f==1 ){
            delta_freq_mhz = (freq_mhz-first_freq);
            printf("AUTODETECTED FREQUENCY RESOLUTION = %.8f [MHz]\n",delta_freq_mhz);
         }
      }
      
      if( gFreqStart<=freq_mhz && freq_mhz<=gFreqEnd ){
         printf("------------------------------------------------------ %.2f [MHz] ------------------------------------------------------\n",freq_mhz);fflush(stdout);
         printf("Analysing input file = |%s|\n",ni_file);
      
         time_t start=get_dttm();
         double Ta = calc_spec( ni_file, t, (i==0), bFirstFreq );
         printf("Integration of sky image took %d [sec] , T_a = %.2f [K]\n",(int)(get_dttm()-start),Ta);
         printf("------------------------------------------------------\n");fflush(stdout);
         bFirstFreq=0;
         spectrum.push_back(Ta);
      }else{
         if( gDebug > 0 ){
            printf("Frequency %.2f MHz is outside required regime of %.2f - %.2f [MHz]\n",freq_mhz,gFreqStart,gFreqEnd);
         }
      }
    }    
    printf("Loop over frequencies took %d [sec]\n",(get_dttm()-start));
    
    if( strlen(gFitsFile.c_str()) > 0 ){
       if( !pOutFits ){
          pOutFits = new CBgFits(spectrum.size(),(time_interval/dt)+1);          
       }       
       pOutFits->add_line( spectrum );
    }
    
    t = t + dt;
    i++;
    
    if( pOutFits ){
       if( n_steps > 100 ){
          if( (i%100) == 0 ){
             pOutFits->PrepareBigHornsHeader( start_ux_time, dt, first_freq, delta_freq_mhz );
             pOutFits->WriteFits( gFitsFile.c_str(), 1 );
          }
       }
    }
  }

  if( pHealpix ){
     delete pHealpix;
  }

  if( pOutFits ){
     pOutFits->PrepareBigHornsHeader( start_ux_time, dt, first_freq, delta_freq_mhz );
     pOutFits->WriteFits( gFitsFile.c_str(), 1 );
     delete pOutFits;
  }

  printf("Program execution took : %d sec\n",(get_dttm()-program_start));

  return 0;
}  

double parse_ni( const char* infile )
{
  // parse FREQ from input file name :
  char szNI[64];
  double gNI=0.00;
  char* ni = (char*)strcasestr(infile,"ni");
  const char* dot = strstr(ni,".");
  if( ni && dot ){
    ni += 2;
    memset(szNI,'\0',64);
    int i=0,dot_count=0;
    while(is_digit(ni[i]) || (ni[i]=='.' && !dot_count) ){
      szNI[i]=ni[i];
      if( ni[i]=='.' )
        dot_count++;
      i++;
    }
    gNI = atof(szNI);
  }

  return gNI;
}

double calc_spec( const char* infile, double unix_time, int save_map, int first_freq )
{  
  // parse FREQ from input file name :
  double gNI = parse_ni( infile );
    

  string szUT=get_gmtime_string_bg((time_t)unix_time);
  char tmp_str[1024],outfile_buffer[8192];
  string outfile;
  getbasename_new(infile,outfile);
//  sprintf(tmp_str,"%s_%s.ang",outfile.c_str(),szUT.c_str());
  sprintf(tmp_str,"%s_%.2fmhz.ang",szUT.c_str(),gNI);  
  outfile=tmp_str;  

  // overwriting save_map in case all maps are required to be saved :
  save_map = (save_map && gDoSaveAllMaps);  // was || but changed to have parameter to disable maps saving ( take too much space on device ! )

  // file .ang with map of galactic radiation map 
  ofstream* ofile=NULL;
  if( save_map ){
    string szOutFileFullPath;
    if( strlen(gOutPutDir.c_str()) ){
      szOutFileFullPath = gOutPutDir.c_str();
      szOutFileFullPath += "/";
    }
    szOutFileFullPath += outfile;
    printf("DEBUG : saving sky map to file |%s|\n",szOutFileFullPath.c_str());
    
    ofile = new ofstream( szOutFileFullPath.c_str() );
    ofile->rdbuf()->pubsetbuf(outfile_buffer,8192);
  }else{
    printf("WARNING : saving of sky map is not required\n");
  }

  char szTimeSnapshotFile[1024];
  string szPlotFile;
  if( strlen(gOutPutDir.c_str()) ){
  
    if( strlen(gOutSpecFileName.c_str()) > 0 ){
       sprintf(szTimeSnapshotFile,"%s/%s",gOutPutDir.c_str(),gOutSpecFileName.c_str());
    }else{
       sprintf(szTimeSnapshotFile,"%s/%s.spec",gOutPutDir.c_str(),szUT.c_str());
    }
  }else{
    if( strlen(gOutSpecFileName.c_str()) > 0 ){
       sprintf(szTimeSnapshotFile,"%s",gOutSpecFileName.c_str());
    }else{
       sprintf(szTimeSnapshotFile,"%s.spec",szUT.c_str());
    }
  }  
  change_ext(szTimeSnapshotFile,"plot",szPlotFile);
  if( gSaveFitsFile > 0 ){
     change_ext(szTimeSnapshotFile,"fits",gFitsFile);
  }

  // reading previous file to auto-continue 
  if( gAutoContinue > 0 ){
     // check if .plot file exists :
     if( MyFile::DoesFileExist(szPlotFile.c_str()) ){
        vector<cValue> out_list;
        out_list.reserve(3000000);
        int prev_results_count = read_file(szPlotFile.c_str(),out_list,0);
        
        printf("WARNING : Auto-continue option is enabled and %d points read from file %s\n",prev_results_count,szPlotFile.c_str());
        if( prev_results_count > 0 ){
           gLastUxTimeProcessed = out_list[out_list.size()-1].z;
           gLastFreqProcessed   = out_list[out_list.size()-1].x;
           
           printf("WARNING : Last processed line %.8f MHz - %.2f K - at uxtime %.2f \n",gLastFreqProcessed,out_list[out_list.size()-1].y,gLastUxTimeProcessed);
        }
        
        gAutoContinue = -gAutoContinue; // used not to re-read the whole results file again and again 
     }
  }
  if( gLastUxTimeProcessed > 0 ){
     if( (time_t)unix_time < (time_t)gLastUxTimeProcessed ){
        printf("WARNING : unix_time %d / freq %.8f MHz skipped due to auto-continue option detected last processed unix_time = %.2f , freq = %.8f MHz\n",unix_time,gNI,gLastUxTimeProcessed,gLastFreqProcessed);
        return -1000; // skipped 
     }else{
        if( (time_t)unix_time == (time_t)gLastUxTimeProcessed ){
           if( gNI <= gLastFreqProcessed ){ // 2014-02-17 : < changed to <= : changed in order to avoid line repetition in the output plot/spec file, it will be tested on HPC and old version will be restored in case something is wrong with the new one
              printf("WARNING : unix_time %d / %.8f MHz skipped due to auto-continue option detected last processed unix_time = %.2f / %.8f MHz\n",unix_time,gNI,gLastUxTimeProcessed,gLastFreqProcessed);
              return -1000.00;
           }else{
              printf("WARNING : unix_time %d / %.8f MHz not yet processed ( last processed unix_time = %.2f / %.8f MHz ) - restarting processing now !\n",unix_time,gNI,gLastUxTimeProcessed,gLastFreqProcessed);
              gLastUxTimeProcessed = -1000;
              gLastFreqProcessed   = -1000;                         
              first_freq = 1; // forcing first_freq flag to be set in order to recalculate horizontal coordinates
           }
        }else{
           printf("WARNING : unix_time %d / %.8f MHz not yet processed ( last processed unix_time = %.2f / %.8f MHz\n",unix_time,gNI,gLastUxTimeProcessed,gLastFreqProcessed);
           gLastUxTimeProcessed = -1000;
           gLastFreqProcessed   = -1000;
           first_freq = 1; // forcing first_freq flag to be set in order to recalculate horizontal coordinates
        }
     }
  }

  double norm_pattern_sum=0.00;
  double sum = 0.00, sun_sum=0.00;
  if( gRunFastVersion ){
    sum = integrate_galactic_spectrum_FAST( infile, unix_time, ofile, norm_pattern_sum, sun_sum, first_freq );
  }else{
    sum = integrate_galactic_spectrum( infile, unix_time, ofile, norm_pattern_sum, sun_sum, first_freq );
  }

  if( ofile ){
    ofile->close();
    delete ofile;
  }
  printf("TOTAL SUM = %.8f [K]\n",sum);  
  
  ofstream spec_file(szTimeSnapshotFile,ios::app);
  ofstream plot_file(szPlotFile.c_str(),ios::app);
  if( first_freq ){
     // for first frequency save current position of Galactic center :
     // calculate (AZIM,ALT) of galactic center :
     double azim_gal_c,alt_gal_c,ra_gal_c,dec_gal_c;
     char szComment[1024];
     gal2hor( 0, 0, geo_long, geo_lat, unix_time, azim_gal_c,alt_gal_c,ra_gal_c,dec_gal_c );
     sprintf(szComment,"# Galactic center (0,0) is at %s UT at (AZIM,ALT) = (%.4f , %.4f) [deg]",szUT.c_str(),azim_gal_c,alt_gal_c);            
     spec_file << szComment << endl;
     
     sprintf(szComment,"# Freq[Mhz]  T_total PatternSum T_antena[K] T_antena*k_b[mW/Hz] T_antena*k_b[dBm/Hz] UnixTime T_ant-sun[K]");
     spec_file << szComment << endl;
          
     printf("%s\n",szComment);    
  }
  
  char szLine[256];
  double T_a = (sum/norm_pattern_sum);
  double T_sun = (sun_sum/norm_pattern_sum);
  double T_point_sources = 0.00;

  // add point sources - to be tested :
  if( gObjectList.size() ){  
     double T_point_sources=0.00;
     for(int i=0;i<gObjectList.size();i++){
        cObjectInfo& object = gObjectList[i];
        
        // calculate (az,h) coordinates from (ra,dec) -> antenna gain in this direction -> antenna effective area -> received flux -> contribution to antenna efficiency :
        double az,alt,z;
        radec2azh( object.ra, object.dec, unix_time, geo_long, geo_lat, az, alt );
        z = (90.00 - alt);

        double freq_mhz = gNI;
        double gain = gAntennaPattern.antenna_gain_from_azh( freq_mhz, az, z );
//        double lamda = 300.00 / freq_mhz;
//        double a_eff = ((lamda*lamda)/(4*M_PI))*gain;
        
        double delta_temp = jy2kelvin_gain( object.flux, freq_mhz, gain );        

        // adding extra temperature from this point source :
// TODO : the extra signal should be weighted by length of the signal vs integration time :
// delta_temp should be a flux averaged over integration time !
// actually every integration should be T_ant = (T_int_start + T_int_end)/2.00 - or rather calculate for (t_start + inttime/2)
        T_a += delta_temp; 
        T_point_sources += delta_temp;
     }
  }
  
  // if T_a < T_ambient -> it is really T_ambient (temperatures do not add ?)
//  if( T_a < gT_amb ){
//     printf("DEBUG : T_a = %.2f [K] < ambient temperature, set to T_ambient = %.2f [K]\n",T_a,gT_amb);
//     T_a = gT_amb;
// }
  
  if( gNorm4PI ){
     printf("WARNING : normalization to 4PI required - does it work ???!!!\n");
     T_a = sum/(4*M_PI);
  }
  
  if( gAddAntEffCorrection>0 ){
     // take account ambient temp and ant. efficiency :
     double ant_eff = get_ce300e_eff(gNI);
     if( gAntEffTable.size() > 0 ){
        ant_eff = interpolate( gAntEffTable, gNI );
     }else{
        printf("WARNING : no antenna efficiency data provided, using default values for BICONICAL ANTENNA !\n");
     }
     
     double T_amb = gT_amb;
     if( gWeather.GetCount() ){
        T_amb = gWeather.GetAmbTempK(unix_time); 
     }
     printf("ANT-EFF CORRECTION : AntEff Corr(%.2f MHz) = %.2f -> T_a := %.2f * %.2f + %.2f*(1-%.2f)\n",gNI,ant_eff,ant_eff,T_a,T_amb,ant_eff);
     T_a = ant_eff*T_a + T_amb*(1.00-ant_eff);
  }
  
  double power_per_hz = T_a*k_b; // was sum
  double power_mW_per_hz = power_per_hz*1000.00;
  sprintf(szLine,"%.8f %.8f %.8f %.8f %E %.8f %.4f %.8f\n",gNI,sum,norm_pattern_sum,T_a,power_mW_per_hz,mW2dbm(power_mW_per_hz),unix_time,T_sun); // %E
  spec_file << szLine;
  
  sprintf(szLine,"%.8f %.8f %.4f %.8f %.8f\n",gNI,T_a,unix_time,get_local_sidereal_time(unix_time),T_sun);
  plot_file << szLine;

  return T_a;
}

double integrate_galactic_spectrum_FAST( const char* infile , double unix_time, ofstream* ofile, double& out_pattern_sum, double& out_sun_sum, int first_freq  )
{
  // uses or not cache - depending on cache_on=1 or =0 parameter (by default disabled)
  vector<cSkyInfo>* pSkyIntensities = gSkyMapCache.get_skyintensities( infile );  
  if( !pSkyIntensities ){
     printf("ERROR : could not get skyintensities from in-file = %s\n",infile);
     exit(-1);     
  }
  vector<cSkyInfo>& sky_intensities = (*pSkyIntensities);

//  vector<cSkyInfo>& sky_intensities = CSkyMapCache::m_SkyIntensities;
//  CSkyMapCache::read_infile( infile );

  string line;

  MollweideSkyMap map2(100);   
 
  double JD_time = ux2jd( unix_time );
  double freq_mhz = parse_ni( infile );
//  double ant_gain_num = get_ce300e_gain(freq_mhz);
  double ant_gain_num = 1.00;

// MODELLING 2012-06-27
//  double gm_filter = get_gm_transmiss(freq_mhz);
  double gm_filter = 1.00;

  double t_amb = gT_amb;
  if( gWeather.GetCount() > 0 ){
     t_amb = gWeather.GetAmbTempK(unix_time);
  }
       
  printf("First freq flag = %d\n",first_freq);
  printf("Npix=%d -> nside=%d -> norder=%d -> npix=%d , dOMEGA=%.10Lf [srad] (OK?)\n",(12*512*512),(int)nside,(int)norder,(int)pHealpix->Npix(),dOMEGA);
  printf("Freq = %.2f [MHz], ant_gain = %.4f, gm_filter = %.4f , t_amb = %.2f [K]\n",freq_mhz,ant_gain_num,gm_filter,t_amb);
        
  long double sum=0.0000;  
  char outline[1024],outline_short[256];
  long double antenna_pattern_sum = 0.00;
  long double sun_sum=0.0000;

  int bAddSun = gAddSun;
  double sun_ra, sun_dec, sun_az, sun_alt;
  if( bAddSun ){
     // check if above horizon :
     get_sun_all( unix_time, geo_long, geo_lat, sun_ra, sun_dec, sun_az, sun_alt );
     
     if( sun_alt <= -1.00 ){
         // ignore Sun if deep below horizon :
         bAddSun = 0;
     }   
  }

  double min_iono_att_db = 1000000.00;
  double max_iono_att_db = -1000000.00;

  int sky_intensities_count = sky_intensities.size();

/*  if( first_freq > 0 ){
     printf("INFO : first frequency for given time (JD_time=%.8f) -> initializing (AZIM,ALT) data in sky_intensities array ...",JD_time);fflush(stdout);
     time_t start = get_dttm();
     for(vector<cSkyInfo>::iterator it=sky_intensities.begin();it!=sky_intensities.end();it++){
        cSkyInfo& skyinfo = (*it);
        
       // at first frequency for given UNIX_TIME fill coorindates info as sky changes in time :
       // using libnova - means in convention of AZIM : S=0, W=90, N=180, E=270 deg 
       //       gal2hor( skyinfo.glon, skyinfo.glat, geo_long, geo_lat, unix_time, skyinfo.azim, skyinfo.alt, skyinfo.ra, skyinfo.dec );
       gal2hor_fromJD( skyinfo.glon, skyinfo.glat, geo_long, geo_lat, JD_time, skyinfo.azim, skyinfo.alt, skyinfo.ra, skyinfo.dec );
     } 
     printf("took %d sec\n",(get_dttm()-start));fflush(stdout);
  }*/
  
  // get antenna pattern array (map) :
  time_t start = get_dttm();
  printf("Getting antenna pattern map ...");fflush(stdout);
//  double* antenna_pattern_array = gAntennaPattern.antenna_pattern( infile,  freq_mhz, sky_intensities );
/*  if( gDebug >= 3 ){
     printf("Antenna gain at %.2f MHz debug :",freq_mhz);
     for(int ll=0;ll<sky_intensities.size();ll+=10000){
        printf(" %.20f",antenna_pattern_array[ll]);       
     }
     printf("\n");
  }*/
  int freq_index = gAntennaPattern.FindFreq(freq_mhz);
  CFreqPattern* pPatternFreq = gAntennaPattern.GetFreqPattern( freq_index );
  CFreqPattern* pPatternFreqNext = gAntennaPattern.GetFreqPattern( freq_index+1 );                           
  if( !pPatternFreq || !pPatternFreqNext ){
     printf("ERROR in code : could not find antenna pattern for frequency index %d or %d (freq_mhz=%.2f MHz)\n",freq_index,(freq_index+1),freq_mhz);
     exit(-1);
  }                                       
  printf("took %d sec\n",(get_dttm()-start));fflush(stdout);

  
  start = get_dttm();

  // healpix structure for small number of points, very crude granulation :
  int npix_fast=12288; // 768 , 3072, only 12 * N^2 allowed N - 2^k = 2,4,8,16,32 etc etc 
//  int npix_fast = 3145728;
  int nside_fast = Healpix_Base::npix2nside(npix_fast);
  int norder_fast = Healpix_Base::nside2order( nside_fast );
  Healpix_Base healpix_fast(norder_fast,RING);
  printf("Healpix_fast(%d,%d)\n",norder_fast,nside_fast);
  printf("Starting loop over %d elements ...",sky_intensities_count);fflush(stdout);  
  
  
  int npix=0;
  while( npix < npix_fast ){
    pointing coord = healpix_fast.pix2ang(npix);
    double phi_deg = coord.phi * rad2deg_const;
    double theta_deg = coord.theta * rad2deg_const;
    double glon = phi_deg; // ???
    double glat = (90.00 - theta_deg);
    
    double azim,alt,ra,dec;
    gal2hor_fromJD( glon, glat, geo_long, geo_lat, JD_time, azim, alt, ra, dec );      
                               
    if( theta_deg >= 180.00 || npix>=3145728 ){	
       printf("WARNING : last line skipped !\n");
       phi_deg = 0.00;
       break;
    }
    double zen_dist = (90.00-alt);


    int npix_global = pHealpix->ang2pix(coord);
    double delta_sky = sky_intensities[npix_global].intensity;
    if( gOnlyAboveHorizont >= 0 ){
       if( alt < gOnlyAboveHorizont /* && alt>0 TEST with signal from underground added - according to antenna-over-gnd-screen-pattern */ ){
          // skip below horizont if required 
//          delta_sky = 0.00;
          delta_sky = t_amb; // 290K below minium elevation 
       }

       if( CAntPatternParser::m_SimulMaxTheta < 180 ){
         if( zen_dist > CAntPatternParser::m_SimulMaxTheta ){
            delta_sky = 0;
         }
       }
    }

    // OPTIMIZED Solar code (should be fast and OK) :
    if( bAddSun >0 && alt>0 ){ // if Sun is enabled in the code and pixel Alt>0 :
       double fabs_dec = fabs_inline( dec );
       if( fabs_dec < 26.00 ){ // if close enough from the celestial Equator 
          // (RA1-RA2)=10deg -> angular distance = 8.90773 deg at DEC=27 deg 
          // (RA1-RA2)= 5deg -> angular distance = 4.45474 deg at DEC=27 deg          
          // this preliminary check is to avoid calculation of angular distance using sin/cos (CPU expensive) calls :
          double ra_diff = fabs_inline(ra - sun_ra);
          double dec_diff = fabs_inline(dec - sun_dec);
          
          if( ra_diff<5.00 && dec_diff<5.00 ){
             // only in case we are close to Sun call get_sun_brightness (expensive sin/cos calculations)
             double T_sun = get_sun_brightness( ra, dec, azim, alt, sun_ra, sun_dec );
             delta_sky += T_sun;
             sun_sum += T_sun;
          }
       }
    }
    
    delta_sky += gExternalNoiseT;


    // INTEGRAL ~= SUM :
//    double pattern = antenna_pattern_array[npix_global];
    double pattern = gAntennaPattern.antenna_pattern( freq_mhz, pPatternFreq, pPatternFreqNext, azim, zen_dist );
    double delta = pattern*ant_gain_num*gm_filter*delta_sky;
    sum += delta;    
    antenna_pattern_sum += pattern;    
        
    npix++;
  }

  sun_sum = sun_sum;
  sum = sum;
  antenna_pattern_sum = antenna_pattern_sum;

  // setting values of output parameters :
  out_pattern_sum = antenna_pattern_sum;  
  out_sun_sum = sun_sum;
  
  printf("npix=%d , sum = %.8Lf, pattern_sum = %.8Lf , sun_norm = (sum/pattern_sum) = %.8Lf, gain_num=%.2f\n",npix,sum,antenna_pattern_sum,(sum/antenna_pattern_sum),ant_gain_num);
  printf("Pure loop took %d [sec] , uxtime = %d\n",(get_dttm()-start),(int)get_dttm());
  if( gIncludeIonosphere ){
     printf("Min/max ionospheric attenuation = %.8f / %.8f\n",min_iono_att_db,max_iono_att_db);
  }

  if( sun_sum > 0 ){
     printf("Sun contribution = %.8Lf [K] = %.20Lf/%.20Lf\n",(sun_sum/antenna_pattern_sum),sun_sum,antenna_pattern_sum);
  }

  return sum;
 
}

void AutoDetectExtensions( vector<string>& list_of_files )
{
   int bIsBinary=1;
   for(int i=0;i<list_of_files.size();i++){
      const char* file = list_of_files[i].c_str();
      
      if(!strstr(file,".dat")){
         bIsBinary=0;
      }
   }

   CSkyMapCache::m_bBinaryFile = bIsBinary;
   printf("INFO : auto-detected file format is %s\n",( CSkyMapCache::m_bBinaryFile ? "BINARY" : "TEXT" ));
}

// calculation based on difference of y - x and cosnines theorem in the 2015-03 notebook 
double delta_s( double theta_deg )
{
   double theta_rad = theta_deg*(M_PI/180.00);
   double cos_theta = cos(theta_rad);
   double ds = gDeltaHeightD * (1.00 + gHeightD/gRadiusEarth) / sqrt(cos_theta*cos_theta + (2.00*gHeightD)/gRadiusEarth );
         
   return ds;
}

double tau_g( double theta_deg )
{
   // variation of opacity with zenith angle :
   /*double gHeightD = 75e3; // [m] = 75km 60-90km
   double gDeltaHeightD = 30e3; // [m] 30km
   double gRadiusEarth = 6371e3; // [m] = 6300 km*/
 
   double Re = gRadiusEarth;
   double H = gHeightD;
   double h = gDeltaHeightD;

//   double theta_deg = x[0];
   double theta_rad = theta_deg*(M_PI/180.00);
   double sqrt1 = sqrt( (Re*cos(theta_rad))*(Re*cos(theta_rad)) + 2.00*Re*H + H*H );
   double ratio1 = (2.00*(Re + h) + h)/((Re*cos(theta_rad))*(Re*cos(theta_rad)) + 2.00*Re*H + H*H );
   double sqrt2 = sqrt( 1.00 + h*ratio1 );

   double ratio_out = (2.00*(Re+H) + h)/( sqrt1 * (1.00 + sqrt2) );
   return ratio_out;
}

double iono_loss_simple( double freq, double theta_deg )
{
   double loss_db_100MHz = gTypicalLossAt100MHz;
   double loss_db = loss_db_100MHz * power( (100.00/freq) , 2 );
   
   if( gIonoLossWithPathLength > 0 ){
      loss_db = loss_db * tau_g( theta_deg );
   }
   double loss_num = db2num( -loss_db );
   
   return loss_num;
}
            
double iono_loss( double freq, double theta_deg )
{
   double freq_hz = freq*1e6;
   double up = 2*M_PI*gPlasmaFreq*gPlasmaFreq*gCollisionFreq*delta_s(theta_deg);
   double down = c*(gCollisionFreq*gCollisionFreq + freq_hz*freq_hz);
                           
   double loss = 1.00 - up/down;
// double loss_db = 10.0*TMath::Log10(loss);
   return loss;   
}

// Evans & Hagfors 1968 : https://archive.org/stream/RadarAstronomy/EvansHagfors-RadarAstronomy#page/n143/mode/2up/search/ionosphere
// A. Datta paper 2014 
double iono_loss_adatta( double freq, double theta_deg )
{
   double freq_hz = freq*1e6;
   double TEC_D   = gTEC_total * gDLayerTEC_fraction;
   double n_e = TEC_D/gDeltaHeightD;
   double te_32 = power(gElectronTemp,1.5);
   double v_c = 3.65*(n_e/te_32)*(19.8 + log(te_32/freq_hz));
   double v_c_harish = gCollisionFreq;
   
   double loss_db = 1.16e-6 * ( (TEC_D*v_c)/(freq_hz*freq_hz) );
   double loss_num = 1.00/db2num(loss_db);

   double loss_num_harish = iono_loss(freq,theta_deg);
   
   if( gDebug>=5 ){
      printf("DEBUG(iono_abs) %.2f [MHz] - HARISH : v_c_H=%e [Hz], loss_H = %.4f , ABHI : v_c_A=%e [Hz], loss_A = %.4f\n",freq,v_c_harish,loss_num_harish,v_c,loss_num);
   }
   
   return loss_num;
}


// baed on / compare to plot_ionospheric_refraction_simple.C
double ThetaDepBailey( double zenith_distance_deg )
{
   double theta_deg = 90.00 - zenith_distance_deg;
   double elevation_of_typical_refraction_scale_deg = 90.00 - gTypicalRefractionZenithAngle;
   double Re = gRadiusEarth;
   double hm = 400e3; // height of max density in F-layer
   
   double theta_rad = theta_deg*(M_PI/180.00);
   double sin_theta = sin(theta_rad);   
   
   double theta_rad_45deg = elevation_of_typical_refraction_scale_deg*(M_PI/180.00);
   double sin_theta_45deg = sin(theta_rad_45deg);                  
   
   double geo_factor = power( (sin_theta*sin_theta + (2.0*hm/Re)) , -1.5 )*cos(theta_rad);
   double geo_factor_45deg = power( (sin_theta_45deg*sin_theta_45deg + (2.0*hm/Re)) , -1.5 )*cos(theta_rad_45deg);
      
   return (geo_factor / geo_factor_45deg);   
}

double refraction_angle( double freq, double theta_deg )
{
   if( gIncludeIonoRefraction == eIonoRefractionSimple ){    
      // simple 
      double angle_at_100MHz = gTypicalRefractionAngleArcMin;
      double angle_arcmin = angle_at_100MHz * power( (100.00/freq) , 2 );
      return (angle_arcmin/60.00); // return in deg
   }
   if( gIncludeIonoRefraction == eIonoRefractionSimpleWithThetaDep ){    
      // simple 
      double angle_at_100MHz = gTypicalRefractionAngleArcMin;
      double angle_arcmin = angle_at_100MHz * power( (100.00/freq) , 2 ) * ThetaDepBailey( theta_deg );
      return (angle_arcmin/60.00); // return in deg
   }

   double freq_hz = freq*1e6;
   double n_e = gTEC_total/gFlayerD; // e/m^3 
   double v_p = 9.5152 * sqrt(n_e);
   double theta_rad = deg2rad(theta_deg);
   double h_m = gFlayerHeight;
   double d   = gFlayerD;
   double bracket = sin(theta_rad)*sin(theta_rad) + (2.00*h_m)/gRadiusEarth;
   double bracket_to_power = 1.00/( bracket*sqrt(bracket) ); // bracket ^-3/2
   
   double delta_theta_rad = ((2.00*d)/(3.00*gRadiusEarth)) * mysqr(v_p/freq_hz) * (1.00 + h_m/gRadiusEarth) * bracket_to_power * cos(theta_rad);

   double ret = rad2deg(delta_theta_rad);
//   if( fabs(theta_deg) >= 90 ){
      // artificially ZERO shifts below the horizon (we have a ground screen so they shouldn't matter) :
//      ret = 0.00;
//   }
                        
   if( gDebug>=20 ){
      printf("DEBUG(refraction) at %.2f [Hz] : n_e = %e [e/m^3], v_p = %.2f [Hz], theta = %.2f [deg], h_m=%.3f [m], d=%.3f [m], bracket=%.8f, delta_thata = %.4f [arcsec]\n",
               freq_hz,n_e,v_p,theta_deg,h_m,d,bracket,ret*3600.00);
   }

   return ret;
}                                    

void test_refraction_angle()
{
   double theta_deg=0;
   
   while( theta_deg < 90 ){
      double angle_arcmin = refraction_angle( 100 , theta_deg )*60.00;
      printf("test_refraction_angle : %.4f %.8f [arcmin]\n",theta_deg,angle_arcmin);
      
      theta_deg += 1.00;
   }      
}

double tropo_refraction_angle( double theta_deg )
{
   double P = gAtmPressure; // milibars
   double a_deg = (90 - theta_deg);
   double a_rad = a_deg*(M_PI/180.00);
   double T_K = gT_amb;
            
   double R_deg = (0.00452* P) / ((T_K) * tan(a_rad));
   
   if( a_deg <= 15 ){
      R_deg  = P*(0.1594 + 0.0196*a_deg + 0.00002*a_deg*a_deg)/( (T_K)*(1 + 0.505*a_deg + 0.0845*a_deg*a_deg) );
   }
   
   return R_deg;               
}
