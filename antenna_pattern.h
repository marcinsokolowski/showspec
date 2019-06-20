#ifndef _ANTANNA_PATTERN_PARSER_H__
#define _ANTANNA_PATTERN_PARSER_H__

#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "cvalue_vector.h"
using namespace std;

#include "skymap_cache.h"

// 100 back to 500 
// #define MAX_THETA_COUNT 500
// #define MAX_PHI_COUNT   500
// TEST for Bradley :
// #define MAX_THETA_COUNT 9001
// #define MAX_PHI_COUNT 36001
extern int MAX_THETA_COUNT;
extern int MAX_PHI_COUNT;

#define DEFAULT_ANT_PATTERN_FILE "~/bighorns/software/analysis/showspec/EOR_bicone_rad_pat.dat"
#define DEFAULT_ANT_BEAM_WIDTH 50.00

#ifndef _NO_ROOT_
class TGraph2D;
#endif
class CAntPatternParser;

enum eAntennaPatternType { eAntPattDipolOverGndScreen=0, eAntPattIsotropic=1, eAntPattGaussian=2, eAntPattUserDefined=3, eMWA2016=2016 };

class CFreqPattern
{
/*   cFreqPattern( double _freq, int theta_samples, int phi_samples )
   : freq(_freq), m_ThetaSamples(theta_samples), m_PhiSamples(phi_samples)
   {
     ant_pattern = new double*[m_ThetaSamples];
     for(int t=0;t<m_ThetaSamples;t++){
        ant_pattern[t] = new double[m_PhiSamples];
        for(int p=0;p<m_PhiSamples;p++){
           ant_pattern[t][p] = 0;
        }
     }
   }*/
public :
   CFreqPattern( double _freq, int theta_samples, int phi_samples );
   CFreqPattern( const CFreqPattern& right );
   CFreqPattern& operator=( const CFreqPattern& right );
   ~CFreqPattern();   


   double freq; // MHz
   int m_ThetaSamples;
   int m_PhiSamples;
//   double ant_pattern[MAX_THETA_COUNT][MAX_PHI_COUNT];
   double** ant_pattern;
#ifndef _NO_ROOT_
   TGraph2D* pRootGraph; // mainly for interpolation purposes 
#else
   void* pRootGraph;
#endif

   double GetGain( double phi_deg, double theta_deg, double simul_max_phi, double simul_max_theta );   
/*   inline double GetGain( double phi_deg, double theta_deg, double simul_max_phi, double simul_max_theta )
   {
      double theta_step = (simul_max_theta/(m_ThetaSamples-1));
      double phi_step = (simul_max_phi/(m_PhiSamples-1));
         
      int phi_idx = (int)(phi_deg/phi_step);
      int theta_idx = (int)(theta_deg/theta_step);

      if( theta_deg > CAntPatternParser::m_SimulMaxTheta && CAntPatternParser::m_SimulMaxTheta<180 ){
         // if simulation only up to 90 deg -> gain below horizon is 0.00 , but if we are requesting theta > 180 -> error in code !!!
         return 0.00;
      }else{
         if( phi_idx <0 || phi_idx>=MAX_PHI_COUNT || theta_idx<0 || theta_idx>=MAX_THETA_COUNT ){
            printf("ERROR in CFreqPattern::GetGain angle index out of range of (%.2f,%.2f) [deg]\n",phi_deg,theta_deg);
            printf("phi_idx   = %d vs. allowed range (0-%d)\n",phi_idx,MAX_PHI_COUNT);
            printf("theta_idx = %d vs. allowed range (0-%d)\n",theta_idx,MAX_THETA_COUNT);
            printf("ERROR in code - there should not be request for theta and phi outside range [0,180] and [0,360] deg\n");
            exit(-1);
         }     
      }
                                      
      double gain = ant_pattern[theta_idx][phi_idx];
      return gain;                                             
   }*/
};

class CAntPatternParser
{
public :
  vector<CFreqPattern> m_FreqPatterns;
  static int m_bEnableInterpolation;

  // FILE PARSER FORMAT:
  static int m_InputFileFormat_GainColIdx;     
  static int m_InpitFileFormat_ThetaColIdx;
  static int m_InpitFileFormat_PhiColIdx;

  // simulation in FEKO (could be parsed automatically - but difficult at the moment)
  static double m_SimulMaxTheta;
  static double m_SimulMaxPhi;  
  static int m_AutoDetectTheta;
  static int m_BradleyVersion;
  static int m_ElevationInAntPattern;
  static int m_AnglesInRadians;
  static int m_GainInLinearScale;

  CAntPatternParser();
  ~CAntPatternParser();

  int FindFreq( double freq_mhz );
  
  double AutoCheckMaxTheta(const char* pattern_file);
  int Read(const char* pattern_file, int n_phi_samples=-1, int n_theta_samples=-1, double param_freq=-1);
  double GetGain(double freq,double phi_deg,double theta_deg);
  double GetGain(double freq,int freq_idx,double phi_deg,double theta_deg);

  double GetGain(double freq,CFreqPattern& pat_f1,CFreqPattern& pat_f2,double phi_deg,double theta_deg);
  CFreqPattern* GetFreqPattern( int freq_idx );
    
};

class CAntPatMap
{
public :
   CAntPatMap();
   CAntPatMap( const CAntPatMap& right );   
   ~CAntPatMap();
   CAntPatMap& operator=( const CAntPatMap& right );

   string m_szFreqFile;
   int m_Count;
   double* map;
};


class CAntPattern : public CAntPatternParser
{
protected :
   CAntPatternParser* m_AntennaPattern;
   CAntPatMap         m_CurrentMap;

public :
   static int          gMWAGridpoint;   
   static eAntennaPatternType m_AntennaPatternType;   
   double m_AntennaBeamWidth;
   static int m_bCacheON;
   static double m_RotFromAxis; // for Muresk : 20-30 deg rotation of South direction from the X-axis
                                // for NS 90 (Caiguna data)
   static string m_AntennaPatternFile;
   static double gSuppressAntPatternAboveZenAngle;
   static int gAntennaPatternFormula;
//   static double m_AntennaPatternIntegraldOMEGA;
   static CValueVector m_AntennaPatternIntegralsdOMEGA;
         

   CAntPattern();
   ~CAntPattern();
   void Init( vector<double>& ni_list, int n_phi_samples=-1, int n_theta_samples=-1, double param_freq=-1 );
   
   CAntPatternParser* GetAntParsedFile(){ return m_AntennaPattern; }
  
   double antenna_pattern_formula( double freq_mhz, double phi_deg, double theta_deg );   
   double antenna_pattern_dipol( double freq_mhz, double phi_deg, double theta_deg );

   double antenna_pattern( double freq_mhz,CFreqPattern* pPatternFreq , CFreqPattern* pPatternFreqNext,
                           double azim_deg, double zenith_dist_deg );
   double antenna_gain_from_azh( double freq_mhz, double azim_deg, double zenith_dist_deg, int bNormalised=1 );                           
   double antenna_gain_from_azz( double freq_mhz, double azim_deg, double zenith_dist_deg, int bNormalised=1 )
   {
      return antenna_gain_from_azh( freq_mhz, azim_deg, zenith_dist_deg, bNormalised );
   }
                                                            

   double* antenna_pattern( const char* infile, double freq_mhz, vector<cSkyInfo>& sky_intensities );
   
   // calculates antenna integral :
   double calculate_antenna_integral( double freq_mhz );
   
   // WARNING : should be 90-azim , but here azim is calculated from S towards W (libnova convention !!!)
   static inline double azim2phi( double azim_deg, double alpha_deg ){
      double phi_simul_deg = 270 - azim_deg + alpha_deg;
      if( phi_simul_deg >= 360 ){
         phi_simul_deg = phi_simul_deg - 360;
      }
      if( phi_simul_deg < 0 ){
         phi_simul_deg = 360 + phi_simul_deg;
      }

      return phi_simul_deg;                                             
   }
   
   // USER DEFINED PATTERN FORMULA :
   // double CAntPattern::antenna_pattern_formula( double freq_mhz, double phi_deg, double theta_deg )
   // static double (*m_UserDefinedPatternFormula)( double freq_mhz, double phi_deg, double theta_deg );
};

#endif
