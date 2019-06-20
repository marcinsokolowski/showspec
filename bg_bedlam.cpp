#include <math.h>
#include "bg_bedlam.h"

int CBedlamSpectrometer::m_SizeX = 4096;
double CBedlamSpectrometer::start_freq = 0.00;
double CBedlamSpectrometer::stop_freq  = 480.00;
int CBedlamSpectrometer::m_nSpectralModelPolyN=10;

double CBedlamSpectrometer::gBedlamQualityLimitDBM_BAD=-37; // see document logbook/201403/tone_test_vs_ORBCOMM.pdf (-100dBm + 63dB gain ~= -37dBm at the BEDLAM input see also bedlam_calibration_tone.odt/pdf )
// The more strict limit would be -110dBm which is super safe then -110 + 63 = -47dBm :
double CBedlamSpectrometer::gBedlamQualityLimitDBM_START=-47; // -110 + 63 = -47 dBm (problems start somewhere in between -47 - -37 dBm )



CBedlamSpectrometer::CBedlamSpectrometer()
{}

// cd /media/BighornsData/labtest/tone_test/20140311/ZVL_comparison_100MHz_60sec/ZVL_100mhz_avg10000
// modelspec_poly acc32000_ch0_20140311_023648_000005.ZVL_CALIB_FIT
// MORE DETAILS in tone_test__vs_ORBCOMM.odt 
// this was done by :
//    measuring power with BEDLAM  -> P_bedlam(freq)
//    measuring power with the ZBL -> P_zvl(freq) -> P_zvl_mW(freq)
//  and -> calib(freq) = P_bedlam(freq) / P_zvl_mW(freq)
//  later polynomial (as below) was fited to calib(freq) and later is used to calibrated any BEDLAM power to dBm by :
//  P_bedlam_dBm(freq) = P_bedlam(freq)/calib(freq)
// this means the method should be robust - no matter how the system is changed it should always relay only on the P_bedlam !!!  
double CBedlamSpectrometer::spectrum_response_model( double freq )
{
   // filters suppress any power outsie ~20-380 MHz range :
   if( freq<=20 || freq>=380 )
      return 1e20;
            

   double par[50];
   for(int i=0;i<50;i++){
      par[i] = 0.00;
   }

/* par[0]                        = -3.52364e+12   ; // +/-    1.041e+12
   par[1]                        =  4.13935e+11   ; // +/-    6.51364e+10
   par[2]                        = -1.02092e+10   ; // +/-    1.7242e+09 
   par[3]                        =  1.23619e+08   ; // +/-    2.53996e+07
   par[4]                        =      -775448   ; // +/-    230140     
   par[5]                        =      2199.82   ; // +/-    1334.28    
   par[6]                        =     0.427662   ; // +/-    4.96585    
   par[7]                        =   -0.0193973   ; // +/-    0.0114757  
   par[8]                        =  4.62558e-05   ; // +/-    1.49857e-05
   par[9]                        = -3.58946e-08   ; // +/-    8.44799e-09 */

// cd /media/BighornsData/labtest/tone_test/20140311/ZVL_comparison_100MHz_60sec/ZVL_100mhz_avg10000
// modelspec_poly acc32000_ch0_20140311_023648_000005.ZVL_CALIB_FIT
par[0] = 4722515362706.85351562500000000000000000000000000000000000000000;
par[1] = -95027966093.02354431152343750000000000000000000000000000000000; 
par[2] = 2858485415.97554826736450195312500000000000000000000000000000;   
par[3] = -60083923.19348197430372238159179687500000000000000000000000;    
par[4] = 788743.02736626658588647842407226562500000000000000000000;       
par[5] = -6199.39397118529541330644860863685607910156250000000000;        
par[6] = 28.97252980372121555774356238543987274169921875000000;           
par[7] = -0.07877571963827538492619595444921287707984447479248;           
par[8] = 0.00011499439127208088760705162467701256900909356773;            
par[9] = -0.00000006966022305926903079173408313284898696338132;           


   double ret=0.00;
   int n = m_nSpectralModelPolyN;

   double pow=1.00;
   for(int i=0;i<n;i++){
      double next = pow*par[i];
      
      ret += next;
      pow = pow*freq;
   }

   return ret;
}

double CBedlamSpectrometer::ch2freq( int ch )
{
   double ch_res = (stop_freq-start_freq)/(m_SizeX-1);

   return start_freq + ch*ch_res;   
}

double CBedlamSpectrometer::freq2ch( double freq_mhz )
{
   int ch = ( (freq_mhz-start_freq)/(stop_freq-start_freq) ) * m_SizeX;
   return ch;   
}
          
double CBedlamSpectrometer::power2mW( double freq_mhz, double bedlam_power )
{
   double calib_const = spectrum_response_model( freq_mhz );
   double power_mW = bedlam_power / calib_const;
   
   return power_mW;
}

double CBedlamSpectrometer::power2dbm( double freq_mhz, double bedlam_power )
{
   double power_mW = power2mW( freq_mhz, bedlam_power );
   double power_dBm = 10.00 * log10( power_mW );
   
   return power_dBm;
}
                
double CBedlamSpectrometer::mW2bedlampower( double freq_mhz, double power_mW )
{
   double calib_const = spectrum_response_model( freq_mhz );
   double bedlam_power = power_mW * calib_const;
      
   return bedlam_power;         
}

double CBedlamSpectrometer::dbm2bedlampower( double freq_mhz, double power_dbm )
{
   double in_mW = exp(log10(10.00)*power_dbm/10.00);
   double bedlam_power = mW2bedlampower( freq_mhz, in_mW );
      
   return bedlam_power;
}

         