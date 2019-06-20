#include <math.h>
#include <stdio.h>
#include "bg_globals.h"

// constants :
const double k_b=1.3806488E-23; // J/K
const double c=299792458.00;  // m/s

double log10(double x)
{
   return log(x)/log(10);
}

double mW2dbm(double in_mW)
{
   double dbm = 10.00*log10(in_mW);
   return dbm;   
}

double w2dbm( double in_w )
{
   return 10.00*log10( in_w*1000 );
}

double dbm2mW( double in_dbm )
{
   double l = in_dbm/10.00;
   double mW = power(10,l);

   return mW;
}

double dbm2mWperHz( double in_dbm, double rbw )
{
   double in_dbm_per_hz = in_dbm - 10.00*log10(rbw);
   double l = in_dbm_per_hz/10.00;
   double mW = power(10,l);
   double inW = mW/1000.00;
   double inK = inW/k_b;
               
   return mW;
}

double dbm2k( double in_dbm, double rbw )
{
   double inW = dbm2mWperHz(in_dbm,rbw)/1000.00;
   double inK = inW/k_b;
               
   return inK;
}

double k2dbm( double in_k )
{
   double in_W_per_hz = in_k * k_b;
   double in_mW_per_hz = in_W_per_hz*1000.00;
   double in_dBm_per_hz = 10.00 * log10( in_mW_per_hz );
   
   return in_dBm_per_hz;
}                                    

double kelvin2mWperHz( double in_k )
{
   double in_W_per_hz = in_k * k_b;
   double in_mW_per_hz = in_W_per_hz*1000.00;
   
   return in_mW_per_hz;
}                                    

double kelvin2dbmperHz( double in_k )
{
   return k2dbm( in_k ); 
}

double kelvin2mW( double in_k , double delta_freq_hz )
{
   double mW_per_hz = kelvin2mWperHz( in_k );
   double mW = mW_per_hz*delta_freq_hz;
   
   return mW;
}

double kelvin2dbm( double in_k , double delta_freq_hz )
{
   double mW = kelvin2mW( in_k , delta_freq_hz );
   double dbm = 10.00 * log10( mW );
   
   return dbm;
}



double num2db(double num)
{
   return (10.00*log10(num));   
}

// gain in given direction can be taken from FEKO
double calc_A_eff(double freq_mhz, double gain, double ant_efficiency ){
   double A_eff = (7161.97/(freq_mhz*freq_mhz))*(gain/ant_efficiency); // (7161.97 = (300)^2/(4 x PI) )
   return A_eff;
}

double jy2brigthnesstemp( double I_jy, double freq_mhz )
{
   double freq_hz = freq_mhz * 1e6;
   double jy_div_kb = 0.000724638;
   double jy_div_kb_times_c_div_2 = 32608695652173.913043478; // 10^-26 * c^2 / (2 * kb)
   double tb =                      32608695652173.913043478*(I_jy/(freq_hz*freq_hz));
   return tb;
}

double jy2kelvin( double S_jy /* in Jy */ , double A_eff /* in m^2 */ )
{
   // 0.5
   double t_ant = 0.000724297*S_jy*A_eff; // 0.000724297 = 10^-26/k_b = (10^-26)/(1.38x10^-23) = ...
   t_ant = 0.50*t_ant;
   return t_ant;
}

double kelvin2jy( double kelvin, double A_eff /* in m^2 */ )
{
   // 2.00
   double jy = (kelvin)/(0.000724297*A_eff); // 0.000724297 = 10^-26/k_b = (10^-26)/(1.38x10^-23) = ...
   jy = 2.00*jy;
// perhaps more familiar form is : = (1380.6488*kelvin)/A_eff
   return jy;
}

double jy2kelvin_gain(double S_jy /* in Jy */ , double freq_mhz, double gain, double ant_efficiency )
{
   double A_eff = calc_A_eff( freq_mhz, gain, ant_efficiency );
   printf("A_eff = %.2f [m^2]\n",A_eff);
   double t_ant = jy2kelvin( S_jy , A_eff );
   
   return t_ant;       
}

double kelvin2jy_gain( double kelvin, double freq_mhz, double gain, double ant_efficiency )
{
   double A_eff = calc_A_eff( freq_mhz, gain, ant_efficiency );
   double jy = kelvin2jy( kelvin, A_eff );
   
   return jy;
}
