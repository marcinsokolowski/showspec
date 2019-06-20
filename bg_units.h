#ifndef _BG_UNITS_H__
#define _BG_UNITS_H__

// extern const double k_b; // J/K
// extern const double c;  // m/s


double mW2dbm(double in_mW);

double dbm2mW( double in_dbm );
double dbm2mWperHz( double in_dbm, double rbw );

double dbm2k( double in_dbm, double rbw );
double k2dbm( double in_k ); // Kelvin -> dBm/Hz
double num2db( double num );
double w2dbm( double in_w ) ;

// K-> mW/Hz, dBm/Hz:
double kelvin2mWperHz( double in_k );
double kelvin2dbmperHz( double in_k );
double kelvin2mW( double in_k , double delta_freq_hz );
double kelvin2dbm( double in_k , double delta_freq_hz );

// Jy -> K :
double calc_A_eff(double freq_mhz, double gain=1.00, double ant_efficiency=1.00 );
double jy2kelvin( double S_jy /* in Jy */ , double A_eff /* in m^2 */ );
double kelvin2jy( double kelvin, double A_eff /* in m^2 */ );
double jy2kelvin_gain(double S_jy /* in Jy */ , double freq_mhz, double gain=1.00, double ant_efficiency=1.00 );

// Jy -> brightness temperature (not antenna temperature ! )
double jy2brigthnesstemp( double I_jy, double freq_mhz );

#endif
