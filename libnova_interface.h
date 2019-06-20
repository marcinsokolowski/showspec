#ifndef _LIBNOVA_INTERFACE_H__
#define _LIBNOVA_INTERFACE_H__

#include <libnova/solar.h>

// global variables 
extern double gObsLongInDeg;
extern double gObsLatInDeg;

// SUN :
// out_rise_ux - next sunrise
// out_set_ux  - next sunset
void get_sun_info( time_t unix_time, double geo_long_deg, double geo_lat_deg, time_t& out_rise_ux, time_t& out_set_ux, time_t& out_transit_ux, int bVerb=1, double horizon=LN_SOLAR_STANDART_HORIZON );
void get_sun_azh( time_t unix_time, double geo_long_deg, double geo_lat_deg, double&  out_az, double& out_alt );
void get_sun_all( time_t unix_time, double geo_long_deg, double geo_lat_deg,
                  double& out_ra, double& out_dec, double&  out_az, double& out_alt );
void get_sun_rise_set( time_t unix_time, double geo_long_deg, double geo_lat_deg, time_t& out_sun_rise_ux, time_t& out_sun_set_ux, double horizon=LN_SOLAR_STANDART_HORIZON, int step=1 );
                  

// GALAXY : 
void get_galaxy_info( time_t unix_time, double geo_long_deg, double geo_lat_deg,
                   time_t& out_rise_ux, time_t& out_set_ux, time_t& out_transit_ux );

void get_galaxy_azh( time_t unix_time, double geo_long_deg, double geo_lat_deg, double&  out_az, double& out_alt );                   

void get_galaxy_up_down( time_t unix_time, double geo_long_deg, double geo_lat_deg,
                         time_t& max_ux, double& h_max,
                         time_t& min_ux, double& h_min );

// MOON : 
void get_moon_info( time_t unix_time, double geo_long_deg, double geo_lat_deg,
                   time_t& out_rise_ux, time_t& out_set_ux, time_t& out_transit_ux, 
                   double& phase_in_transit, double& phase_angle_in_transit, double& moon_radius_deg );

void get_moon_info( time_t unix_time, double geo_long_deg, double geo_lat_deg,
                    double& out_ra, double& out_dec,
                    double& out_az, double& out_alt,
                    double& phase, double& phase_angle, double& moon_radius_deg );

// TIME / LST ect :
double get_local_sidereal_time(double uxtime_d);
double get_local_sidereal_time(double uxtime_d,double& jd_out);
double get_local_sidereal_time(double uxtime_d,double geo_long_deg,double& jd_out); 
double get_jd_plus_ndays_sidalligned(double start_ux_d,double n_days);
double ux2jd( time_t unix_time );
double uxd2jd( double unix_time );


// general functions :
double get_ang_distance_deg( double ra1, double dec1, double ra2, double dec2 );

void gal2hor( double glon_deg, double glat_deg, double geo_long_deg, double geo_lat_deg, time_t unix_time,
              double& out_azim, double& out_alt, double& out_ra, double& out_dec );
void gal2hor_fromJD( double glon_deg, double glat_deg, double geo_long_deg, double geo_lat_deg, double JD,
              double& out_azim, double& out_alt, double& out_ra, double& out_dec );              
              
void hor2gal_fromJD( double azim, double alt, double geo_long_deg, double geo_lat_deg, double JD, double& out_glon_deg, double& out_glat_deg, double& out_ra, double& out_dec );              

void radec2azh( double ra, double dec, time_t unix_time, double& out_azim, double& out_alt );
              
void get_ymd_hms_ut( time_t ut_time, int& year, int& month, int& day, int& hour, int& minute, double& sec );

void radec2azh( double ra, double dec, time_t unix_time, double geo_long_deg, double geo_lat_deg, double& out_az, double& out_alt );                  
void azh2radec( double az, double alt, time_t unix_time, double geo_long_deg, double geo_lat_deg, double& out_ra, double& out_dec );
double hour_angle( double ra, time_t unix_time );

#endif
