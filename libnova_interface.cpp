// LIBNOVA HELP :
//  local   : file:///opt/caastro/ext/dload/libnova/doc/html/index.html
//  webpage : 
#include <libnova/utility.h>
#include <libnova/solar.h>
#include <libnova/julian_day.h>
#include <libnova/rise_set.h>
#include <libnova/transform.h>
#include <libnova/lunar.h>
#include <libnova/sidereal_time.h>
#include <libnova/julian_day.h>
#include <libnova/angular_separation.h>

#include <stdio.h>
#include <stdlib.h>

#include "libnova_interface.h"
#include "bg_globals.h"

#include <string>

using namespace std;

const double gGalaxyCenterRA=266.41683333;
const double gGalaxyCenterDEC=-29.00780556;

// double geo_long=116.68609722; // default Muresk
// double geo_lat=-31.74625000;  // default Muresk
   

void get_ymd_hms_ut( time_t ut_time, int& year, int& month, int& day,
                  int& hour, int& minute, double& sec )
{
	struct tm gmtm;
	gmtime_r( &ut_time , &gmtm );
	year = (gmtm.tm_year-100)+2000;
   month = (gmtm.tm_mon+1);
   day = gmtm.tm_mday;
	hour = gmtm.tm_hour;
	minute = gmtm.tm_min;
	sec = gmtm.tm_sec;
}


/*string get_gmtime_string( time_t ut_time )
{
   struct tm gmtime_tm;
   string szRet;
   if(gmtime_r( &ut_time, &gmtime_tm )){
      char tempstring[64];
              
      // bug ??? first %.2u -> %.4u ???
      sprintf(tempstring,"%.2u%.2u%.2u_%.2u%.2u%.2u",
                          gmtime_tm.tm_year+1900,(gmtime_tm.tm_mon+1),gmtime_tm.tm_mday,
                          gmtime_tm.tm_hour,gmtime_tm.tm_min,gmtime_tm.tm_sec);
      szRet = tempstring;
   }  
   return szRet;
}*/
                                                                                             

// galatic center : RA 17h45m40.04s, Dec -29° 00' 28.1"
// 266.41683333 deg , -29.00780556 deg 
void get_galaxy_info( time_t unix_time, double geo_long_deg, double geo_lat_deg,
                   time_t& out_rise_ux, time_t& out_set_ux, time_t& out_transit_ux )
{
   struct ln_lnlat_posn observer;
   double JD;
   struct ln_date date;   
   
   get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
   JD = ln_get_julian_day (&date);
   printf("Getting galaxy info for unix_time=%d -> jd=%.8f\n",(int)unix_time,JD);

   // position of observatory :
   observer.lng = geo_long_deg;
   observer.lat = geo_lat_deg;
   
   struct ln_equ_posn galaxy_coord;
   galaxy_coord.ra = gGalaxyCenterRA; 
   galaxy_coord.dec = gGalaxyCenterDEC;
                   
   struct ln_rst_time galaxy_info;
   ln_get_object_rst (JD, &observer, &galaxy_coord, &galaxy_info);

   ln_get_timet_from_julian( galaxy_info.rise , &out_rise_ux );
   ln_get_timet_from_julian( galaxy_info.set , &out_set_ux );
   ln_get_timet_from_julian( galaxy_info.transit , &out_transit_ux );

   double max_alt,min_alt;
   time_t max_alt_ux, min_alt_ux;
   get_galaxy_up_down( unix_time, geo_long_deg, geo_lat_deg, max_alt_ux, max_alt, min_alt_ux, min_alt );

   string szRise,szTransit,szSet; 
   szRise = get_localtime_string(out_rise_ux);
   szTransit = get_localtime_string(out_transit_ux);
   szSet = get_localtime_string(out_set_ux);            
   
   printf("Galaxy rise    = %.8f JD = %d UX = %s LOCAL\n",galaxy_info.rise,(int)out_rise_ux,szRise.c_str());
   printf("Galaxy transit = %.8f JD = %d UX = %s LOCAL\n",galaxy_info.transit,(int)out_transit_ux,szTransit.c_str());
   printf("Galaxy set     = %.8f JD = %d UX = %s LOCAL\n",galaxy_info.set,(int)out_set_ux,szSet.c_str());
   printf("Galaxy min alt = %.2f [deg] at ux = %d ( max alt = %.2f [deg] at %d - should be same as transit )\n",min_alt,(int)min_alt_ux,max_alt,(int)max_alt_ux);   
}

void get_galaxy_azh( time_t unix_time, double geo_long_deg, double geo_lat_deg, double&  out_az, double& out_alt )
{
   struct ln_lnlat_posn observer;
   double JD;
   struct ln_date date;   
   
   get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
   JD = ln_get_julian_day (&date);
//   printf("Getting galaxy info for unix_time=%d -> jd=%.8f\n",(int)unix_time,JD);

   // position of observatory :
   observer.lng = geo_long_deg;
   observer.lat = geo_lat_deg;
   
   struct ln_equ_posn galaxy_coord;
   galaxy_coord.ra = gGalaxyCenterRA; 
   galaxy_coord.dec = gGalaxyCenterDEC;

   struct ln_hrz_posn hrz_pos;
   ln_get_hrz_from_equ( &galaxy_coord, &observer, JD, &hrz_pos );
   out_az  = hrz_pos.az;
   out_alt = hrz_pos.alt;
}

void get_sun_azh( time_t unix_time, double geo_long_deg, double geo_lat_deg, double&  out_az, double& out_alt )
{
   struct ln_lnlat_posn observer;
   double JD;
   struct ln_date date;   
   
   get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
   JD = ln_get_julian_day (&date);
//   printf("Getting galaxy info for unix_time=%d -> jd=%.8f\n",(int)unix_time,JD);

   ln_equ_posn sun_coord;
   ln_get_solar_equ_coords(JD,&sun_coord);

   // position of observatory :
   observer.lng = geo_long_deg;
   observer.lat = geo_lat_deg;
   
   struct ln_hrz_posn hrz_pos;
   ln_get_hrz_from_equ( &sun_coord, &observer, JD, &hrz_pos );
   out_az  = hrz_pos.az;
   out_alt = hrz_pos.alt;
}

void get_sun_rise_set( time_t unix_time, double geo_long_deg, double geo_lat_deg, time_t& out_sun_rise_ux, time_t& out_sun_set_ux,
                       double horizon, int step )   
{
   out_sun_rise_ux = 0;
   out_sun_set_ux = 0;
   
   double curr_sun_az, curr_sun_alt;
   get_sun_azh( unix_time, geo_long_deg, geo_lat_deg, curr_sun_az, curr_sun_alt );
   
   // looking BACKWARD for sign change in ALT 
   time_t ux = unix_time;
   double sun_alt = curr_sun_alt,sun_az=curr_sun_az;
   while (ux >= (unix_time-86400) && (sun_alt*curr_sun_alt)>0){
      ux = ux - step;
      get_sun_azh( ux, geo_long_deg, geo_lat_deg, sun_az, sun_alt );
   }
   if( sun_alt<0 && curr_sun_alt>0 ){
      // sun-rise identified, but go down until sun_alt is below required elevation :
      while (sun_alt<0 && sun_alt>horizon ){
         ux = ux - step;
         get_sun_azh( ux, geo_long_deg, geo_lat_deg, sun_az, sun_alt );
      }
      out_sun_rise_ux = ux;
   }else{
      // curr_sun_alt<0 && sun_alt>0  -> sun-set identified 
      // go forward until below SunElev < horizon :
      while ( sun_alt>horizon && ux<=(unix_time+86400) ){
         ux = ux + step;
         get_sun_azh( ux, geo_long_deg, geo_lat_deg, sun_az, sun_alt );
      }
      out_sun_set_ux = ux;      
   }
   
/*   if( (sun_alt*curr_sun_alt) <= 0 ){
     // sign changes :
     if( curr_sun_alt > 0 ){
        // sun rise detected :        
        out_sun_rise_ux = ux;
     }else{
        if( curr_sun_alt < 0 ){
           out_sun_set_ux = ux;
        }else{
          // equal -> previous step decides :          
          if( sun_alt > 0 ){
             out_sun_set_ux = ux;
          }else{
              out_sun_rise_ux = ux;
          }
        }
     }
   }*/
   
   // looking FORWARD for sign change in ALT
   ux = unix_time;sun_alt = curr_sun_alt,sun_az=curr_sun_az;
   while (ux <= (unix_time+86400) && (sun_alt*curr_sun_alt)>0){
      ux = ux + step;
      get_sun_azh( ux, geo_long_deg, geo_lat_deg, sun_az, sun_alt );
   }
   if( sun_alt<0 && curr_sun_alt>0 ){
      // sun-set identified, but go down until sun_alt is below required elevation :
      while( sun_alt<0 && sun_alt>horizon ){
         ux = ux + step;
         get_sun_azh( ux, geo_long_deg, geo_lat_deg, sun_az, sun_alt );
      }
      out_sun_set_ux = ux;
   }else{
      // curr_sun_alt<0 && sun_alt>0 -> sun rise
      // go back until hit SunElev < horzion
      while( sun_alt>horizon && ux>=(unix_time-86400) ){
         ux = ux - step;
         get_sun_azh( ux, geo_long_deg, geo_lat_deg, sun_az, sun_alt );
      }
      out_sun_rise_ux = ux;      
   }

/*   if( (sun_alt*curr_sun_alt) <= 0 ){
     // sign changes :
     if( curr_sun_alt > 0 ){
        // sun rise detected :        
        out_sun_set_ux = ux;
     }else{
        if( curr_sun_alt < 0 ){
           out_sun_rise_ux = ux;
        }else{
          // equal -> next step decides :          
          if( sun_alt > 0 ){
             out_sun_rise_ux = ux;
          }else{
              out_sun_set_ux = ux;
          }
        }
     }
   }*/

   if( out_sun_rise_ux == 0 || out_sun_set_ux == 0 ){
      printf("ERROR in code : could not determine SunSet or SunRise for unix_time = %d , geo_long = %.8f , geo_lat = %.8f\n",(int)unix_time,geo_long,geo_lat);
      exit(-1);
   }
}


void get_sun_all( time_t unix_time, double geo_long_deg, double geo_lat_deg, 
                  double& out_ra, double& out_dec, double&  out_az, double& out_alt )
{
   struct ln_lnlat_posn observer;
   double JD;
   struct ln_date date;   
   
   get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
   JD = ln_get_julian_day (&date);
//   printf("Getting galaxy info for unix_time=%d -> jd=%.8f\n",(int)unix_time,JD);

   ln_equ_posn sun_coord;
   ln_get_solar_equ_coords(JD,&sun_coord);
   
   out_ra = sun_coord.ra;
   out_dec = sun_coord.dec;

   // position of observatory :
   observer.lng = geo_long_deg;
   observer.lat = geo_lat_deg;
   
   struct ln_hrz_posn hrz_pos;
   ln_get_hrz_from_equ( &sun_coord, &observer, JD, &hrz_pos );
   out_az  = hrz_pos.az;
   out_alt = hrz_pos.alt;
}



void get_galaxy_up_down( time_t unix_time, double geo_long_deg, double geo_lat_deg,
                         time_t& max_ux, double& h_max,
                         time_t& min_ux, double& h_min )
{
   int dt=60;
   h_max = -10000.00;
   h_min = 10000.00;
   max_ux = 0;
   min_ux = 0;
   
   time_t t = unix_time;
   while( t < (unix_time+24*3600) ){
      double az,h;
      
      get_galaxy_azh( t, geo_long_deg, geo_lat_deg, az ,h );
      if( h > h_max ){
         h_max = h;
         max_ux = t;
      }
      if( h < h_min ){
         h_min = h;
         min_ux = t;
      }
      
      t += dt;
   }    
}                         

void get_sun_info( time_t unix_time, double geo_long_deg, double geo_lat_deg,
                   time_t& out_rise_ux, time_t& out_set_ux, time_t& out_transit_ux,
                   int bVerb, double horizon )
{
   struct ln_lnlat_posn observer;
   double JD;
   struct ln_date date;   
   
   get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
   JD = ln_get_julian_day (&date);
   if( bVerb > 0 ){
      printf("Getting solar info for unix_time=%d -> jd=%.8f\n",(int)unix_time,JD);
   }

   // position of observatory :
   observer.lng = geo_long_deg;
   observer.lat = geo_lat_deg;
             
   struct ln_rst_time sun_info;
//   ln_get_solar_rst( JD, &observer, &sun_info );
   ln_get_solar_rst_horizon( JD, &observer, horizon, &sun_info );

   ln_get_timet_from_julian( sun_info.rise , &out_rise_ux );
   ln_get_timet_from_julian( sun_info.set , &out_set_ux );
   ln_get_timet_from_julian( sun_info.transit , &out_transit_ux );

   string szRise,szTransit,szSet;
   szRise = get_localtime_string(out_rise_ux);
   szTransit = get_localtime_string(out_transit_ux);
   szSet = get_localtime_string(out_set_ux);                        
   
   if( bVerb > 0 ){
      double jd;
      double lst_sunrise = get_local_sidereal_time( (double)out_rise_ux, jd );      
      double lst_transit = get_local_sidereal_time( (double)out_transit_ux, jd );      
      double lst_sunset = get_local_sidereal_time( (double)out_set_ux, jd );      

      double lt_sunrise = get_local_hour_decimal( out_rise_ux );
      double lt_transit = get_local_hour_decimal( out_transit_ux );
      double lt_sunset = get_local_hour_decimal( out_set_ux );
            
   
      printf("Solar rise    = %.8f JD = %d UX = %s LOCAL (%.4f) = %.4f (%02d:%02d) LST, LT-LST = %.4f\n",sun_info.rise,(int)out_rise_ux,szRise.c_str(),lt_sunrise,lst_sunrise,((int)lst_sunrise),(int)((lst_sunrise-(int)lst_sunrise)*60),(lt_sunrise-lst_sunrise));
      printf("Solar transit = %.8f JD = %d UX = %s LOCAL (%.4f) = %.4f (%02d:%02d) LST, LT-LST = %.4f\n",sun_info.transit,(int)out_transit_ux,szTransit.c_str(),lt_transit,lst_transit,((int)lst_transit),(int)((lst_transit-(int)lst_transit)*60),(lt_transit-lst_transit));
      printf("Solar set     = %.8f JD = %d UX = %s LOCAL (%.4f) = %.4f (%02d:%02d) LST, LT-LST = %.4f\n",sun_info.set,(int)out_set_ux,szSet.c_str(),lt_sunset,lst_sunset,((int)lst_sunset),(int)((lst_sunset-(int)lst_sunset)*60),(lt_sunset-lst_sunset));
      printf("Sunrise (AWST) & Sunrise (LST) & Sunset (AWST) & Sunset (LST) & LST - Local time = %.2f & %.2f & %.2f & %.2f & %.2f \n",lt_sunrise,lst_sunrise,lt_sunset,lst_sunset,(lst_transit-lt_transit));
   }
}

void get_moon_info( time_t unix_time, double geo_long_deg, double geo_lat_deg,
                   time_t& out_rise_ux, time_t& out_set_ux, time_t& out_transit_ux, 
                   double& phase_in_transit , double& phase_angle_in_transit, double& moon_radius_deg )
{
   struct ln_lnlat_posn observer;
   double JD;
   struct ln_date date;   
   
   get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
   JD = ln_get_julian_day (&date);
   printf("Getting solar info for unix_time=%d -> jd=%.8f\n",(int)unix_time,JD);

   // position of observatory :
   observer.lng = geo_long_deg;
   observer.lat = geo_lat_deg;
             
   struct ln_rst_time moon_info;
   ln_get_lunar_rst( JD, &observer, &moon_info );

   phase_in_transit =  ln_get_lunar_disk( moon_info.transit );
   phase_angle_in_transit = ln_get_lunar_phase( moon_info.transit );

   ln_get_timet_from_julian( moon_info.rise , &out_rise_ux );
   ln_get_timet_from_julian( moon_info.set , &out_set_ux );
   ln_get_timet_from_julian( moon_info.transit , &out_transit_ux );

   moon_radius_deg = ln_get_lunar_sdiam( JD )/3600.00;
   
   string szRise,szTransit,szSet; 
   szRise = get_localtime_string(out_rise_ux);
   szTransit = get_localtime_string(out_transit_ux);
   szSet = get_localtime_string(out_set_ux);                    
   
   printf("Moon rise    = %.8f JD = %d UX = %s LOCAL\n",moon_info.rise,(int)out_rise_ux,szRise.c_str());
   printf("Moon transit = %.8f JD = %d UX = %s LOCAL ( phase = %.2f %% , phase_angle = %.2f [deg] )\n",moon_info.transit,(int)out_transit_ux,szTransit.c_str(),phase_in_transit*100.00,phase_angle_in_transit);
   printf("Moon set     = %.8f JD = %d UX = %s LOCAL\n",moon_info.set,(int)out_set_ux,szSet.c_str());
   printf("Moon angular size = %.4f [deg]\n",moon_radius_deg);
}

void get_moon_info( time_t unix_time, double geo_long_deg, double geo_lat_deg, 
                   double& out_ra, double& out_dec,
                   double& out_az, double& out_alt, 
                   double& phase, double& phase_angle, double& moon_radius_deg )
{
   struct ln_lnlat_posn observer;
   double JD;
   struct ln_date date;   
   
   get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
   JD = ln_get_julian_day (&date);
//   printf("Getting galaxy info for unix_time=%d -> jd=%.8f\n",(int)unix_time,JD);

   // position of observatory :
   observer.lng = geo_long_deg;
   observer.lat = geo_lat_deg;
   
   struct ln_equ_posn moon_coord;
   ln_get_lunar_equ_coords( JD, &moon_coord );   
   out_ra = moon_coord.ra;
   out_dec = moon_coord.dec;

   struct ln_hrz_posn hrz_pos;
   ln_get_hrz_from_equ( &moon_coord, &observer, JD, &hrz_pos );
   out_az  = hrz_pos.az;
   out_alt = hrz_pos.alt;

   phase =  ln_get_lunar_disk( JD );
   phase_angle = ln_get_lunar_phase( JD );    
   moon_radius_deg = ( ln_get_lunar_sdiam( JD ) / 3600.0 ) / 2.00;
   
   printf("Moon angular radius = %.4f [deg]\n",moon_radius_deg);
}                   

double ux2jd( time_t unix_time )
{
   double JD;
   struct ln_date date;

   /* UT date and time */
   get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
   JD = ln_get_julian_day (&date);

   return JD;                      
}

double uxd2jd( double unix_time )
{
   double JD;
   struct ln_date date;
   time_t unix_time_int = (time_t)unix_time;

   /* UT date and time */
   get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
   JD = ln_get_julian_day (&date);

   JD += double(unix_time-unix_time_int)/86400.00;

   return JD;                      
}


void gal2hor( double glon_deg, double glat_deg, double geo_long_deg, double geo_lat_deg, time_t unix_time, 
              double& out_azim, double& out_alt, double& out_ra, double& out_dec )
{
     struct ln_lnlat_posn observer;
     double JD;
     struct ln_date date;
	  struct ln_gal_posn gal_pos;
	  struct ln_equ_posn equ_pos;
	  struct ln_hrz_posn azim_alt_pos;

	  // position of observatory :
	  observer.lng = geo_long_deg;
	  observer.lat = geo_lat_deg;
        
	  // assign input values 
	  gal_pos.l = glon_deg;
	  gal_pos.b = glat_deg;

	  // calculate equatorial coordinates :
//	  ln_get_equ_from_gal( &gal_pos, &equ_pos );
     ln_get_equ2000_from_gal( &gal_pos, &equ_pos );
	  
	  out_ra = equ_pos.ra;
	  out_dec = equ_pos.dec;

     /* UT date and time */
     get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
     JD = ln_get_julian_day (&date);

	  // calculate horizontal coordinates :
	  ln_get_hrz_from_equ(&equ_pos, &observer, JD, &azim_alt_pos );

	  out_azim = azim_alt_pos.az;
	  out_alt = azim_alt_pos.alt;
}

void gal2hor_fromJD( double glon_deg, double glat_deg, double geo_long_deg, double geo_lat_deg, double JD, 
              double& out_azim, double& out_alt, double& out_ra, double& out_dec )
{
     struct ln_lnlat_posn observer;
	  struct ln_gal_posn gal_pos;
	  struct ln_equ_posn equ_pos;
	  struct ln_hrz_posn azim_alt_pos;

	  // position of observatory :
	  observer.lng = geo_long_deg;
	  observer.lat = geo_lat_deg;
        
	  // assign input values 
	  gal_pos.l = glon_deg;
	  gal_pos.b = glat_deg;

	  // calculate equatorial coordinates :
//	  ln_get_equ_from_gal( &gal_pos, &equ_pos );
     ln_get_equ2000_from_gal( &gal_pos, &equ_pos );
	  
	  out_ra = equ_pos.ra;
	  out_dec = equ_pos.dec;

	  // calculate horizontal coordinates :
	  ln_get_hrz_from_equ(&equ_pos, &observer, JD, &azim_alt_pos );

	  out_azim = azim_alt_pos.az;
	  out_alt = azim_alt_pos.alt;
}

void hor2gal_fromJD( double azim, double alt, double geo_long_deg, double geo_lat_deg, double JD, double& out_glon_deg, double& out_glat_deg, double& out_ra, double& out_dec )
{
     struct ln_lnlat_posn observer;
	  struct ln_gal_posn gal_pos;
	  struct ln_equ_posn equ_pos;
	  struct ln_hrz_posn azim_alt_pos;

	  // position of observatory :
	  observer.lng = geo_long_deg;
	  observer.lat = geo_lat_deg;
        
	  // assign input values 
	  azim_alt_pos.az = azim;
	  azim_alt_pos.alt = alt;

	  // calculate equatorial coordinates from azimuthal :
     ln_get_equ_from_hrz( &azim_alt_pos, &observer, JD, &equ_pos );
	  
	  out_ra = equ_pos.ra;
	  out_dec = equ_pos.dec;

	  // calculate galactic coordinates :
//	  ln_get_hrz_from_equ(&equ_pos, &observer, JD, &azim_alt_pos );
     ln_get_gal_from_equ2000( &equ_pos, &gal_pos );

     out_glon_deg = gal_pos.l;
     out_glat_deg = gal_pos.b;
}



double cut_to_range(double& sid_local_h)
{
   if( sid_local_h > 24.00 ){
      sid_local_h = sid_local_h-24.00;
   }
   
   if( sid_local_h < 0 ){
      sid_local_h = sid_local_h+24.00;
   }
   
   return sid_local_h;
}

double get_local_sidereal_time(double uxtime_d,double geo_long_deg,double& jd_out)
{
  time_t uxtime = (time_t)uxtime_d;
  jd_out = ln_get_julian_from_timet( &uxtime );
  if( gBGPrintfLevel >= 3 ){
     printf("DEBUG : jd(%d = floor(%.8f)) = %.8f\n",(int)uxtime,uxtime_d,jd_out);
  }
  jd_out += (uxtime_d - uxtime)/(24.00*3600.00);
  if( gBGPrintfLevel >= 3 ){
     printf("DEBUG : jd_prim = %.8f\n",jd_out);
  }
  

  // which to use ???
  double sid_greenwich_h = ln_get_apparent_sidereal_time(jd_out);
//  double sid_greenwich_h = ln_get_mean_sidereal_time (jd_out);
  if( gBGPrintfLevel >= 3 ){
     printf("DEBUG : sid_greenwich = %.8f [h] = %.8f [deg]\n",sid_greenwich_h,sid_greenwich_h*15.00);
  }
  
  double sid_local_h = sid_greenwich_h + geo_long_deg/15.00;
  if( gBGPrintfLevel >= 3 ){
     printf("DEBUG : sid_local = %.8f + %.8f = %.8f [h] = %.8f [deg]\n",sid_greenwich_h,geo_long_deg/15.00,sid_local_h,sid_local_h*15.00);
  }

  cut_to_range(sid_local_h);
     
  return sid_local_h;  
}
    
double get_local_sidereal_time(double uxtime_d,double& jd_out)
{    
  return get_local_sidereal_time(uxtime_d,geo_long,jd_out);
}  

double get_local_sidereal_time(double uxtime_d)
{
   double jd;
   return get_local_sidereal_time(uxtime_d,geo_long,jd);
}

double get_local_sidereal_time_from_jd(double jd,double geo_long_deg)
{
  double sid_greenwich_h = ln_get_apparent_sidereal_time(jd);
//  double sid_greenwich_h = ln_get_mean_sidereal_time (jd_out);
//  printf("DEBUG : sid_greenwich = %.8f [h] = %.8f [deg]\n",sid_greenwich_h,sid_greenwich_h*15.00);
  
  double sid_local_h = sid_greenwich_h + geo_long_deg/15.00;
//  printf("DEBUG : sid_local = %.8f + %.8f = %.8f [h] = %.8f [deg]\n",sid_greenwich_h,geo_long_deg/15.00,sid_local_h,sid_local_h*15.00);
  
  cut_to_range(sid_local_h);
  
  return sid_local_h;     
}

double get_jd(double uxtime_d)
{
   time_t uxtime = (time_t)uxtime_d;
   double jd = ln_get_julian_from_timet( &uxtime );
// printf("DEBUG : jd(%d = floor(%.8f)) = %.8f\n",(int)uxtime,uxtime_d,jd_out);
   jd += (uxtime_d - uxtime)/(24.00*3600.00);
// printf("DEBUG : jd_prim = %.8f\n",jd_out);
  
   return jd;         
}

double get_jd_plus_ndays_sidalligned(double start_ux_d,double n_days)
{
   double jd_step=(0.10/(24.00*3600.00)); // 1 sec 
   double jd_start = get_jd(start_ux_d);   
   double jd_end   = get_jd(start_ux_d+n_days*(24*3600));


   double sid_start_h = get_local_sidereal_time_from_jd(jd_start,geo_long);
   double sid_end_h   = get_local_sidereal_time_from_jd(jd_end,geo_long);
   
   double jd=jd_end;
   if( sid_end_h > sid_start_h ){
      while( sid_end_h > sid_start_h ){
         sid_end_h   = get_local_sidereal_time_from_jd(jd,geo_long);
         jd = jd - jd_step;
      }      
   }else{
      if( sid_end_h < sid_start_h ){
         while( sid_end_h < sid_start_h ){
            sid_end_h   = get_local_sidereal_time_from_jd(jd,geo_long);
            
            jd = jd + jd_step;
         }         
      }
   }
   
   if( gBGPrintfLevel >= 3 ){
      printf("DEBUG : %.8f JD + %.2f days = %.8f JD -> SID alligned JD = %.8f (sid_end=%.8f [h] VS sid_start=%.8f [h])\n",jd_start,n_days,jd_end,jd,sid_end_h,sid_start_h);
   }      
   
   return jd;
}


double get_ang_distance_deg( double ra1, double dec1, double ra2, double dec2 )
{
   struct ln_equ_posn obj1,obj2;
   
   obj1.ra = ra1;
   obj1.dec = dec1;
   
   obj2.ra = ra2;
   obj2.dec = dec2;
      
   double ret = ln_get_angular_separation( &obj1, &obj2 );
   
   return ret;   
}
   
void radec2azh( double ra, double dec, time_t unix_time, double geo_long_deg, double geo_lat_deg, double& out_az, double& out_alt )
{
   ln_equ_posn radec_coord;
   radec_coord.ra = ra;
   radec_coord.dec = dec;
   
   // unix_time -> JD :
   double JD;
   struct ln_date date;         
   get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
   JD = ln_get_julian_day (&date);
                                    
   // position of observatory :
   struct ln_lnlat_posn observer;
   observer.lng = geo_long_deg;
   observer.lat = geo_lat_deg;
                
   struct ln_hrz_posn hrz_pos;
   ln_get_hrz_from_equ( &radec_coord, &observer, JD, &hrz_pos );
   
   out_az  = hrz_pos.az;
   out_alt = hrz_pos.alt;                            
}   

void azh2radec( double az, double alt, time_t unix_time, double geo_long_deg, double geo_lat_deg, double& out_ra, double& out_dec )
{
   struct ln_hrz_posn hrz_pos;
   hrz_pos.az = az;
   hrz_pos.alt = alt;
   
   // unix_time -> JD :
   double JD;
   struct ln_date date;         
   get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
   JD = ln_get_julian_day (&date);
                                    
   // position of observatory :
   struct ln_lnlat_posn observer;
   observer.lng = geo_long_deg;
   observer.lat = geo_lat_deg;
                
   ln_equ_posn radec_coord;
   ln_get_equ_from_hrz( &hrz_pos, &observer, JD, &radec_coord );
   
   out_ra  = radec_coord.ra;
   out_dec = radec_coord.dec;                            
}   

   
double hour_angle( double ra, time_t unix_time )
{
   double uxtime_d = unix_time;
   double jd;
   
   double lst = get_local_sidereal_time( uxtime_d, jd );
   printf("lst = %.8f [h]\n",lst);
   
   double ha = lst - ra/15.00; 
   printf("ha = %.8f [h]\n",ha);
   return ha;
}   

   