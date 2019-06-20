#include <libnova/transform.h>
#include <libnova/julian_day.h>
#include <libnova/utility.h>
#include <stdio.h>

/* void get_ymd_hms_ut( time_t ut_time, int& year, int& month, int& day,
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

        // UT date and time 
        get_ymd_hms_ut( unix_time, date.years, date.months, date.days, date.hours, date.minutes, date.seconds );
        JD = ln_get_julian_day (&date);

	  // calculate horizontal coordinates :
	  ln_get_hrz_from_equ(&equ_pos, &observer, JD, &azim_alt_pos );

	  out_azim = azim_alt_pos.az;
	  out_alt = azim_alt_pos.alt;

}*/
