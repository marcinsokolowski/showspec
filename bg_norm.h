#ifndef _BG_NORM_H__
#define _BG_NORM_H__


enum eTimeStampType_T  { eNormalNonUniqueTimeStampFile=1, eUniqueTimeStampFile=2, eUniqueFlagFile=3 };

struct cIntInfo
{
   int int_idx;
   double uxtime;
   double ux_start;
   double ux_end;
   
   double lst; // local sidreal time
   double local_solar_time; // SID - RA_sun
   double sun_alt_deg;      // Sun elevation
   double galaxy_alt_deg;   // Galaxy elevation 
   double t_amb;            // Ambient temperature
   double delta_time;
   
   // solar activity :
   double solar_activity_pre_ao;
   double solar_activity_post_ao;
   int is_sun_ok;
   int is_ok;
};

cIntInfo& find_closest_same_lst_integration( const char* timestamp_fname, cIntInfo& out_info );

int FindClosestMedianInt( cIntInfo& refIntInfo, vector<cIntInfo>& int_info_tab , int size=-1, int min_i=-1, int max_i=-1, double limit=(50.00/86400.00) );

int read_file(const char* file,vector<cIntInfo>& out_list, int bFixUT=0, eTimeStampType_T in_file_type=eNormalNonUniqueTimeStampFile, double minUxTime=-1, double maxUxTime=-1  );

double calc_diff_time( double ra, double ra0 );

int find_norm_int( vector<cIntInfo>& lst_list, double lst, double& min_dist );

#endif
