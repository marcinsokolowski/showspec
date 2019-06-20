#include "bg_globals.h"
#include "bg_norm.h"
#include "libnova_interface.h"

#include <myfile.h>
#include <myparser.h>
#include <mystrtable.h>

#include <math.h>

double fix_ut( double ux )
{
   if( 1349049600 <= ux && ux <= 1352505600 ){
      // Muresk run 201210-201211 :
      double drift = 21.5;
      double dt = 711+((ux-1351753823.00)/86400.00)*drift;
               
      double fixed_ux = (ux - dt);
      return fixed_ux;
   }
   
   return ux;
}

int parse_non_unique_file( CMyStrTable& items, cIntInfo& tmp, int bFixUT )
{
      tmp.int_idx = atof(items[2].c_str());
      tmp.uxtime  = atof(items[3].c_str());
      if( bFixUT ){
         tmp.uxtime = fix_ut( tmp.uxtime );
      }        
                 
                 
      if( items.size() >= 5 ){
         tmp.lst     = atof(items[4].c_str());
         if( items.size() >= 6 ){
            tmp.galaxy_alt_deg = atof(items[5].c_str());
            if( items.size() >= 7 ){
               tmp.sun_alt_deg = atof(items[6].c_str());
               if( items.size() >= 8 ){
                  tmp.t_amb = atof(items[7].c_str());
               }
            }   
         }      
                
         // calculation of Sun's coordinates :
         double sun_ra, sun_dec, sun_az, sun_alt;
         get_sun_all( tmp.uxtime, geo_long, geo_lat, sun_ra, sun_dec, sun_az, sun_alt );
         
         tmp.local_solar_time = tmp.lst - (sun_ra/15.00) + 12.00; // Local solar time 
      }else{
         printf("WARNING : no extra info provided in timestamp file\n");
      }
}


int parse_unique_file( CMyStrTable& items, cIntInfo& tmp, int bFixUT )
{
      tmp.int_idx = atof(items[1].c_str());
      tmp.uxtime  = atof(items[3].c_str());
      if( bFixUT ){
         tmp.uxtime = fix_ut( tmp.uxtime );
      }        
                 
                 
      if( items.size() >= 5 ){
         tmp.lst     = atof(items[4].c_str());
         tmp.ux_start = atof(items[5].c_str());
         tmp.ux_end   = atof(items[6].c_str());
         tmp.delta_time = atof(items[7].c_str());
                  
         // calculation of Sun's coordinates :
         double sun_ra, sun_dec, sun_az, sun_alt;
         get_sun_all( tmp.uxtime, geo_long, geo_lat, sun_ra, sun_dec, sun_az, sun_alt );
         
         tmp.local_solar_time = tmp.lst - (sun_ra/15.00) + 12.00; // Local solar time 
      }else{
         printf("WARNING : no extra info provided in timestamp file\n");
      }
}

int parse_unique_flag_file( CMyStrTable& items, cIntInfo& tmp )
{
   tmp.int_idx = atol(items[0].c_str());
   tmp.is_ok   = atol(items[1].c_str());
   tmp.t_amb   = atof(items[4].c_str());
   tmp.sun_alt_deg = atof(items[5].c_str());
   tmp.galaxy_alt_deg = atof(items[6].c_str());
   tmp.uxtime = get_unixtime_from_local_string2( items[3].c_str() , "%Y%m%d_%H%M%S" );
   
   return 5;
}
                     
int read_file(const char* file,vector<cIntInfo>& out_list, int bFixUT, eTimeStampType_T in_file_type, double minUxTime, double maxUxTime  )
{
   MyFile infile(file);
   const char* pLine=NULL;
   out_list.clear();

   printf("Reading file %s of type %d\n",file,in_file_type);
   while( pLine = infile.GetLine( TRUE ) ){
      if( mystring::get_first_non_white( pLine )=='#' )
         continue;
                  
      MyParser pars=pLine;
      CMyStrTable items;  
      pars.GetItems( items );
      
      cIntInfo tmp;

      switch( in_file_type )
      {
         case eNormalNonUniqueTimeStampFile :
//             printf("Reading file %s of type eNormalNonUniqueTimeStampFile\n",file);
             parse_non_unique_file( items, tmp, bFixUT );
             break;
         case eUniqueTimeStampFile :
//             printf("Reading file %s of type eUniqueTimeStampFile\n",file);
             parse_unique_file( items, tmp, bFixUT );
             break;
         case eUniqueFlagFile:
             parse_unique_flag_file( items, tmp );
             break;
         default :
             printf("ERROR : unknown type of file %s\n",file);
      }

      if( minUxTime>0 && tmp.uxtime < minUxTime ){ // skip to small unixtime < minUxTime
         printf("Timestamps before minimal unix time = %.2f -> skipped\n",minUxTime);
         continue;
      }       
      if( maxUxTime>0 && tmp.uxtime >= maxUxTime ){ // skip if past maximum unixtime  maxUxTime
         printf("Timestamps above maximum allowed unixtime = %.2f ignored\n",maxUxTime);
         break;
      }
       
      out_list.push_back(tmp);
   }
   printf("Read %d values from file %s\n",(int)out_list.size(),file);
   
   return out_list.size();
}                         
                          


double calc_diff_time( double ra, double ra0 )
{
  double diff_ra=(ra-ra0);
  if( ra > 12.00 && ra0<6.00 ){
     // ra=23.43 ra0=0.23
     diff_ra = (ra-(ra0+24.00));
  }
  if( ra<6.00 && ra0>12.00 ){
     // ra=0.23 ra0=23.4
     diff_ra = ( ra+24.00 - ra0 );
  }
  return fabs(diff_ra);
}
                              

double GetCriteriaDistance( cIntInfo& refIntInfo, cIntInfo& currIntInfo )
{
   double dist = 10000000.00;
   dist = calc_diff_time(refIntInfo.lst , currIntInfo.lst);
   return dist;
}                             

int FindClosestMedianInt( cIntInfo& refIntInfo, vector<cIntInfo>& int_info_tab , int size, int min_i, int max_i, double limit )
{
   if( size <= 0 ){
      size = int_info_tab.size();
   }         

   if( max_i < 0 ){
      max_i = size;
   }
    
   cIntInfo* pRetInt = NULL;
   int ret_int=-1;
   double minDist = 100000000.00;
   double prev_dist=-1;
   for(int i=0;i<size;i++){
      cIntInfo& currIntInfo = int_info_tab[i];
      if (i>=min_i && i<=max_i && (currIntInfo.uxtime-refIntInfo.uxtime)>(12*3600) ){  // separated by at least 12hours - added 20141114
         double dist = GetCriteriaDistance( refIntInfo , currIntInfo );
         if( gBGPrintfLevel >= 3 ){
            printf("Integration %d : distance %.2f\n",i,dist);
         }
         if( dist < minDist ){
            minDist = dist;   
            pRetInt = &currIntInfo;
            ret_int = i;
            printf("Found minimum distance %.2f at integration = %d\n",dist,i);
         }else{
            if( dist > prev_dist ){
               // distance increases again, so we passed first minimum 
               if( ret_int > 0 ){
                  // new minimum found :
                  if( minDist < limit ){
                     printf("Found 1st closest LST\n");
                     return ret_int;
                  }
               }   
            }      
         }
         prev_dist = dist;
      }
   }   
       
   return ret_int;
}
 
double Dist2Sec( double dist )
{
   double dist_in_sec = dist*3600.00;
   return dist_in_sec;   
}


cIntInfo& find_closest_same_lst_integration( const char* timestamp_fname, cIntInfo& out_info )
{
  vector<cIntInfo> int_info_tab;
  printf("Reading timestamp file %s ...",timestamp_fname);fflush(stdout);
  if( read_file(timestamp_fname,int_info_tab) <= 0 ){
     printf("ERROR : could not read integration timestamp file\n");
     exit(-1);
  }
  printf("OK\n");fflush(stdout);

  // find 
  cIntInfo& intZeroInfo = int_info_tab[0];
  int next_int_same_lst = FindClosestMedianInt( intZeroInfo, int_info_tab, int_info_tab.size(), 100 );
  cIntInfo& sameLstInt = int_info_tab[next_int_same_lst];
  double dist = fabs(intZeroInfo.lst-sameLstInt.lst);
  printf("LST(int 0) = %.8f -> the next closest one is at line = %d, integration = %d, uxtime=%.8f, lst=%.8f (distance = %.8f days = %.2f sec)\n",intZeroInfo.lst,next_int_same_lst,sameLstInt.int_idx,sameLstInt.uxtime,sameLstInt.lst,dist,dist*86400);
  fflush(stdout);
  
  out_info = sameLstInt;
  return out_info;
}


int find_norm_int( vector<cIntInfo>& lst_list, double lst, double& min_dist )
{
   min_dist=10000000.00;
   int best_int=-1;
   
   
   for(int i=0;i<lst_list.size();i++){
      if( fabs( lst_list[i].lst - lst ) < min_dist ){
         min_dist = fabs( lst_list[i].lst - lst );
         best_int = lst_list[i].int_idx;
      }   
   }
   
   return best_int;
}
