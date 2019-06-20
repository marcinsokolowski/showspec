#include "bg_geo.h"
#include <string.h>

//double geo_long=127.885437;
//double geo_lat=-31.826232; 
double geo_long=116.670815;
double geo_lat=-26.703319;

string gSiteName="MWA";


void set_geo_location( const char* szSiteName )
{
    if( strcasecmp(szSiteName,"mwa") == 0 || strcasecmp(szSiteName,"mro") == 0 ){
       geo_long = 116.670815;
       geo_lat  = -26.703319;
       gSiteName = "MWA";
    }
    if( strcmp(szSiteName,"muresk") == 0 || strcmp(szSiteName,"MURESK") == 0 ){
       geo_long=116.68609722;
       geo_lat=-31.74625000;
       gSiteName = "MURESK";
    }
    if( strcmp(szSiteName,"ebo") == 0 || strcmp(szSiteName,"EBO") == 0 || strcasecmp(szSiteName,"ebo201309") == 0 ){
       // 2013-09
       geo_long = 126.30268889; // was 126.3013;  
       geo_lat  = -32.25256389; // was -32.246391;
       gSiteName = "EBO201309";
    }
    if( strcasecmp(szSiteName,"ebo201312") == 0  ){
       // 2013-12 :
       geo_long = 126.29777778; // was 126.3013;  
       geo_lat  = -32.25250000; // was -32.246391;
       gSiteName = "EBO201312";
    }
    if( strcasecmp(szSiteName,"wond20140406") == 0  || strcasecmp(szSiteName,"wond")==0 ){
       // 2014-04-06 :
       geo_long = 118.43999167; // - 27deg 51'10.31''
       geo_lat  = -27.85286389; // 118deg 26' 23.97''
       gSiteName = "WONDINONG_20140406";
    }
    if( strcasecmp(szSiteName,"wond20140405") == 0  ){
       // 2014-04-05 :
       // TODO: change according to Randall's info ! 
       geo_long = 118.43999167; // - 27deg 51'10.31''
       geo_lat  = -27.85286389; // 118deg 26' 23.97''
       gSiteName = "WONDINONG_20140406";
    }    
    if( strcasecmp(szSiteName,"carmel") == 0  ){
       // 2014-04-05 :
       // TODO: change according to Randall's info ! 
       geo_long = 116.096;
       geo_lat  = -32.021;
       gSiteName = "CARMEL";
    }    
}
    

