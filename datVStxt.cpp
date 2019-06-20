#include "skymap_cache.h"
#include "MollweideSkyMap.h"
#include <stdio.h>

class Healpix_Base;
Healpix_Base* pHealpix=NULL;
MollweideSkyMap gMap(100);

int main(int argc,char* argv[])
{
   string datfile = argv[1];
   string txtfile = argv[2];

   vector<double> sky_intensities_dat,sky_intensities_txt;

   CSkyMapCache datcache,txtcache;

   printf("-------------------------------------- BINARY --------------------------------------\n");   
   printf("Reading binary file %s\n",datfile.c_str());
   int datcache_cnt = datcache.read_infile_binary( datfile.c_str(),sky_intensities_dat);
   printf("Read %d elements\n",datcache_cnt);
   printf("\n\n\n");
   
   printf("-------------------------------------- TEXT --------------------------------------\n");
   printf("Reading text file %s\n",txtfile.c_str());
   int txtcache_cnt = txtcache.read_infile( txtfile.c_str(),sky_intensities_txt);
   printf("Read %d elements\n",txtcache_cnt);

   printf("-------------------------------------- COMPARISON --------------------------------------\n");
   int bCompRes=1;
   if( datcache_cnt == txtcache_cnt ){
      printf("Number of elements is equal both = %d\n",datcache_cnt);
   }else{
      printf("ERROR : different number of elements binary file = %d vs txt file = %d\n",datcache_cnt,txtcache_cnt);      
      bCompRes=0;
   }
   
   double max_diff=-1.00;
   int max_diff_i=-1;
   for(int i=0;i<datcache_cnt;i++){
      if( fabs(sky_intensities_txt[i] - sky_intensities_dat[i]) > max_diff ){
         max_diff = fabs(sky_intensities_txt[i] - sky_intensities_dat[i]);
         max_diff_i=i;
      }
   }
   
   printf("Max difference = %e [K] at index = %d (T_dat=%.15f vs T_txt=%.15f)\n",max_diff,max_diff_i,sky_intensities_dat[max_diff_i],sky_intensities_dat[max_diff_i]);
   if( max_diff > 0.000001 ){
      bCompRes=0;
   }
   printf("Comparison final result = %d\n",bCompRes);
}
