#include "skymap_cache.h"
#include "skyinfo.h"
#include <stdio.h>
#include <stdlib.h>
#include <myfile.h>
#include <myparser.h>
#include <mystring.h>
#include <mystrtable.h>
#include <basestructs.h>
#include <bg_globals.h>


// HEALPIX includes :
#include <healpix_base.h>
#include <healpix_map.h>
#include <arr.h>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>
#include "MollweideSkyMap.h"


// binary file reading :
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>


extern Healpix_Base* pHealpix;
extern MollweideSkyMap gMap;  


vector<cSkyInfo> CSkyMapCache::m_SkyIntensities;
int CSkyMapCache::m_bBinaryFile=0;

CSkyMapCache::CSkyMapCache( int bCacheON )
{
   m_bCacheON = bCacheON;
}

vector<cSkyInfo>* CSkyMapCache::get_skyintensities( const char* infile )
{
   double size_b=0.00;
   for(int i=0;i<size();i++){     
      cSkyMap& map = (*this)[i];
      
      size_b += map.sky_intensities.size()*sizeof(double);      
   }
   printf("MEMORY REPORT : amount of memory used for sky maps = %.2f MB\n",size_b/1e6);


   if( m_bCacheON > 0 ){         
      for(int i=0;i<size();i++){
         cSkyMap& map = (*this)[i];
      
         if( strcmp(map.szNiFile.c_str(),infile) == 0 ){
            printf("INFO : data from infile %s found in memory cache\n",infile);
            double2skyinfo( map.sky_intensities , m_SkyIntensities );
            return &m_SkyIntensities;
         }
      }

      printf("INFO : infile %s not yet read, reading file now ...\n",infile);   
      cSkyMap new_skymap_freq;
      new_skymap_freq.szNiFile = infile;
      int n_read_count=0;
      
      if( CSkyMapCache::m_bBinaryFile > 0 ){
         n_read_count = read_infile_binary( infile, new_skymap_freq.sky_intensities );
      }else{
         n_read_count = read_infile( infile, new_skymap_freq.sky_intensities );
      }
      
      if( n_read_count > 0 ){
         push_back(new_skymap_freq);
         double2skyinfo( new_skymap_freq.sky_intensities , m_SkyIntensities );
         return &m_SkyIntensities;
      }            
   }else{
/*      printf("WARNING : caching of skymodel is disabled\n");
      if( size() <= 0 ){
         cSkyMap new_skymap_freq;
         new_skymap_freq.szNiFile = "any";
         push_back(new_skymap_freq);
      }
      read_infile( infile, (*this)[0].sky_intensities );
      double2skyinfo( (*this)[0].sky_intensities , m_SkyIntensities );*/
      
      if( CSkyMapCache::m_bBinaryFile > 0 ){
         read_infile_binary( infile );
      }else{
         read_infile( infile );
      }
      return &m_SkyIntensities;      
   }            
   
   return NULL;
}


int CSkyMapCache::double2skyinfo( vector<double>& sky_temp, vector<cSkyInfo>& sky_intensities )
{
   cSkyInfo skyinfo;

   int bFillSkyMap=0;
   if( sky_intensities.size() == 0 ){
      printf("INFO : filling of sky map is required\n");fflush(stdout);
      bFillSkyMap = 1;
   }

   time_t start = get_dttm();
   for(int npix=0;npix<sky_temp.size();npix++){
      skyinfo.intensity = sky_temp[npix];      
      
      if( bFillSkyMap ){
         skyinfo.coord = pHealpix->pix2ang(npix);
         int i = gMap.project(skyinfo.coord,skyinfo.xp,skyinfo.yp);         

         // TODO : how (phi_deg,theta_deg) are related to (GLON,GLAT) ???
         // HOW they are related to galactic coordinates ???             
         double phi_deg = (skyinfo.coord.phi/M_PI)*180.00;
         double theta_deg = (skyinfo.coord.theta/M_PI)*180.00;
         skyinfo.glon = phi_deg; // ???
         skyinfo.glat = (90.00 - theta_deg);                     
         skyinfo.pPointSource = NULL;
         
         sky_intensities.push_back(skyinfo);
      }

      if( npix >= sky_intensities.size() ){
         printf("ERROR : number of pixels in sky map file %d different than the first map (%d), optimized code does not support such situation !\n",npix,(int)sky_intensities.size());
         exit(-1);
      }
       
      sky_intensities[npix].intensity = skyinfo.intensity;      
      
      if( npix >= 3145720 ){
         printf("%d = %.8f\n",npix,skyinfo.intensity);
      }
   }
   double memory_used_MB = sizeof(cSkyInfo)*sky_intensities.size() / 1000000.00;
   printf("Converted vector<double> -> vector<cSkyInfo> in %d [sec] , memory used = %.2f MB\n",(int)(get_dttm()-start),memory_used_MB);fflush(stdout);
   
   return sky_intensities.size();
}


// Read file containing sky map for given frequency :
int CSkyMapCache::read_infile( const char* infile, vector<double>& sky_intensities )
{
   printf("CSkyMapCache::read_infile(PARAM1,PARAM2) : Reading file %s ...\n",infile);fflush(stdout);
//   ifstream ifile( infile );
   MyIFile ifile( infile );
   char mybuffer [8192];
   string line;


//   ifile.rdbuf()->pubsetbuf(mybuffer,8192);

   int bFillSkyMap=0;
   if( sky_intensities.size() == 0 ){
      printf("INFO : filling of sky map is required\n");fflush(stdout);
      bFillSkyMap = 1;
   }

   time_t start = get_dttm();
// commented out in optimized code :  sky_intensities.clear();
   int npix=0;
//   while ( ifile.good() ){
   const char* pLine=NULL;
   while( pLine = ifile.GetLine( TRUE ) ){
//      getline(ifile,line);
//      sscanf(line.c_str(),"%lf",&(skyinfo.intensity));
      if( mystring::get_first_non_white( pLine )=='#' )
         continue; 
                                      
      MyParser pars=pLine;
      CMyStrTable items;  
      pars.GetItems( items );                                                        
      double intensity = atof( items[0].c_str() );
      
      
      if( bFillSkyMap ){
         sky_intensities.push_back( intensity );
      }

      if( npix >= sky_intensities.size() ){
         printf("ERROR : number of pixels in sky map file %d different than the first map (%d), optimized code does not support such situation !\n",npix,(int)sky_intensities.size());
         exit(-1);
      }
       
      sky_intensities[npix] = intensity;      
      
      npix++;
      
      if( npix >= 3145720 ){
         printf("%d = %.8f\n",npix,intensity);
      }
   }
//   ifile.close();
   double memory_used_MB = sizeof(double)*sky_intensities.size() / 1000000.00;
   printf("Read file %s in %d [sec] into vector<double>, memory used = %.2f MB\n",infile,(int)(get_dttm()-start),memory_used_MB);fflush(stdout);
   
   return sky_intensities.size();
}

// Read file containing sky map for given frequency :
int CSkyMapCache::read_infile_binary( const char* infile, vector<double>& sky_intensities )
{
   time_t start = get_dttm();
   printf("Reading binary (double) file %s ...\n",infile);fflush(stdout);
   int fd = open( infile, O_RDONLY );
   if( fd>0 ){
      int bFillSkyMap=0;
      if( sky_intensities.size() == 0 ){
         printf("INFO : filling of sky map is required\n");fflush(stdout);
         bFillSkyMap = 1;
      }
                        

      int n_bytes_read=0;
      char buffer_char[1000*sizeof(double)];
                    
      int npix=0;                           
      while ( (n_bytes_read = read(fd,buffer_char,1000*sizeof(double)))>0 ){
         double* buffer_double = (double*)buffer_char;
         for(int i=0;i<(n_bytes_read/sizeof(double));i++){
            if( bFillSkyMap ){
               sky_intensities.push_back(buffer_double[i]);
            }

            if( npix >= sky_intensities.size() ){
               printf("ERROR : number of pixels in sky map file %d different than the first map (%d), optimized code does not support such situation !\n",npix,(int)sky_intensities.size());
               exit(-1);
            }
                                                                                     
            sky_intensities[npix] = buffer_double[i];
            npix++;
         }
      }
      close(fd);
   }else{
      printf("ERROR : could not read binary file %s\n",infile);
      exit(-1);
   }
                                                                                    


   double memory_used_MB = sizeof(double)*sky_intensities.size() / 1000000.00;
   printf("Read file %s in %d [sec] into vector<double>, memory used = %.2f MB\n",infile,(int)(get_dttm()-start),memory_used_MB);fflush(stdout);
   
   return sky_intensities.size();
}


// vector<cSkyInfo> sky_intensities; // SKY PIXELS WITH TEMPERATURES AND COORDINATES (optimized code)
// Read file containing sky map for given frequency :
int CSkyMapCache::read_infile( const char* infile )
{
   printf("CSkyMapCache::read_infile(PARAM) Reading file %s ...\n",infile);fflush(stdout);
//   ifstream ifile( infile );
   MyIFile ifile( infile );
   char mybuffer [8192];
   string line;
   cSkyInfo skyinfo;


//   ifile.rdbuf()->pubsetbuf(mybuffer,8192);

   int bFillSkyMap=0;
   if( CSkyMapCache::m_SkyIntensities.size() == 0 ){
      printf("INFO : filling of sky map is required\n");fflush(stdout);
      bFillSkyMap = 1;
   }

   time_t start = get_dttm();
// commented out in optimized code :  CSkyMapCache::m_SkyIntensities.clear();
   int npix=0;
//   while ( ifile.good() ){
   const char* pLine=NULL;
   while( pLine = ifile.GetLine( TRUE ) ){
//      getline(ifile,line);
//      sscanf(line.c_str(),"%lf",&(skyinfo.intensity));
      if( mystring::get_first_non_white( pLine )=='#' )
         continue; 
                                      
      MyParser pars=pLine;
      CMyStrTable items;  
      pars.GetItems( items );                                                        
      skyinfo.intensity = atof( items[0].c_str() );
      
      
      if( bFillSkyMap ){
         skyinfo.coord = pHealpix->pix2ang(npix);
         int i = gMap.project(skyinfo.coord,skyinfo.xp,skyinfo.yp);         

         // TODO : how (phi_deg,theta_deg) are related to (GLON,GLAT) ???
         // HOW they are related to galactic coordinates ???             
         double phi_deg = (skyinfo.coord.phi/M_PI)*180.00;
         double theta_deg = (skyinfo.coord.theta/M_PI)*180.00;
         skyinfo.glon = phi_deg; // ???
         skyinfo.glat = (90.00 - theta_deg);                     
         skyinfo.pPointSource = NULL;
         
         CSkyMapCache::m_SkyIntensities.push_back(skyinfo);
      }

      if( npix >= CSkyMapCache::m_SkyIntensities.size() ){
         printf("ERROR : number of pixels in sky map file %d different than the first map (%d), optimized code does not support such situation !\n",npix,(int)CSkyMapCache::m_SkyIntensities.size());
         exit(-1);
      }
       
      CSkyMapCache::m_SkyIntensities[npix].intensity = skyinfo.intensity;      
      
      npix++;
      
      if( npix >= 3145720 ){
         printf("%d = %.8f\n",npix,skyinfo.intensity);
      }
   }
//   ifile.close();
   double memory_used_MB = sizeof(cSkyInfo)*CSkyMapCache::m_SkyIntensities.size() / 1000000.00;
   printf("Read file %s in %d [sec] , memory used = %.2f MB\n",infile,(int)(get_dttm()-start),memory_used_MB);fflush(stdout);
   
   return CSkyMapCache::m_SkyIntensities.size();
}

// vector<cSkyInfo> sky_intensities; // SKY PIXELS WITH TEMPERATURES AND COORDINATES (optimized code)
// Read file containing sky map for given frequency :
int CSkyMapCache::read_infile_binary( const char* infile )
{
   printf("Reading binary file %s ...\n",infile);fflush(stdout);
   vector<double> sky_intensities;
   int n_count = CSkyMapCache::read_infile_binary(infile,sky_intensities);

   cSkyInfo skyinfo;

   int bFillSkyMap=0;
   if( CSkyMapCache::m_SkyIntensities.size() == 0 ){
      printf("INFO : filling of sky map is required\n");fflush(stdout);
      bFillSkyMap = 1;
   }

   time_t start = get_dttm();
   for(int i=0;i<sky_intensities.size();i++){
      int npix=i;
      skyinfo.intensity = sky_intensities[i];      
      
      if( bFillSkyMap ){
         skyinfo.coord = pHealpix->pix2ang(npix);
         int i = gMap.project(skyinfo.coord,skyinfo.xp,skyinfo.yp);         

         // TODO : how (phi_deg,theta_deg) are related to (GLON,GLAT) ???
         // HOW they are related to galactic coordinates ???             
         double phi_deg = (skyinfo.coord.phi/M_PI)*180.00;
         double theta_deg = (skyinfo.coord.theta/M_PI)*180.00;
         skyinfo.glon = phi_deg; // ???
         skyinfo.glat = (90.00 - theta_deg);                     
         skyinfo.pPointSource = NULL;
         
         CSkyMapCache::m_SkyIntensities.push_back(skyinfo);
      }

      if( npix >= CSkyMapCache::m_SkyIntensities.size() ){
         printf("ERROR : number of pixels in sky map file %d different than the first map (%d), optimized code does not support such situation !\n",npix,(int)CSkyMapCache::m_SkyIntensities.size());
         exit(-1);
      }
       
      CSkyMapCache::m_SkyIntensities[npix].intensity = skyinfo.intensity;      
      
      npix++;
      
      if( npix >= 3145720 ){
         printf("%d = %.8f\n",npix,skyinfo.intensity);
      }
   }
//   ifile.close();
   double memory_used_MB = sizeof(cSkyInfo)*CSkyMapCache::m_SkyIntensities.size() / 1000000.00;
   printf("Read file %s in %d [sec] , memory used = %.2f MB\n",infile,(int)(get_dttm()-start),memory_used_MB);fflush(stdout);
   
   return CSkyMapCache::m_SkyIntensities.size();
}
