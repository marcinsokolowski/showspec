#ifndef _SKYMAP_CACHE_H__
#define _SKYMAP_CACHE_H__

#include <vector>
#include <string>
using namespace std;


// HEALPIX includes :
/*#include <healpix_base.h>
#include <healpix_map.h>
#include <arr.h>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>*/

struct cObjectInfo
{
   string name;
   double ra,dec;    // deg   
   double glon,glat; // deg
   double azim,alt;  // deg
   double radius_arcsec; // size of the source 

   double flux;     // Jy
   double err_flux; // Jy
   
   // if periodic (pulsars or others)
   double period; // seconds 
   double pulse_length; // length of pulse in seconds 
   double flux_index; // frequency scaling of flux 
   double flux_ref_freq_mhz; // reference frequency if flux scales with index flux_index 
};

/*struct cSkyInfo
{
   double intensity;
   pointing coord; // healpix.pix2ang(npix);
   double xp,yp;   // map.project(coord,xp,yp)
   double glon,glat;
   double azim, alt;
   double ra, dec;

   cObjectInfo* pPointSource;
};*/

class cSkyInfo;

struct cSkyMap
{
   string szNiFile;
   vector<double> sky_intensities;
};


class CSkyMapCache : public vector<cSkyMap> 
{
// protected :
public :
   static int read_infile( const char* infile, vector<double>& sky_intensities );   
   static int read_infile_binary( const char* infile, vector<double>& sky_intensities );
   int double2skyinfo( vector<double>& sky_temp, vector<cSkyInfo>& sky_intensities );

   int m_bCacheON;
   static int m_bBinaryFile;

public :
   CSkyMapCache( int bCacheON=0 );
   void SetCache( int bCacheON ){ m_bCacheON = bCacheON; }
   int GetCacheONOFF(){ return m_bCacheON; }
   vector<cSkyInfo>* get_skyintensities( const char* infile );
   
   // old version of reading - cache disabled etc, etc 
   static int read_infile( const char* infile );
   static int read_infile_binary( const char* infile );
   

   static vector<cSkyInfo> m_SkyIntensities;   
};

#endif
