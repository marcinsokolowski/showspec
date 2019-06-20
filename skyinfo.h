#ifndef _SKYINFO_H__
#define _SKYINFO_H__

// HEALPIX includes :
#include <healpix_base.h>
#include <healpix_map.h>
#include <arr.h>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>


struct cSkyInfo
{
   double intensity;
   pointing coord; // healpix.pix2ang(npix);
   double xp,yp;   // map.project(coord,xp,yp)
   double glon,glat;
   double azim, alt;
   double ra, dec;
                 
   cObjectInfo* pPointSource;
};


                    
                    
#endif
                     