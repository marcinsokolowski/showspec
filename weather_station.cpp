#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>


#include <mystring.h>
#include <myfile.h>
#include <myparser.h>
#include <mystrtable.h>
#include "bg_globals.h"

#include "weather_station.h"

CWeatherStation::CWeatherStation( const char* szWeatherFile )
: m_bIsSorted(0)
{
   if( szWeatherFile && szWeatherFile[0] ){
      m_szWeatherFile = szWeatherFile;
   }
}

int CWeatherStation::CheckSorted()
{
   double prev_uxtime=0;

   for(int i=0;i<m_WeatherList.size();i++){
      cWeatherInfo& winfo = m_WeatherList[i];
      
      if( winfo.uxtime < prev_uxtime ){
         m_bIsSorted = 0;
         return 0;
      }
      
      prev_uxtime = winfo.uxtime;
   }   
   
   m_bIsSorted = 1;
   return 1;
}

int CWeatherStation::ReadTempFile( const char* szWeatherFile )
{
   if( szWeatherFile && szWeatherFile[0] ){
      m_szWeatherFile = szWeatherFile;
   }

   if( strlen(m_szWeatherFile.c_str())<=0 ){
      printf("ERROR : file name was not provided for CWeatherStation::ReadTempFile nor in constructor !\n");
      return -1;
   }

   MyFile infile(m_szWeatherFile.c_str());
   const char* pLine=NULL;
   m_WeatherList.clear();

   int cnt=0;
   while( pLine = infile.GetLine( TRUE ) ){
      if( mystring::get_first_non_white( pLine )=='#' )
         continue;  
                  
      MyParser pars=pLine;
      CMyStrTable items; 
      pars.GetItems( items );  
 
      if( items.size() >= 2 ){
         cWeatherInfo tmp;
         tmp.uxtime = atof(items[0].c_str());
         tmp.temp_c = atof(items[1].c_str());
         if( items.size() >= 3 ){
            tmp.humidity = atof(items[2].c_str());
         }
         m_WeatherList.push_back(tmp);
         cnt++;
         
//         printf("DEBUG : size() = %d vs cnt = %d\n",m_WeatherList.size(),cnt);
      }
   }
   printf("Read %d values (debug = %d) from file %s\n",(int)m_WeatherList.size(),cnt,m_szWeatherFile.c_str());

   if( CheckSorted() > 0 ){
      printf("INFO weather file %s was time sorted\n",m_szWeatherFile.c_str());
   }

   return m_WeatherList.size();
}

double CWeatherStation::GetAmbTempK( double uxtime )
{
   cWeatherInfo prev_winfo=m_WeatherList[0];
   prev_winfo.uxtime -= 10;

   for(int i=0;i<m_WeatherList.size();i++){
      cWeatherInfo& winfo = m_WeatherList[i];
      
      if( fabs(winfo.uxtime - uxtime )<=1 ){
         return (winfo.temp_c + 273.15); // C -> K 
      }else{
         if( prev_winfo.uxtime<=uxtime && uxtime<=winfo.uxtime ){
            double u1 = prev_winfo.uxtime;
            double t1 = prev_winfo.temp_c;
            double u2 = winfo.uxtime;
            double t2 = winfo.temp_c;
            double u  = uxtime;
         
            double interpol_temp = t1 + (u-u1)*(t2-t1)/(u2-u1);
            if( gBGPrintfLevel ){
               printf("Out_temp = %.2f [K]\n",(interpol_temp + 273.15));
            }
            return (interpol_temp + 273.15); // C -> K ;
         }      
      }

      if( m_bIsSorted ){
         if( (winfo.uxtime-uxtime) > 100 ){
            // already passed required uxtime :
            break;
         }
      }

      prev_winfo = winfo;
   }

   printf("ERROR : could not find temperature for uxtime=%.2f !!!\n",uxtime);

   // not found :   
   return -1000;
}


int CWeatherStation::GetTempList( double ux_start, double ux_end, vector<cWeatherInfo>& temp_list )
{
   temp_list.clear();
   for(int i=0;i<m_WeatherList.size();i++){
      cWeatherInfo& winfo = m_WeatherList[i];

      if( ux_start <= winfo.uxtime && winfo.uxtime <= ux_end ){
         temp_list.push_back( winfo );
      }
   }
   
   return temp_list.size();         
}