#ifndef _WEATHER_STATION_H__
#define _WEATHER_STATION_H__

#include <string>
#include <vector>
using namespace std;

struct cWeatherInfo
{
   double uxtime;
   double temp_c; // temp in Celsius
   double humidity;
};

class CWeatherStation
{
public :
   string m_szWeatherFile;
   vector<cWeatherInfo> m_WeatherList;   
   int m_bIsSorted;

   CWeatherStation( const char* szWeatherFile=NULL );
   
   int ReadTempFile( const char* szWeatherFile=NULL );
   double GetAmbTempK( double uxtime );
   int CheckSorted();
   int GetCount(){ return m_WeatherList.size(); }
   int GetTempList( double ux_start, double ux_end, vector<cWeatherInfo>& temp_list );
};

#endif
