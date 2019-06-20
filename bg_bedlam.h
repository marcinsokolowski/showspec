#ifndef _BEDLAM_SPECTROMETER_H__
#define _BEDLAM_SPECTROMETER_H__

#define BIGHORNS_MAX_FREQ_MHZ 360.00 

class CBedlamSpectrometer 
{
public : 
   CBedlamSpectrometer();

   // spectrometer parameters :
   static int    m_SizeX;
   static double stop_freq;
   static double start_freq;
   static int    m_nSpectralModelPolyN;
   static double gBedlamQualityLimitDBM_BAD;
   static double gBedlamQualityLimitDBM_START;
   
   static double ch2freq( int ch );
   static double freq2ch( double freq_mhz );
   
   // BEDLAM UNITS -> mW , dBm 
   static double power2mW( double freq_mhz, double bedlam_power );
   static double power2dbm( double freq_mhz, double bedlam_power );

   // mW,dBm -> BEDLAM UNITS
   static double dbm2bedlampower( double freq_mhz, double power_dbm );
   static double mW2bedlampower( double freq_mhz, double power_mW );
   
   // conversion curve (in Freq)
   static double spectrum_response_model( double freq );
};

#endif
