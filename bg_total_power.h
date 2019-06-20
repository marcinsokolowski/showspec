#ifndef _BG_TOTAL_POWER_H__
#define _BG_TOTAL_POWER_H__

#include "bg_defines.h"
#include "cvalue_vector.h"
#include <vector>

using namespace std;

struct cTotPowerCut
{
   double min_uxtime;
   double max_uxtime;
   double threshold;
};
         

class CTotalPowerList : public vector<cTotalPowerInfo>
{
public :

   CTotalPowerList();

   double GetLocalMedian( double curr_uxtime, double time_interval, double& sigma_iqr );   
   double GetLocalMedianBothSides( double curr_uxtime, double time_interval, double& sigma_iqr, double maxTotalPower=1e20  );

   // cuts :
   // MAX POWER IN a Channel in dBm :
   static double gMaxChannelPower;   
   static double gTotalPowerMaxValue;
   static vector<cValue> gTotalPowerCuts;
   static vector<cTotPowerCut> gTotalPowerThreshRanges;
   static double GetTotalPowerThreshold( double uxtime );
   
   // median envelope :
   static CValueVector gLocalMedianSigma;
   static double gLocalTotalPowerCutThreshold;
};

#endif
