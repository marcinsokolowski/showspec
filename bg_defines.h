#ifndef _BG_DEFINES_H__
#define _BG_DEFINES_H__

#include <string>

using namespace std;

struct cTotalPowerInfo
{
   string fits_file;
   int fits_int;

   double t;
   double t_int;
   long double total_sum;
   long double total_sum_bedlam;
   long double rms_total_sum;   
   double diff_ratio;
   double rms_local;
   double avg_local;
   double kurtosis;
   double kurtosis_milsat;
   double fit_chi2;

   int ok; // state 

   // temporary value - used for any processing purposes :
   int tmp_value;   
   
   // maximum power :
   double maxPowerSingleChannel;
   double maxPowerFreq;
   
   // orbcomm :
   double orbcommPower;
};

struct cAntIntRange
{
   string fits_file;
   int ant_start;
   int ant_end;
   int state;
   int total_power_idx;
};

struct cRangeDef
{
   double min_total_power_val;
   double max_total_power_val;
   int state; // 0,1,2,3 ...
   string m_szName; // state name TERM,ANT, etc
   int n_integrations; // optional how long is expected state 
   
   cRangeDef() : min_total_power_val(0),max_total_power_val(0),state(0),n_integrations(-1)
   {
   }
};



#endif
