#include "bg_total_power.h"
#include "bg_stat.h"
#include "bg_globals.h"
#include <math.h>

// total power cut threshold specified for different time ranges :
vector<cTotPowerCut> CTotalPowerList::gTotalPowerThreshRanges;
vector<cValue> CTotalPowerList::gTotalPowerCuts;

// RFI rejection :  
double CTotalPowerList::gTotalPowerMaxValue=4e10; // 4e10 is ok for Muresk Data 
// MAX POWER IN a Channel in dBm :
double CTotalPowerList::gMaxChannelPower=1000000.00;
   


CValueVector CTotalPowerList::gLocalMedianSigma;
double CTotalPowerList::gLocalTotalPowerCutThreshold=1;
   

CTotalPowerList::CTotalPowerList()
{

}


double CTotalPowerList::GetLocalMedian( double curr_uxtime, double time_interval, double& sigma_iqr  )
{
   int last = size()-1;
   vector<double> tmp_list;
   
   int i = last;
   while( (*this)[i].t_int >= (curr_uxtime - time_interval) ){
//   while( fabs( (*this)[i].t_int - curr_uxtime ) <= time_interval ){   
      tmp_list.push_back( (*this)[i].total_sum );      
      i--;
   }
   
   double* tmp_list2 = new double[tmp_list.size()];
   for(i=0;i<tmp_list.size();i++){
      tmp_list2[i] = tmp_list[i];
   }   
   my_sort_float(  tmp_list2, tmp_list.size() );
   
   int trim_upper=1;
//   double median = GetAvgEstimator( tmp_list2, tmp_list.size(), 3, sigma_iqr, 0, trim_upper );
//   double median = GetAvgEstimator( tmp_list2, tmp_list.size(), 0, sigma_iqr, 0, trim_upper );
   double median = GetAvgEstimator( tmp_list2, tmp_list.size(), 10, sigma_iqr, 0, trim_upper );
   
   delete [] tmp_list2;
   return median;
}

double CTotalPowerList::GetLocalMedianBothSides( double curr_uxtime, double time_interval, double& sigma_iqr, double maxTotalPower  )
{
   vector<double> tmp_list;
   
   for(int i=0;i<size();i++){
      if ( fabs( (*this)[i].t_int - curr_uxtime ) <= time_interval ){   
         if( (*this)[i].total_sum_bedlam < maxTotalPower ){ // check in BEDLAM units, but save units as saved in the list. As thresholds are defined in arbitrary units 
            tmp_list.push_back( (*this)[i].total_sum );
         }
      }else{
         if( tmp_list.size() > 0 ){
            // break the loop in case something already found and we are out of range 
            break;
         }
      }
   }
   
   double* tmp_list2 = new double[tmp_list.size()];
   for(int i=0;i<tmp_list.size();i++){
      tmp_list2[i] = tmp_list[i];
   }   
   my_sort_float(  tmp_list2, tmp_list.size() );
   
   int trim_upper=1;
//   double median = GetAvgEstimator( tmp_list2, tmp_list.size(), 3, sigma_iqr, 0, trim_upper );
//   double median = GetAvgEstimator( tmp_list2, tmp_list.size(), 0, sigma_iqr, 0, trim_upper );
   double median = GetAvgEstimator( tmp_list2, tmp_list.size(), 10, sigma_iqr, 0, trim_upper );
   
   delete [] tmp_list2;
   return median;
}

double CTotalPowerList::GetTotalPowerThreshold( double uxtime )
{
   if( gTotalPowerThreshRanges.size() > 0 ){
      for(int i=0;i<gTotalPowerThreshRanges.size();i++){
         cTotPowerCut& cut = gTotalPowerThreshRanges[i];

         if( cut.min_uxtime<=uxtime && uxtime<=cut.max_uxtime ){
            return cut.threshold;
         }
      }
   }

   return gTotalPowerMaxValue;
}



