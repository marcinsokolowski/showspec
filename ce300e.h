#ifndef _CE300E_H__
#define _CE300E_H__

#include "bg_globals.h"

#include <vector>
#include <string>
using namespace std;

/*struct sAntGain { 
   double freq;   
   double gain_num;
};*/

#define CE300E_GAIN_POINTS 57      
// #define CE300E_GAIN_POINTS 21
extern string gAntCorrFile;
extern cValue gCE300EGainTab_DataSheet[CE300E_GAIN_POINTS];
extern vector<cValue> gCE300EGainTab;

// int read_file(const char* file,vector<cValue>& out_list, int x_col=0, int y_col=1 );
double get_ce300e_gain(double freq);
double ce300e_gain_fit( double freq );

// #define CE300E_EFF_POINTS 22
#define CE300E_EFF_POINTS 38
extern cValue gCE300E_Efficiency[CE300E_EFF_POINTS];

double get_ce300e_eff(double freq);
            
#endif
