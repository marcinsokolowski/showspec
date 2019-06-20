#ifndef _BG_COMPONENTS_H__
#define _BG_COMPONENTS_H__

struct sCompChar
{
  double freq;
  double trans;
};

#define GAIN_MODULE_POINTS 72
extern sCompChar gGainModuleTransmiss[GAIN_MODULE_POINTS];

double get_gm_transmiss(double freq);


#endif
