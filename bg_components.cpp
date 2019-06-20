#include "bg_components.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double get_gm_transmiss(double freq)
{
   double g_db=0;
   int found=0;

   if( freq < gGainModuleTransmiss[0].freq ){
      g_db = gGainModuleTransmiss[0].trans;
      found=1;
   }
   if( freq > gGainModuleTransmiss[GAIN_MODULE_POINTS-1].freq ){
      g_db = gGainModuleTransmiss[GAIN_MODULE_POINTS-1].trans;
      found=1;
   }
   
   for(int i=1;i<GAIN_MODULE_POINTS;i++){
      if( freq>=gGainModuleTransmiss[i-1].freq && freq<=gGainModuleTransmiss[i].freq ){
         double f1 = gGainModuleTransmiss[i-1].freq;
         double g1 = gGainModuleTransmiss[i-1].trans;
         double f2 = gGainModuleTransmiss[i].freq;
         double g2 = gGainModuleTransmiss[i].trans;
         
         double g = ((g2-g1)/(f2-f1))*(freq-f1) + g1;
         g_db = g;
         found=1;
         break;
      }
   }
   
   if( !found ){
      printf("ERROR in code - could not find freq=%.2f in gGainModuleTransmiss !\n",freq);
      exit(-1);   
   }
   
   double g_num = exp(log(10)*(g_db/10.00));
   return g_num;
}

sCompChar gGainModuleTransmiss[GAIN_MODULE_POINTS] = 
{
   { 10 , 7.0933 },
   { 17.3333 , 10.3487 },
   { 21.7333 , 11.6512 },
   { 26.1333 , 19.8951 },
   { 32 , 30.0916 },
   { 40.8 , 39.638 },
   { 46.6667 , 41.8085 },
   { 64.2667 , 42.4631 },
   { 80.4 , 42.4666 },
   { 109.733 , 42.473 },
   { 147.867 , 41.6137 },
   { 208 , 40.9761 },
   { 252 , 40.335 },
   { 257.867 , 39.4686 },
   { 278.4 , 37.087 },
   { 306.267 , 33.4055 },
   { 340 , 29.9421 },
   { 372.267 , 26.6954 },
   { 406 , 23.2321 },
   { 430.933 , 20.2007 },
   { 451.467 , 16.5175 },
   { 466.133 , 13.9177 },
   { 476.4 , 11.9677 },
   { 485.2 , 10.0174 },
   { 494 , 8.7178 },
   { 501.333 , 8.2855 },
   { 505.733 , 8.2865 },
   { 516 , 8.5057 },
   { 529.2 , 8.5086 },
   { 539.467 , 7.86 },
   { 548.267 , 6.5605 },
   { 557.067 , 4.827 },
   { 564.4 , 2.8764 },
   { 570.267 , 0.9254 },
   { 576.133 , -1.4595 },
   { 582 , -3.6274 },
   { 586.4 , -4.711 },
   { 592.267 , -4.4928 },
   { 604 , -3.8395 },
   { 612.8 , -3.6206 },
   { 624.533 , -3.835 },
   { 636.267 , -4.4832 },
   { 646.533 , -5.1317 },
   { 658.267 , -5.346 },
   { 675.867 , -5.1253 },
   { 693.467 , -5.5552 },
   { 711.067 , -7.7206 },
   { 716.933 , -9.8885 },
   { 724.267 , -12.273 },
   { 727.2 , -14.0077 },
   { 731.6 , -11.4037 },
   { 737.467 , -8.5825 },
   { 747.733 , -6.6279 },
   { 758 , -5.5411 },
   { 774.133 , -4.6699 },
   { 787.333 , -3.5824 },
   { 807.867 , -2.9271 },
   { 841.6 , -3.1367 },
   { 851.867 , -3.1344 },
   { 863.6 , -3.9995 },
   { 881.2 , -5.0803 },
   { 892.933 , -5.2946 },
   { 901.733 , -5.2927 },
   { 910.533 , -5.9415 },
   { 920.8 , -5.7223 },
   { 931.067 , -6.5878 },
   { 939.867 , -6.5858 },
   { 947.2 , -7.0181 },
   { 951.6 , -7.451 },
   { 963.333 , -7.2315 },
   { 979.467 , -7.2279 },
   { 1000 , -6.5727 }
};
