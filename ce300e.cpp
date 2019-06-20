#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ce300e.h"

#include <mystring.h>
#include <myfile.h>
#include <myparser.h>
#include <mystrtable.h>

static double sqr(double x)
{
   return (x*x);
}

double get_ce300e_gain(double freq)
{
   if( gCE300EGainTab.size() <= 0 ){
      printf("WARNING : gain array is empty, filling with default values from the datasheet\n");
      for(int i=0;i<CE300E_GAIN_POINTS;i++){
         gCE300EGainTab.push_back(gCE300EGainTab_DataSheet[i]);
      }
   }

   int size = gCE300EGainTab.size();

   if( freq < gCE300EGainTab[0].x ){
      return gCE300EGainTab[0].y;
   }
   if( freq > gCE300EGainTab[size-1].x ){
      return gCE300EGainTab[size-1].y;
   }
   
   for(int i=1;i<gCE300EGainTab.size();i++){
      if( freq>=gCE300EGainTab[i-1].x && freq<=gCE300EGainTab[i].x ){
         double f1 = gCE300EGainTab[i-1].x;
         double g1 = gCE300EGainTab[i-1].y;
         double f2 = gCE300EGainTab[i].x;
         double g2 = gCE300EGainTab[i].y;
         
         double g = ((g2-g1)/(f2-f1))*(freq-f1) + g1;
         return g;
         
      }
   }
   
   printf("ERROR in code - could not find freq=%.2f in gCE300EGainTab !\n",freq);
   exit(-1);   
}

double ce300e_gain_fit( double freq )
{
   double ret = 1.45*sqr(sin(2.0*M_PI*(0.351/300.0)*freq));

   return ret;   
}   

string gAntCorrFile;
vector<cValue> gCE300EGainTab;

cValue gCE300EGainTab_DataSheet[CE300E_GAIN_POINTS] = 
{
   cValue( 20 , 0.01 , 0 ),
   cValue( 25 , 0.01 , 0 ),
   cValue( 30 , 0.02 , 0 ),
   cValue( 35 , 0.03 , 0 ),
   cValue( 40 , 0.05 , 0 ),
   cValue( 45 , 0.08 , 0 ),
   cValue( 50 , 0.12 , 0 ),
   cValue( 55 , 0.17 , 0 ),
   cValue( 60 , 0.22 , 0 ),
   cValue( 65 , 0.29 , 0 ),
   cValue( 70 , 0.37 , 0 ),
   cValue( 75 , 0.45 , 0 ),
   cValue( 80 , 0.55 , 0 ),
   cValue( 85 , 0.66 , 0 ),
   cValue( 90 , 0.73 , 0 ),
   cValue( 95 , 0.83 , 0 ),
   cValue( 100 , 0.87 , 0 ),
   cValue( 105 , 0.90 , 0 ),
   cValue( 110 , 0.97 , 0 ),
   cValue( 115 , 0.98 , 0 ),
   cValue( 120 , 0.98 , 0 ),
   cValue( 125 , 1.01 , 0 ),
   cValue( 130 , 1.02 , 0 ),
   cValue( 135 , 1.01 , 0 ),
   cValue( 140 , 1.03 , 0 ),
   cValue( 145 , 1.06 , 0 ),
   cValue( 150 , 1.08 , 0 ),
   cValue( 155 , 1.16 , 0 ),
   cValue( 160 , 1.20 , 0 ),
   cValue( 165 , 1.31 , 0 ),
   cValue( 170 , 1.39 , 0 ),
   cValue( 175 , 1.37 , 0 ),
   cValue( 180 , 1.36 , 0 ),
   cValue( 185 , 1.34 , 0 ),
   cValue( 190 , 1.29 , 0 ),
   cValue( 195 , 1.32 , 0 ),
   cValue( 200 , 1.30 , 0 ),
   cValue( 205 , 1.30 , 0 ),
   cValue( 210 , 1.28 , 0 ),
   cValue( 215 , 1.31 , 0 ),
   cValue( 220 , 1.40 , 0 ),
   cValue( 225 , 1.47 , 0 ),
   cValue( 230 , 1.43 , 0 ),
   cValue( 235 , 1.39 , 0 ),
   cValue( 240 , 1.39 , 0 ),
   cValue( 245 , 1.35 , 0 ),
   cValue( 250 , 1.34 , 0 ),
   cValue( 255 , 1.24 , 0 ),
   cValue( 260 , 1.15 , 0 ),
   cValue( 265 , 1.14 , 0 ),
   cValue( 270 , 1.16 , 0 ),
   cValue( 275 , 1.23 , 0 ),
   cValue( 280 , 1.25 , 0 ),
   cValue( 285 , 1.21 , 0 ),
   cValue( 290 , 1.14 , 0 ),
   cValue( 295 , 1.13 , 0 ),
   cValue( 300 , 1.16 , 0 )
};

// from file /home/msok/Desktop/CAASTRO/bighorns/doc/hardware/antenna/biconical/simulation/EOR_bicone_efficiency.csv
cValue gCE300E_Efficiency[CE300E_EFF_POINTS] = 
{
/* Aziz 50-250 MHz (2012-07)  
   { 0  , 0.000000 },
   { 50 , 0.717177 },
   { 60 , 0.833593 },
   { 70 , 0.884053 },
   { 80 , 0.912624 },
   { 90 , 0.931636 },
   { 100 , 0.943663 },
   { 110 , 0.951783 },
   { 120 , 0.956075 },
   { 130 , 0.951613 },
   { 140 , 0.959683 },
   { 150 , 0.964269 },
   { 160 , 0.965972 },
   { 170 , 0.966614 },
   { 180 , 0.967048 },
   { 190 , 0.967095 },
   { 200 , 0.966608 },
   { 210 , 0.965603 },
   { 220 , 0.964064 },
   { 230 , 0.962111 },
   { 240 , 0.959854 },
   { 250 , 0.957447 }   */
   
// MS 30-400 MHz (2013-09) 
  cValue( 30 , 0.154206 , 0 ),
  cValue( 40 , 0.550876 , 0 ),
  cValue( 50 , 0.717155 , 0 ),
  cValue( 60 , 0.833604 , 0 ),
  cValue( 70 , 0.884048 , 0 ),
  cValue( 80 , 0.912623 , 0 ),
  cValue( 90 , 0.931638 , 0 ),
  cValue( 100 , 0.943661 , 0 ),
  cValue( 110 , 0.951785 , 0 ),
  cValue( 120 , 0.956075 , 0 ),
  cValue( 130 , 0.951617 , 0 ),
  cValue( 140 , 0.959685 , 0 ),
  cValue( 150 , 0.964269 , 0 ),
  cValue( 160 , 0.965973 , 0 ),
  cValue( 170 , 0.966614 , 0 ),
  cValue( 180 , 0.967047 , 0 ),
  cValue( 190 , 0.967096 , 0 ),
  cValue( 200 , 0.966609 , 0 ),
  cValue( 210 , 0.965605 , 0 ),
  cValue( 220 , 0.964063 , 0 ),
  cValue( 230 , 0.962111 , 0 ),
  cValue( 240 , 0.959855 , 0 ),
  cValue( 250 , 0.957447 , 0 ),
  cValue( 260 , 0.954105 , 0 ),
  cValue( 270 , 0.951482 , 0 ),
  cValue( 280 , 0.948335 , 0 ),
  cValue( 290 , 0.945365 , 0 ),
  cValue( 300 , 0.941507 , 0 ),
  cValue( 310 , 0.939129 , 0 ),
  cValue( 320 , 0.936803 , 0 ),
  cValue( 330 , 0.933675 , 0 ),
  cValue( 340 , 0.929477 , 0 ),
  cValue( 350 , 0.92655 , 0 ),
  cValue( 360 , 0.923564 , 0 ),
  cValue( 370 , 0.92041 , 0 ),
  cValue( 380 , 0.91685 , 0 ),
  cValue( 390 , 0.912201 , 0 ),
  cValue( 400 , 0.907542, 0 )
};

double get_ce300e_eff(double freq)
{
   if( freq < gCE300E_Efficiency[0].x ){
      return gCE300E_Efficiency[0].y;
   }
   if( freq > gCE300E_Efficiency[CE300E_EFF_POINTS-1].x ){
      return gCE300E_Efficiency[CE300E_EFF_POINTS-1].y;
   }
   
   for(int i=1;i<CE300E_EFF_POINTS;i++){
      if( freq>=gCE300E_Efficiency[i-1].x && freq<=gCE300E_Efficiency[i].x ){
         double f1 = gCE300E_Efficiency[i-1].x;
         double g1 = gCE300E_Efficiency[i-1].y;
         double f2 = gCE300E_Efficiency[i].x;
         double g2 = gCE300E_Efficiency[i].y;
         
         double g = ((g2-g1)/(f2-f1))*(freq-f1) + g1;
         return g;
         
      }
   }
   
   printf("ERROR in code - could not find freq=%.2f in gCE300E_Efficiency !\n",freq);
   exit(-1);   
}


/*
cValue gCE300EGainTab[CE300E_GAIN_POINTS] = 
{

{ 50 , 0.69351 },
{ 60 , 0.806918 },
{ 70 , 0.863631 },
{ 80 , 0.90724 },
{ 90 , 0.927257 },
{ 100 , 0.938945 },
{ 110 , 0.949308 },
{ 120 , 0.955406 },
{ 130 , 0.951518 },
{ 140 , 0.959683 },
{ 150 , 0.964269 },
{ 160 , 0.965972 },
{ 170 , 0.966421 },
{ 180 , 0.966274 },
{ 190 , 0.965548 },
{ 200 , 0.963805 },
{ 210 , 0.961644 },
{ 220 , 0.958665 },
{ 230 , 0.957685 },
{ 240 , 0.955823 },
{ 250 , 0.953617 }

};*/

/*int read_file(const char* file,vector<cValue>& out_list, int x_col, int y_col )
{
   MyFile infile(file);
   const char* pLine=NULL;
   out_list.clear();

   while( pLine = infile.GetLine( TRUE ) ){
      if( mystring::get_first_non_white( pLine )=='#' )
         continue;

      MyParser pars=pLine;
      CMyStrTable items;
      pars.GetItems( items );

      int max_col = MAX(x_col,y_col);

      if( items.size() > max_col ){
         cValue tmp;
         tmp.x = atof(items[x_col].c_str());
         tmp.y = atof(items[y_col].c_str());
         out_list.push_back(tmp);
      }
   }
   printf("Read %d values from file %s\n",(int)out_list.size(),file);

   return out_list.size();
}*/


