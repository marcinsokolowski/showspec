#include "bg_stat.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bg_globals.h"

double get_trim_median( double n_sigma_iqr, double* tab, int& cnt, double& sigma_iqr )
{   
   double* newtab = new double[cnt];

   int q75= (int)(cnt*0.75);
   int q25= (int)(cnt*0.25);
   double iqr = tab[q75]-tab[q25];
   sigma_iqr = iqr/1.35;
   double range = sigma_iqr*n_sigma_iqr;
   double median = tab[(int)cnt/2];
   
   int newcnt=0;
   for(int i=0;i<cnt;i++){
      if( fabs(tab[i]-median) <= range ){
         newtab[newcnt] = tab[i];
         newcnt++;
      }
   }
                                          
   double ret=newtab[newcnt/2];
   if( gBGPrintfLevel>=2 ){
      printf("\tmedian = %.4f, iqr = %.4f -> sigma_iqr = %.4f -> range = %.4f -> new_median = %.4f [diff = %.4f]\n",median,iqr,sigma_iqr,range,ret,fabs(median-ret));
   }

   // returning smaller array :
   cnt = newcnt;
   for(int i=0;i<newcnt;i++){
      tab[i] = newtab[i];
   }   

   delete [] newtab;
   return ret;
}
                                                
double get_trim_median_up( double n_sigma_iqr, double* tab, int& cnt, double& sigma_iqr )
{   
   double* newtab = new double[cnt];

   int q75= (int)(cnt*0.75);
   int q25= (int)(cnt*0.25);
   double iqr = tab[q75]-tab[q25];
   sigma_iqr = iqr/1.35;
   double range = sigma_iqr*n_sigma_iqr;
   double median = tab[(int)cnt/2];
   
   int newcnt=0;
   for(int i=0;i<cnt;i++){
      if( (tab[i]-median) < range ){
         newtab[newcnt] = tab[i];
         newcnt++;
      }
   }
                                          
   double ret=newtab[newcnt/2];
   if( gBGPrintfLevel>=2 ){
      printf("\tmedian = %.4f, iqr = %.4f -> sigma_iqr = %.4f -> range = %.4f -> new_median = %.4f [diff = %.4f]\n",median,iqr,sigma_iqr,range,ret,fabs(median-ret));
   }

   // returning smaller array :
   cnt = newcnt;
   for(int i=0;i<newcnt;i++){
      tab[i] = newtab[i];
   }   

   delete [] newtab;
   return ret;
}
                                                
                                                
// INPUT  : 
// intab  : sorted table                                                 
// cnt    : number of elements in a table
// n_iter : number of iterations 
double GetAvgEstimator( double* intab, int cnt, int n_iter, double& sigma_iqr, int x, int trim_up )
{
   double* tab = new double[cnt]; // working copy 
   for(int i=0;i<cnt;i++){
      tab[i] = intab[i];
   }

/*   int q75= (int)(cnt*0.75);
   int q25= (int)(cnt*0.25);
   double iqr = tab[q75]-tab[q25];
   double sigma_iqr = iqr/1.35;
   double median = tab[(int)cnt/2];
   printf("median = %.4f, iqr = %.4f -> sigma_iqr = %.4f\n",median,iqr,sigma_iqr);*/

   if( gBGPrintfLevel>=2 ){
      printf("\n\nChannel %d / %.4f [MHz]\n",x,ch2freq(x));
   }                  
   int newcnt=cnt;
   double out_median = tab[cnt/2];
   for(int i=0;i<n_iter;i++){
      if( gBGPrintfLevel>=2 ){
         printf("\t------ Iteration = %d -----\n",i);
      }

      if( trim_up > 0 ){
         out_median = get_trim_median_up( 5.00, tab, newcnt, sigma_iqr );
      }else{
         out_median = get_trim_median( 5.00, tab, newcnt, sigma_iqr );
      }
   }   
   
   delete [] tab;
   
   return out_median;   
}
