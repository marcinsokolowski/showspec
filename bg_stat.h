#ifndef _BG_STAT_H__
#define _BG_STAT_H__

double get_trim_median( double n_sigma_iqr, double* tab, int& cnt, double& sigma_iqr );
double get_trim_median_up( double n_sigma_iqr, double* tab, int& cnt, double& sigma_iqr );
double GetAvgEstimator( double* intab, int cnt, int n_iter, double& sigma_iqr, int x=0, int trim_up=0 );

#endif
