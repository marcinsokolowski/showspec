#ifndef _CVALUE_VECTOR_H__
#define _CVALUE_VECTOR_H__

#include "bg_globals.h"

class CBgArray;


class CValueVector : public vector<cValue>
{
public :
   CValueVector();
//   CValueVector& operator=( const CValueVector& right );

   void init(int size);
   int IsFreqExcluded( double freq );
   int read_file(const char* file, int db2num=0, int x_col=0, int y_col=1, int min_col=1 );
   int read_file(const char* file, const char* szComment, int db2num=0, int x_col=0, int y_col=1, int min_col=1 );
   double interpolate( double x, int bRe=1 );
   double interpolate_both( double x, double& out_val1, double& out_val2 );
   void vector2array( CBgArray& out_array );
   cValue* find_value( double x, double precision=0.01 );
   void add_value( double x, double y, double z=0.00 );
      
   int get_list_around_radius( CValueVector& out_list, double x, double radius );
   
   // when INT UXTIME data :
   double get_uxtime( int integr );
   double get_uxtime_fast( int integr );
   cValue* find_int_info( int integr );
   
   // if solar info files :
   double find_max_value( double ux_start=0, double ux_end=1e20 );
   
   // arithmetic operations:
   int subtract( CValueVector& right, CValueVector& diff );
   int count_zeros( double min_x=-1e20, double max_x=1e20, double radius=0.00001);
   
   int SaveToFile( const char* outfile );
};

#endif
