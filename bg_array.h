#ifndef _BG_ARRAY_H__
#define _BG_ARRAY_H__

#include <stdlib.h>
#include <vector>
using namespace std;

class CValueVector;

class CBgFits;

class CBgArray : public vector<double>
{
public :
  CBgArray( int _size=0, double _values=0 );  
  ~CBgArray();

  void alloc(int _size, double _values=0 );
  double amax();  
  void set_value(double value=0);
  
  CBgArray& Subtract( CBgArray& right );
  CBgArray& Multiply( double val );
  CBgArray& Divide( double val );
  void SaveToFile( const char* outfile , CBgFits* pFits=NULL );
  int find_min_max( double& min_val, double& max_val, double& avg_val );
  void array2vector( CValueVector& out_vector );
};

#endif
