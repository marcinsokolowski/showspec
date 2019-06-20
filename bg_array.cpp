#include "bg_array.h"
#include <stdlib.h>
#include <stdio.h>
#include "bg_fits.h"
    
CBgArray::CBgArray( int _size, double _values )
{
    if( _size > 0 ){
       alloc(_size,_values);
    }
}

CBgArray::~CBgArray()
{
}
          

double CBgArray::amax()
{
 double max=-1000000000.00;
 for(int i=0;i<size();i++){
    if( (*this)[i] > max ){
       max = (*this)[i];
    }
 }
 
 return max;
}

void CBgArray::alloc(int _size, double _values)
{
   clear();
   reserve(_size);
   for(int i=0;i<_size;i++){
      push_back(_values);
   }                    
}

void CBgArray::set_value(double value)
{
   for(int i=0;i<size();i++){
      (*this)[i] = value;
   }
}

CBgArray& CBgArray::Subtract( CBgArray& right )
{
   if( size() != right.size() ){
      printf("ERROR : sized of arrays to be subtracted differ %d != %d\n",(int)size(),(int)right.size());
      return (*this);
   }
   
   for(int i=0;i<size();i++){
      (*this)[i] = (*this)[i] - right[i];      
   }
   
   return (*this);
}

CBgArray& CBgArray::Multiply( double val )
{
   for(int i=0;i<size();i++){
      (*this)[i] = (*this)[i] * val; 
   }
   
   return (*this);
}

CBgArray& CBgArray::Divide( double val )
{
   for(int i=0;i<size();i++){
      (*this)[i] = (*this)[i] / val; 
   }
   
   return (*this);
}


void CBgArray::SaveToFile( const char* outfile, CBgFits* pFits  )
{
   FILE* outf = fopen(outfile,"w");
   
   for(int i=0;i<size();i++){
      if( pFits ){
         fprintf(outf,"%.3f %.8f\n",pFits->ch2freq(i),(*this)[i]);         
      }else{
         fprintf(outf,"%d %.8f\n",i,(*this)[i]);
      }
   }
   
   fclose(outf);  
}

int CBgArray::find_min_max( double& min_val, double& max_val, double& avg_val )
{
   min_val=1e20;
   max_val=-1e20;
   avg_val=0;
   
   int cnt=0;
   for(int i=0;i<size();i++){
      double val = (*this)[i];
      
      if( val > max_val ){
         max_val = val;
      }
      if( val < min_val ){
         min_val = val;
      }
       
      avg_val += val;
   }
   
   avg_val = avg_val / size();
   
   return size();
}

void CBgArray::array2vector( CValueVector& out_vector )
{
   out_vector.clear();
//   out_vector.assign(size());
   
   for(int i=0;i<size();i++){ 
      cValue tmp;   
   
      tmp.x = 0;
      tmp.z = 0;
      tmp.y = (*this)[i];
      
      out_vector.push_back(tmp);
   }
}
