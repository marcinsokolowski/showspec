#include "cvalue_vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bg_array.h"

CValueVector::CValueVector()
{}

/*CValueVector& CValueVector::operator=( const CValueVector& right )  
{
   clear();
   for(int i=0;i<right.size();i++){
      push_back( right[i] );
   }
}*/

void CValueVector::init(int size)
{
  //   vec.clear();
  reserve(size);
  cValue val;
  val.x = 0.00;
  val.y = 0.00;
            
  for(int i=0;i<size;i++){
     push_back( val );
  }   
}                           

int CValueVector::IsFreqExcluded( double freq )
{
   for(int i=0;i<size();i++){
      cValue& excluded_range = (*this)[i];
      
      if( excluded_range.x <= freq && freq <= excluded_range.y ){
         return 1;
      }
   }
   
   return 0;
}

int CValueVector::read_file(const char* file, const char* szComment, int db2num, int x_col, int y_col, int min_col )
{
   int ret=0;
   if( file  && strlen(file) ){
      ret = read_file(file, db2num, x_col, y_col, min_col );
      printf("Read %d values from file %s (%s)\n",ret,file,szComment);
   }else{
      printf("WARNING : no file name provided to read information about : %s\n",szComment);
   }
   
   return ret;
}

int CValueVector::read_file(const char* file, int db2num, int x_col, int y_col, int min_col )     
{
   return ::read_file(file, (*this), db2num, x_col, y_col, min_col );
}


double CValueVector::interpolate( double x, int bRe )
{
   return ::interpolate( (*this), x, bRe );
}

double CValueVector::interpolate_both( double x, double& out_val1, double& out_val2 )
{
   cValue prev,after;
   ::find_interpol_values( (*this), x, prev, after );
   
   out_val1 = ::interpolate( x, prev.x, prev.y, after.x, after.y );
   out_val2 = ::interpolate( x, prev.x, prev.z, after.x, after.z );
   
}

int CValueVector::get_list_around_radius( CValueVector& out_list, double x, double radius )
{
   out_list.clear();
   for(int i=0;i<size();i++){
      cValue& val = (*this)[i];
      
      if( fabs(val.x-x) <= radius ){
         out_list.push_back(val);
      }
   }
   
   return out_list.size();
}

double CValueVector::get_uxtime( int integr )
{
   for(int i=0;i<size();i++){
      cValue& val = (*this)[i];
      int x = (int)val.x;
      
      if( integr == x ){
          return val.y;
      }
   }
   
   return -1;         
}

static int gi0=0;
double CValueVector::get_uxtime_fast( int integr )
{
   for(int i=gi0;i<size();i++){
      cValue& val = (*this)[i];
      int x = (int)val.x;
      
      if( integr == x ){
          gi0=i;
          return val.y;
      }
   }

   // if not found retry from the beginning and reset a starting point:
   gi0=0;
   for(int i=0;i<size();i++){
      cValue& val = (*this)[i];
      int x = (int)val.x;
      
      if( integr == x ){
          return val.y;
      }
   }
   
   return -1;         
}



cValue* CValueVector::find_int_info( int integr )
{
   for(int i=0;i<size();i++){
      cValue& val = (*this)[i];
      int x = (int)val.x;
      
      if( integr == x ){
          return &val;
      }
   }
   
   return NULL;         

}

double CValueVector::find_max_value( double ux_start, double ux_end )
{
   double max_val = -1e20;
   int found=0;
   for(int i=0;i<size();i++){
      cValue& val = (*this)[i];
      double uxtime = val.x;
      
      if( ux_start<=uxtime && uxtime<=ux_end ){
         if( val.y > max_val ){
            max_val = val.y;
         }
         found = 1;
      }
   }    

   if( found > 0 ){   
      return max_val;
   }
   
   return -100000000;
}

int CValueVector::subtract( CValueVector& right, CValueVector& diff )
{
   if( size() != right.size() ){
      printf("ERROR in code : trying to subtract arrays of different sizes %d != %d\n",(int)size(),(int)right.size());
      exit(-1);
   }
   
   int count = size();
   diff.init(count);
   for(int i=0;i<count;i++){
      diff[i].x = (*this)[i].x;
      diff[i].y = (*this)[i].y - right[i].y;
   }
   
   return count;
}

int CValueVector::count_zeros( double min_x, double max_x, double radius )
{
   int zeros = 0;
   
   for(int i=0;i<size();i++){
      double x_val = (*this)[i].x;
      
      if( min_x<=x_val && x_val<=max_x ){
         if( fabs((*this)[i].y) < radius ){
            zeros++;
         }
      }
   }
   
   return zeros;
}

int CValueVector::SaveToFile( const char* outfile )
{
   FILE* out_f = fopen(outfile,"w");
   for(int i=0;i<size();i++){
      fprintf(out_f,"%.8f %.8f\n",(*this)[i].x,(*this)[i].y);
   }
   fclose(out_f);
   
   return 1;
}

void CValueVector::vector2array( CBgArray& out_array )
{
   out_array.alloc(size());   
   for(int i=0;i<size();i++){
      out_array[i] = (*this)[i].y;
   }
}

cValue* CValueVector::find_value( double x, double precision )
{
   for(int i=0;i<size();i++){
      cValue& val = (*this)[i];
      
      double diff = fabs( val.x - x );
      if ( diff <= precision ){
        return &((*this)[i]);
      }
   }
   
   return NULL;
}

void CValueVector::add_value( double x, double y, double z ){
   cValue tmp;
   tmp.x = x;
   tmp.y = y;
   tmp.z = z;
   
   push_back( tmp );
}