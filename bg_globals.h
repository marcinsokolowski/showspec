#ifndef _MY_GLOBALS_H__
#define _MY_GLOBALS_H__

#include <time.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <string.h>

#include "bg_date.h"
#include "bg_geo.h"

class cIntRange;

#define BG_MAX_FREQ 480
#define BG_CHANNELS 4096

// calibration point :
struct cValue {
   double x;

   double y; // can be used as Re
   double z; // and Img
   double v;

   cValue() : x(0),y(0),z(-1000) {}
   cValue( double _x, double _y, double _z ) : x(_x), y(_y), z(_z) {}
   
   double Re(){ return y; }
   double Im(){ return z; }
};
            

void init_cvalue_vector(vector<cValue>& vec, int size);

// GLOBAL FLAGS :
extern int gBGPrintfLevel;

using namespace std;

int is_number(const char* string);
int is_digit( char znak );

// file/directories :
const char* getbasename_new(const char* name,string& out);
const char* change_ext(const char* name,const char* new_ext,string& out);
const char* add_postfix(const char* name,const char* new_ext,string& out);


// reading list files :
int bg_read_list( const char* file, vector<string>& out_list );
void print_cmdline(int argc, char * argv[]);

// angles :
double deg2rad(double in_deg );
double rad2deg(double in_rad);

// sort/median :
void my_sort_float( double* ftab, long cnt );
double calc_rms( double* ftab, long cnt );
double calc_avg_rms( double* ftab, long cnt, double& rms, int bOnlyPositive=0 );

// spectrometer basics :
double ch2freq(int channel);
int freq2channel( double freq );

// files, extensions, dirs etc :
int mkdir(const char* path);
const char* getbasename_new(const string& name,string& out);

double mysqr(double x);

// file, string parsing :
int ParseCommaList( char* szList, vector<string>& out_list, const char * sep="," );
int read_file(const char* file,vector<cValue>& out_list, int db2num=0, int x_col=0, int y_col=1, int min_col=1 );
int read_file(const char* file,vector<cIntRange>& out_list, int x_col=0, int y_col=1, int bDoClear=1 );

// mathematical :
double interpolate( vector<cValue>& file_values, double x, int bRe=1 );


inline double interpolate( double x, double x1, double y1, double x2, double y2 )
{ 
   double delta = 0.00;
   if( x2 != x1 ){
      delta = ((y2-y1)/(x2-x1))*(x-x1);
   }
   double interpol_val = y1 + delta;
   return interpol_val;
}   
         
int find_interpol_values( vector<cValue>& list, double x, double& x1, double& y1, double& x2, double& y2, int bRe=1 );
int find_interpol_values( vector<cValue>& list, double x, cValue& prev, cValue& after );
double power(double x,double y);
double db2num(double in_db);

// vector<cValue> functions :
void alloc(vector<cValue>& arr, int _size );

// vector<int>
int find_value( vector<int>& arr, int value );

#endif
