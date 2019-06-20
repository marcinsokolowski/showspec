#include "bg_globals.h"
#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string>

using namespace std;

#include <myparser.h>
#include <mystrtable.h>
#include <myfile.h>

#include "bg_fits.h"

// GLOBAL FLAGS :
int gBGPrintfLevel=0;


#define REPLACE_ELEMS( tab, pos1, pos2 ) { tmp=tab[pos1]; tab[pos1]=tab[pos2]; tab[pos2]=tmp; }

/*static int is_number(const char* string)
{
   if( string && string[0] ){
      int i=0;
      while( string[i] ){
         if( string[i]!='0' && string[i]!='1' && string[i]!='2' &&
             string[i]!='3' && string[i]!='4' && string[i]!='5' &&
             string[i]!='6' && string[i]!='7' && string[i]!='8' &&
             string[i]!='9' && string[i]!='.' ){
            return 0;
         }
         i++;
      }
      return 1;
   }
   return 0; 
}*/

/*int is_digit( char znak )
{
   if( znak=='0' || znak=='1' || znak=='2' || znak=='3' || znak=='4' ||
       znak=='5' || znak=='6' || znak=='7' || znak=='8' || znak=='9' ){
      return 1;
   }
   return 0;
}*/



const char* change_ext(const char* name,const char* new_ext,string& out)
{
   string tmp_file;
   getbasename_new(name,tmp_file);
   
   out=tmp_file.c_str();
   out += ".";
   out += new_ext;
   
   return out.c_str();
}

const char* add_postfix(const char* name,const char* new_ext,string& out)
{
   string tmp_file;
   getbasename_new(name,tmp_file);
   
   out=tmp_file.c_str();
   out += new_ext;
   
   return out.c_str();
}



const char* getbasename_new(const char* name,string& out)
{
  char outtmp[1024];
  int i=0;
  out.clear();
  memset(outtmp,'\0',1024);

/*  while(name[i]!='.' && name[i]!='\0'){
    out += name[i];   
    i++;
  }
  out += '\0';*/

  i = strlen(name)-1;
  while( name[i]!='.' && i>=0 ){
    i--;  
  }
  if( i>=0 ){
    strncpy(outtmp,name,i);    
    out=outtmp;
  }else{
    out=out.c_str();
  }

  return out.c_str();
}                        

int bg_does_file_exists(const char* fname)
{
   int bRet=0;
   if(fname && fname[0] ){
//      string szFile=fname;
//      szFile.env2str();

      if( access( fname, F_OK ) == 0 ) {
         bRet = 1;
      }
   }
   return bRet;
}



int bg_read_list( const char* file, vector<string>& out_list )
{
 if( !bg_does_file_exists(file) ){
  printf("ERROR : file %s does not exist !\n",file);
  return -1;
 } 
  
 ifstream ifile( file );
 string line;

 while ( ifile.good() ){
    getline(ifile,line);
    if( mystring::get_first_non_white( line.c_str() )=='#' ){
       printf("Commented line skipped: %s\n",line.c_str());
       continue;
    }
    
    if( strlen(line.c_str()) ){
      out_list.push_back(line);
    }
 }
         
 
 return out_list.size();
}


double deg2rad(double in_deg )
{
 return (in_deg/180.00)*M_PI;
}

double rad2deg(double in_rad)
{
 return (in_rad/M_PI)*180.00;
}

double calc_rms( double* ftab, long cnt )
{
   double sum=0.00,sum2=0.00;
   
   for(int i=0;i<cnt;i++){
      sum += ftab[i];
      sum2 += (ftab[i]*ftab[i]);
   }
   
   double mean2 = (sum2/cnt);
   double avg2  = (sum/cnt)*(sum/cnt);
   
    double rms = sqrt( mean2 - avg2 );
    return rms;
}

double calc_avg_rms( double* ftab, long cnt, double& rms, int bOnlyPositive /*=0*/ )
{
   double sum=0.00,sum2=0.00;
   int count=0;
   
   for(int i=0;i<cnt;i++){
      if( bOnlyPositive<=0 || ftab[i] > 0 ){
         sum += ftab[i];
         sum2 += (ftab[i]*ftab[i]);
         count++;
      }
   }
   
   double mean2 = (sum2/count);
   double avg2  = (sum/count)*(sum/count);
   
   rms = sqrt( mean2 - avg2 );
   return (sum/count);
}


void my_sort_float( double* ftab, long cnt )
{
	double divider = ftab[0];
	
	int beg = 1;
	int end = cnt-1;
	double tmp;

	if(cnt){	
		while(beg<=end){
			if(ftab[beg]>divider){
				if(ftab[end]<=divider){
					REPLACE_ELEMS( ftab, beg, end )
					beg++;
					end--;
				}else{
					end--;
				}
			}else{		
				beg++;
				if(ftab[end]>divider)
					end--;
			}
		}
		if(end!=0){
			REPLACE_ELEMS( ftab, end, 0)
		}

		my_sort_float( ftab, end );
		my_sort_float( &(ftab[beg]), cnt-beg );
	}

}


int freq2channel( double freq )
{
   int ch = (int)(freq*(4096/480.00));   
   return ch;
}
            

double ch2freq(int channel)
{
        double ch_width = ((double)BG_MAX_FREQ)/((double)BG_CHANNELS);
        return (ch_width*channel);
        
}

int mkdir(const char* path)
{
        char szCmd[1024];
        sprintf(szCmd,"mkdir -p %s",path);
        int ret = system(szCmd);
        
        return ret;        
}

double mysqr(double x)
{
   return (x*x);
}                 

const char* getbasename_new(const string& name,string& out)
{
   int i=0;
   out.clear();
   while(name[i]!='.' && name[i]!='\0'){
      out += name[i];
      i++;
   }
   out += '\0';
   return out.c_str();
}

                              
                              
                                                                                                      
int ParseCommaList( char* szList, vector<string>& out_list, const char * sep )
{
   char* saveptr=NULL;
   char* str1=szList;   
   
   for (str1 = szList; ; str1 = NULL) {
      char* token = strtok_r(str1, sep, &saveptr);
      if (token == NULL)
         break;

      // trim right :         
      char szToken[64];
      strcpy(szToken,token);
      int i=strlen(szToken)-1;
      while(szToken[i]==' '){
         szToken[i]='\0';
         i--;
      }   
         
      out_list.push_back(string(szToken));
   }
      
   return out_list.size();
}
                                                                                                      

                                                                                                      
/*double interpolate( double x, double x1, double y1, double x2, double y2 )
{ 
   double delta = ((y2-y1)/(x2-x1))*(x-x1);
   double interpol_val = y1 + delta;
   return interpol_val;
} */  

int find_interpol_values( vector<cValue>& list, double x, double& x1, double& y1, double& x2, double& y2, int bRe )
{
   cValue prev_val = list[0];
   
   if( x < prev_val.x ){
      x1 = x2 = prev_val.x;
      if( bRe > 0 ){
         y1 = y2 = prev_val.y;
      }else{
         y1 = y2 = prev_val.z;
      }
      return 0;
   }
   if( x > list[list.size()-1].x ){
      cValue& last = list[list.size()-1];
      x1 = x2 = last.x;
      if( bRe > 0 ){
         y1 = y2 = last.y;
      }else{
         y1 = y2 = last.z;
      }
      return 0;                  
   }
   
   for(int i=0;i<list.size();i++){
      cValue curr_val = list[i];
      if( prev_val.x<=x && x<=curr_val.x ){
         x1 = prev_val.x;
         x2 = curr_val.x;
         
         if( bRe > 0 ){
            y1 = prev_val.y;
            y2 = curr_val.y;
         }else{
            y1 = prev_val.z;
            y2 = curr_val.z;                        
         }
         return 1;
      }
 
      prev_val = list[i];
   }

   return 0;
}

int find_interpol_values( vector<cValue>& list, double x, cValue& prev, cValue& after )
{
   cValue prev_val = list[0];
   
   if( x < prev_val.x ){
      prev = prev_val;
      after = prev_val;
      return 0;
   }
   if( x > list[list.size()-1].x ){      
      cValue& last = list[list.size()-1];
      prev = last;
      after = last;
      return 0;                  
   }
   
   for(int i=0;i<list.size();i++){
      cValue curr_val = list[i];
      if( prev_val.x<=x && x<=curr_val.x ){
         prev = prev_val;
         after = curr_val;
         return 1;
      }
 
      prev_val = list[i];
   }

   return 0;
}


double interpolate( vector<cValue>& file_values, double x, int bRe )
{
   double x1,y1,x2,y2;
   int n_find = find_interpol_values( file_values, x, x1, y1, x2, y2, bRe );
   double val = interpolate( x, x1, y1, x2, y2 );
   
   return val;
}   
            
double power(double x,double y)
{
   return exp(y*log(x));
}                    
                                                                                                                                                      
double db2num(double in_db)
{
   double num_val = power(10.00,(in_db/10.00));
   return num_val;
}
                                                                                                                                                      
                                                                                                                                                      

int read_file(const char* file,vector<cValue>& out_list, int db2num, int x_col, int y_col, int min_col )
{
   MyFile infile(file);
   const char* pLine=NULL;
   out_list.clear();

   while( pLine = infile.GetLine( TRUE ) ){
      if( mystring::get_first_non_white( pLine )=='#' ){
//         printf("Line skipped : %s\n",pLine);
         continue;
      }
                        
      MyParser pars=pLine;
      CMyStrTable items;
      pars.GetItems( items );
      
      int max_col = MAX(x_col,y_col);

      if( items.size() > 0 ){
         cValue tmp;
         
         if( items.size() > x_col ){
            tmp.x = atof(items[x_col].c_str());
         }

         if( items.size() > y_col ){
            tmp.y = atof(items[y_col].c_str());
         }

         if( items.size() > (y_col+1) ){
            tmp.z = atof(items[y_col+1].c_str());
         }
         if( items.size() > (y_col+2) ){
            tmp.v = atof(items[y_col+2].c_str());
         }
         
         if( db2num > 0 ){
            tmp.y = exp( log(10.00) * tmp.y * 0.1 );
         }

         if( items.size() >= min_col ){    
            out_list.push_back(tmp);
         }
      }
   }
   printf("Read %d values from file %s\n",(int)out_list.size(),file);
   
   return out_list.size();
}                   

int read_file(const char* file,vector<cIntRange>& out_list, int x_col, int y_col, int bDoClear )
{
   MyFile infile(file);
   const char* pLine=NULL;

   if( bDoClear > 0 )
      out_list.clear();

   while( pLine = infile.GetLine( TRUE ) ){
      if( mystring::get_first_non_white( pLine )=='#' )
         continue;
                        
      MyParser pars=pLine;
      CMyStrTable items;
      pars.GetItems( items );
      
      int max_col = MAX(x_col,y_col);

      if( items.size() > 0 ){
         cIntRange tmp;
         
         if( items.size() > x_col ){
            tmp.start_int = atol(items[x_col].c_str());
         }

         tmp.end_int=tmp.start_int;
         if( items.size() > y_col ){
            tmp.end_int = atol(items[y_col].c_str());
         }

         out_list.push_back(tmp);
      }
   }
   printf("Read %d values from file %s\n",(int)out_list.size(),file);
   
   return out_list.size();
}                   


void alloc(vector<cValue>& arr, int _size )
{
   arr.clear();
   arr.reserve(_size);

   cValue tmp;
   for(int i=0;i<_size;i++){
      arr.push_back(tmp);
   }
                  
}

void init_cvalue_vector(vector<cValue>& vec, int size)
{
//   vec.clear();
   vec.reserve(size);
   cValue val;
   val.x = 0.00;
   val.y = 0.00;
   
   for(int i=0;i<size;i++){
      vec.push_back( val );
   }   
}

void print_cmdline(int argc, char * argv[])
{
   for(int i=0;i<argc;i++){
      printf("%s ",argv[i]);
   }
   printf("\n");
}


int find_value( vector<int>& arr, int value )
{
   for(int i=0;i<arr.size();i++){
      if( arr[i] == value ){
         return i;
      }
   }
   
   return -1;
}
