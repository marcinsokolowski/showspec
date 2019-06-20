/*
 **************************************************************************\
 **  transientlib : library for identification of short optical flashes 
 **						with the wide field cameras.
 **  This software was written by Marcin Sokolowski ( msok@fuw.edu.pl ) 
 **	it was a substantial part of analysis performed for PHD thesis 
 **  ( Investigation of astrophysical phenomena in short time scales with "Pi of the Sky" apparatus )
 **	it can be used under the terms of GNU General Public License.
 **	In case this software is used for scientific purposes and results are
 **	published, please refer to the PHD thesis submited to astro-ph :
 **
 **		http://arxiv.org/abs/0810.1179
 **
 ** Public distribution was started on 2008-10-31
 **
 ** 
 ** NOTE : some of the files (C files) were created by other developers and 
 **        they maybe distributed under different conditions.
 ** 

 ******************************************************************************
 ** This program is free software; you can redistribute it and/or modify it
 ** under the terms of the GNU General Public License as published by the
 ** Free Software Foundation; either version 2 of the License or any later
 ** version. 
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 ** General Public License for more details. 
 **
 *\**************************************************************************

*/           
#include <stdio.h>
#include "mystring.h"
#include "myutil.h"
#include <sys/utsname.h>
#include "myparser.h"
#include "basestructs.h"
#include "mystrtable.h"

double MAX_func( double x, double y )
{
	if(x>y)
		return x;
	else 
		return y;
}

double MIN_func( double x, double y )
{
   if(x<y)
      return x;
   else
      return y;
}


void my_insertion_sort( LONG_T* tab, LONG_T cnt )
{
	for(int i=0;i<cnt;i++){
		LONG_T minim = tab[i];
		LONG_T j_min = i;
		for(int j=i+1;j<cnt;j++){
			if(tab[j]<minim){
				minim = tab[j];
				j_min = j;
			}				
		}
		if(j_min!=i){
			LONG_T tmp;
			tmp = minim;
			tab[j_min] = tab[i];
			tab[i] = minim;
		}
	}
}


#define REPLACE_ELEMS( tab, pos1, pos2 ) { tmp=tab[pos1]; tab[pos1]=tab[pos2]; tab[pos2]=tmp; } 

int get_median( LONG_T* tab, LONG_T cnt )
{
	int pos = (cnt/2);
	my_qsort( tab , cnt );
	return tab[pos];
}


int get_no_stars_median( LONG_T* tab, LONG_T cnt )
{
	int end_cnt = (int)(0.75*cnt);
	my_qsort( tab , cnt );
	int pos = (end_cnt/2);
	return tab[pos];
}



// to sort table call : my_qsort( tab, 10, 10 )
void my_qsort( LONG_T* tab, LONG_T cnt )
{
	register LONG_T divider = tab[0];
	
	register int beg = 1;
	register int end = cnt-1;
	LONG_T tmp;

	if(cnt){	
		while(beg<=end){
			if(tab[beg]>divider){
				if(tab[end]<=divider){
					REPLACE_ELEMS( tab, beg, end )
					beg++;
					end--;
				}else{
					end--;
				}
			}else{		
				beg++;
				if(tab[end]>divider)
					end--;
			}
		}
		if(end!=0){
			REPLACE_ELEMS( tab, end, 0)
		}

		my_qsort( tab, end );
		my_qsort( &(tab[beg]), cnt-beg );
	}
} 

// to sort table call : my_qsort( tab, 10, 10 )
void my_qsort_int( int* tab, LONG_T cnt )
{
	register LONG_T divider = tab[0];
	
	register int beg = 1;
	register int end = cnt-1;
	LONG_T tmp;

	if(cnt){	
		while(beg<=end){
			if(tab[beg]>divider){
				if(tab[end]<=divider){
					REPLACE_ELEMS( tab, beg, end )
					beg++;
					end--;
				}else{
					end--;
				}
			}else{		
				beg++;
				if(tab[end]>divider)
					end--;
			}
		}
		if(end!=0){
			REPLACE_ELEMS( tab, end, 0)
		}

		my_qsort_int( tab, end );
		my_qsort_int( &(tab[beg]), cnt-beg );
	}
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


BOOL_T my_sorted_table_overlap_check( LONG_T* tab1, LONG_T cnt1,
                                      LONG_T* tab2, LONG_T cnt2 )
{
	long i=0;
	long j=0;
	
	while(i<cnt1 && j<cnt2){
		if(tab1[i]==tab2[j])
			return TRUE;
		if(tab1[i]<tab2[j])
			i++;
		else
			j++;
	}
	return FALSE;
}

BOOL_T my_table_overlap_check( LONG_T* tab1, LONG_T cnt1,
                               LONG_T* tab2, LONG_T cnt2 )
{
	for(int i=0;i<cnt1;i++){
		for(int j=0;j<cnt2;j++){	
			if(tab1[i]==tab2[j])
				return TRUE;
		}
	}
	return FALSE;
}

LONG_T my_find_max_value_long( LONG_T* table, LONG_T cnt )
{
	LONG_T max=-100000;
	for(register int i=0;i<cnt;i++){
		if(table[i]>max)
			max = table[i];
	}
	return max;
}

LONGLONG_T my_find_max_value( LONGLONG_T* table, LONG_T cnt )
{
	LONGLONG_T max=-1000000;
	for(int i=0;i<cnt;i++){
		if(table[i]>max)
			max = table[i];
	}
	return max;
}

BOOL_T check_sort( LONG_T* tab, LONG_T size )
{
	LONG_T prev=tab[0];
	for(int i=0;i<size;i++){
		if(tab[i]<prev){
			printf("Table sorted incorrectly !!!\n");	
			printf("tab[%d]=%d < tab[%d]=%d\n",i,tab[i],i-1,prev);
			return FALSE;
		}
		prev = tab[i];
	}
	return TRUE;
}

void DumpTable( LONG_T* tab, LONG_T size )
{
	mystring szList;
	for(int i=0;i<size;i++){
		szList << tab[i] << ",";
	}
	printf("%s\n",szList.c_str() );
}

/*LONG_T my_round(double x)
{
	register int x_l = (int)floor(x);

	register double h = (x-x_l);
	if(h>0.5)
		return (x_l+1);
	else
		return (x_l);
}*/


LONG_T long_sum( LONG_T* tab, LONG_T cnt, LONG_T start_from /*=0*/ )
{
	register long sum=0;
	for(register long i=start_from;i<cnt;i++){
		sum += tab[i];
	}	
	return sum;
}

int find_value( LONG_T* tab, int tab_size, LONG_T val )
{
	for(register int i=0;i<tab_size;i++){
		if(tab[i] == val)
			return i;
	}
	return NOT_FOUND;
}

int qfind_value( LONG_T* tab, int tab_size, LONG_T val )
{
	register int pos0=0;
	register int pos1=(tab_size-1);
	while( pos1>pos0 && pos1>=0 && pos0>=0 ){
		int pos_new = (pos1+pos0)/2;
		if( tab[pos_new]==val ){
			while( pos_new>=0 && tab[pos_new]==val ){
				pos_new--;
			}
			return (pos_new+1);
		}else{
			if( val < tab[pos_new] ){
				pos1 = pos_new-1;
			}else{
				pos0 = pos_new+1;
			}			
		}
	}
	if( pos1>=0 && pos0>=0 && pos1<tab_size && pos0<tab_size ){
		if( tab[pos1]==val ){
			while( pos1>=0 && tab[pos1]==val ){
				pos1--;
			}
			return (pos1+1);
		}
	}
	return NOT_FOUND;

}

int find_value( vector<int>& tab, int value )
{
	int pos=0;
	for( vector<int>::iterator i=tab.begin();i!=tab.end();i++,pos++){
		if( (*i)==value )
			return pos;
	}
	return -1;
}

int qfind_value( vector<int>& tab, int val )
{
	int tab_size=tab.size();
	register int pos0=0;
	register int pos1=(tab_size-1);
	while( pos1>pos0 && pos1>=0 && pos0>=0 ){
		int pos_new = (pos1+pos0)/2;
		if( tab[pos_new]==val ){
			while( pos_new>=0 && tab[pos_new]==val ){
				pos_new--;
			}
			return (pos_new+1);
		}else{
			if( val < tab[pos_new] ){
				pos1 = pos_new-1;
			}else{
				pos0 = pos_new+1;
			}			
		}
	}
	if( pos1>=0 && pos0>=0 && pos1<tab_size && pos0<tab_size ){
		if( tab[pos1]==val ){
			while( pos1>=0 && tab[pos1]==val ){
				pos1--;
			}
			return (pos1+1);
		}
	}
	return NOT_FOUND;

}


int do_print_usage( int argc, char* argv[], int min_argc_req )
{
	if( argc<min_argc_req || ( strncmp(argv[1],"-h",2)==0 || strncmp(argv[1],"--h",2)==0 ) ){
		return 1;
	}
	return 0;
}

int is_default( const char* param )
{
	return (strlen(param)==0 || strcmp(param,"-")==0);
}


int find_string( const char** tab, const char* szValue )
{
	register int i=0;
	while(tab[i]){
		if(strcmp(tab[i],szValue)==0)
			return i;
		i++;
	}	
	return -1;
}


mystring get_uname()
{
	struct utsname info;
	uname( &info );
	mystring ret = info.nodename;
	return ret;	
}

int find_min_value( int* tab, int count, int& min_pos )
{
	register int min_value=5000000;
	min_pos=-1;
	for(register int i=0;i<count;i++){
		if(tab[i]<min_value){
			min_value=tab[i];
			min_pos = i;
		}
	}
	return min_value;
}

BOOL_T compare_double( double l, double r )
{
	char szLeft[20],szRight[20];
	sprintf(szLeft,"%.8f",l);
	sprintf(szRight,"%.8f",r);
	int ret = strcmp( szLeft, szRight );
	if(ret==0)
		return TRUE;
	return FALSE;		
}

const char* safe_string( const char* szValue )
{
	if( szValue )
		return szValue;
	return "";
}

int safe_atol( const char* szValue )
{
	if(szValue && szValue[0])
		return atol( szValue );
	return 0;
}

double safe_atof( const char* szValue )
{
	if(szValue && szValue[0])
		return atof( szValue );
	return 0;
}


mystring get_param_value( const char* szParam )
{
	mystring szName,szValue;
	if( szParam && szParam[0] ){
		MyParser pars = szParam;
		pars.GetVarAndValue( szName, szValue );


		// NEW - 20050419 due to f.g : -where="idaynight=20050404"
		char* ptr_eq = (char*)strstr( szParam, "=" );
		if( ptr_eq && ptr_eq[0] ){
			szValue = ptr_eq+1;
		}		
	}
	return szValue;
}

BOOL_T check_param( const char* param, const char* name, double& fValue )
{
	if( strncmp( param, name, strlen(name))==0){
		mystring szValue = get_param_value( param );
		fValue = atof(szValue.c_str());
		return TRUE;
	}
	return FALSE;
}

BOOL_T check_param( const char* param, const char* name, float& fValue )
{
	if( strncmp( param, name, strlen(name))==0){
		mystring szValue = get_param_value( param );
		fValue = atof(szValue.c_str());
		return TRUE;
	}
	return FALSE;
}


BOOL_T check_param( const char* param, const char* name, int& nValue )
{
	if( strncmp( param, name, strlen(name))==0){
		mystring szValue = get_param_value( param );
		nValue = atol(szValue.c_str());
		return TRUE;
	}
	return FALSE;
}


BOOL_T check_param( const char* param, const char* name, mystring& szValue )
{
	if( strncmp( param, name, strlen(name))==0){
		szValue = get_param_value( param );
		return TRUE;
	}
	return FALSE;
}

void my_printf_now( const char* msg )
{
	printf("%s",msg);
	fflush(0);
}

int exec_cmd( const char* szCmd, char* szOutput, int size )
{
//	memset(szOutput,'\0',size);
	FILE* a = popen(szCmd, "r");
	int n_read = fread(szOutput, size, 1, a);
   pclose(a);
	if( n_read>size )
		return -1;
	return n_read;	
}

const char* GetOK( BOOL_T bOK )
{
	if( bOK )
		return "OK";
	return "FAILED";
}

double find_max_value( double* values, int cnt, int& pos )
{
	double ret=-1000000.00;
	for(register int i=0;i<cnt;i++){
		if(values[i]>ret){
			ret = values[i];
			pos = i;
		}
	}
	return ret;
}

double calc_average_rms( double* values, int cnt, double& rms )
{
	double sum =0.00;
	double avg2=0.00;
	for(register int i=0;i<cnt;i++){
		sum += values[i];
		avg2 += (values[i]*values[i]);
	}
	double avg = (sum/cnt);
	rms = sqrt( avg2/cnt - avg*avg );
	return avg;
}

double calc_rms( int mean, LONG_T* values, int cnt )
{
	double sum=0.00;
	for(register int i=0;i<cnt;i++){
		sum += (values[i]-mean)*(values[i]-mean);
	}
	double ret = sqrt(sum/cnt);

	return ret;
}


const char* my_stristr( const char* in_string, const char* substr )
{
	mystring szTMP = in_string;
	mystring szLONG = szTMP.getupper();

	mystring szSHORT_tmp = substr;
	mystring szSHORT = szSHORT_tmp.getupper();


	const char* ptr = strstr( szLONG.c_str(), szSHORT.c_str() );
	if( ptr ){
		int pos = ( ptr - szLONG.c_str() );
		return ( in_string + pos );
	}

	return NULL;
}

void CopyCluster( LONG_T* tab_out, LONG_T* tab_in, int cnt )
{
	for(register int i=0;i<cnt;i++){
		tab_out[i] = tab_in[i];
	}
}

void print_line()
{
	printf("--------------------------------------\n");fflush(0);
}

void dump_table( LONG_T* tab ,int cnt )
{
   for(int i=0;i<cnt;i++){
      printf("%d ",tab[i]);
   }
   printf("\n");
}

int check_zero( int value )
{
	if( value<0 )
		return 0;
	return value;
}

float non_zero( float val, float val_default )
{
	if( val>0 )
		return val;
	return val_default;
}


void my_calc_mean_and_sigma( int* tab, int count, double& mean, double& rms,
										int min_val, int max_val )
{
	double sum=0.00;
	double sum2=0.00;
	int used_cnt=0;

	for(int i=0;i<count;i++){
		if( tab[i] >= min_val && tab[i] <= max_val ){
         double newval = ((double)tab[i]);
         double newval2 = (newval*newval);
		
			sum = sum + newval;
			sum2 = sum2 + newval2;
			used_cnt++;
		}
	}
	
	mean = ((double)sum) / ((double)used_cnt);
	double mean2 = ((double)sum2) / ((double)used_cnt);
	double avg2 = mean*mean;
	if( mean2 >= avg2 ){
		rms = sqrt( mean2 - avg2 );
	}else{
		double sigma=0.00;
		for(int i=0;i<count;i++){
	      if( tab[i] >= min_val && tab[i] <= max_val ){
				double diff = (tab[i]-mean);
				sigma = sigma + diff*diff;
			}
		}

		rms =  sqrt( sigma / used_cnt );
		printf("WARNING : mean2=%.4f < avg2=%.4f -> rms=%.4f\n",mean2,avg2,rms);
	}	
}


void my_calc_mean_and_sigma_elem( ELEM_TYPE* tab, int count, double& mean, double& rms,
										int min_val, int max_val )
{
	double sum=0.00;
	double sum2=0.00;
	int used_cnt=0;

	for(int i=0;i<count;i++){
		if( tab[i] >= min_val && tab[i] <= max_val ){
         double newval = ((double)tab[i]);
         		
			sum = sum + newval;
			sum2 = sum2 + (newval*newval);
			used_cnt++;
		}
	}
	
	mean = ((double)sum) / ((double)used_cnt);
	double mean2 = ((double)sum2) / ((double)used_cnt);
	double avg2 = mean*mean;
	if( mean2 >= avg2 ){
		rms = sqrt( mean2 - avg2 );
	}else{
		double sigma=0.00;
		for(int i=0;i<count;i++){
	      if( tab[i] >= min_val && tab[i] <= max_val ){
				double diff = (tab[i]-mean);
				sigma = sigma + diff*diff;
			}
		}

		rms =  sqrt( sigma / used_cnt );
		printf("WARNING : mean2=%.4f < avg2=%.4f -> rms=%.4f\n",mean2,avg2,rms);
	}	
}


BOOL_T ParseCommand( const char* szCmd, sCommand& cmd )
{
	MyParser pars=szCmd;
	CMyStrTable items;  
	if( pars.GetItems( items, "," ) ){
		cmd.command = atol( items[0].c_str() );
		if( items.size() >= 2 ){
			cmd.szParam1 = items[1].c_str();
		}
		if( items.size() >= 3 ){
			cmd.szParam2 = items[2].c_str();
		}
	}   
                                                            
   return ( items.size()>0 );
}


                                                               