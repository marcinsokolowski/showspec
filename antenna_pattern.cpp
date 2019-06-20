 #include "antenna_pattern.h"
#include "myparser.h"
#include "myfile.h"
#include "mystrtable.h"
#include <math.h>

#include "bg_globals.h"
#include "mwa_beam_interface.h"

#ifndef _NO_ROOT_
 #include <TGraph2D.h>
#endif

#include "showspec_math.h"
#include "skyinfo.h"

#define DEFAULT_PATTERN_VALUE -10000

int MAX_THETA_COUNT=9001;
int MAX_PHI_COUNT=36001;

double UserDefinedPatternFormula( double freq_mhz, double phi_deg, double theta_deg );

int CAntPatternParser::m_BradleyVersion=0;
int CAntPatternParser::m_ElevationInAntPattern=0;
int CAntPatternParser::m_AnglesInRadians=0;
int CAntPatternParser::m_GainInLinearScale=0;
int CAntPatternParser::m_InputFileFormat_GainColIdx=8;
int CAntPatternParser::m_InpitFileFormat_ThetaColIdx=0;
int CAntPatternParser::m_InpitFileFormat_PhiColIdx=1;  
  

CFreqPattern::CFreqPattern( double _freq, int theta_samples, int phi_samples )
   : freq(_freq), m_ThetaSamples(theta_samples), m_PhiSamples(phi_samples), pRootGraph(NULL)
   {
      // ant_pattern[MAX_THETA_COUNT][MAX_PHI_COUNT];
      ant_pattern = new double*[MAX_THETA_COUNT];
      for(int i=0;i<MAX_THETA_COUNT;i++){
         ant_pattern[i] = new double[MAX_PHI_COUNT];
      }
   
      for(int t=0;t<m_ThetaSamples;t++){
         for(int p=0;p<m_PhiSamples;p++){
            ant_pattern[t][p] = DEFAULT_PATTERN_VALUE;
         }
      }
   }


CFreqPattern::CFreqPattern( const CFreqPattern& right )
{
   (*this) = right;
}

CFreqPattern::~CFreqPattern()
{
#ifndef _NO_ROOT_
   if( pRootGraph ){   
      delete pRootGraph;
   }
#endif

  if( ant_pattern ){
     for(int i=0;i<MAX_THETA_COUNT;i++){
        delete [] ant_pattern[i];
     }
  }
}

CFreqPattern& CFreqPattern::operator=( const CFreqPattern& right )
{
   freq = right.freq; // MHz
   m_ThetaSamples = right.m_ThetaSamples;
   m_PhiSamples = right.m_PhiSamples;


//   memcpy(ant_pattern,right.ant_pattern, sizeof(ant_pattern) );
   ant_pattern = new double*[MAX_THETA_COUNT];
   for(int i=0;i<MAX_THETA_COUNT;i++){
      ant_pattern[i] = new double[MAX_PHI_COUNT];
      
      for(int x=0;x<MAX_PHI_COUNT;x++){
         ant_pattern[i][x] = right.ant_pattern[i][x];
      }
   }
                        


   if( gBGPrintfLevel>2 ){
      printf("DEBUG : sizeof = %d ?= %d x %d x %d, [0][0]=%.8f vs %.8f\n",(int)sizeof(ant_pattern),MAX_THETA_COUNT,MAX_PHI_COUNT,(int)sizeof(double),ant_pattern[0][0],right.ant_pattern[0][0]);
   }
   
   pRootGraph = NULL;
#ifndef _NO_ROOT_
   if( right.pRootGraph ){
      pRootGraph = new TGraph2D(*right.pRootGraph);
   }
#endif
}
                   

int CAntPatternParser::m_bEnableInterpolation=0;
double CAntPatternParser::m_SimulMaxTheta=180.00;
double CAntPatternParser::m_SimulMaxPhi=360.00;  
int CAntPatternParser::m_AutoDetectTheta=0;  

CAntPatternParser::CAntPatternParser()
{
}


CAntPatternParser::~CAntPatternParser()
{
}

double CAntPatternParser::AutoCheckMaxTheta(const char* pattern_file)
{
   MyIFile in( pattern_file );
   const char* pLine=NULL;
   int nFreq=0;
   double max_theta=-1000;
   
   while( pLine = in.GetLine( TRUE ) ){
      MyParser pars=pLine;
      CMyStrTable items;
      pars.GetItems( items );
                         
                         
      if( mystring::get_first_non_white( pLine )=='#' )
      {
         if( strstr(pLine,"Frequency") ){
            if( nFreq>0 ){
               break;
            }
            nFreq++;
         }
      }else{
         double theta = atof(items[0].c_str());
         if( theta > max_theta ){
            max_theta = theta;
         }
      }
   }
   printf("Maximum theta angle auto-detected = %.4f [deg]\n",max_theta);
   
   return max_theta;
}
    
int CAntPatternParser::Read(const char* pattern_file, int n_phi_samples, int n_theta_samples, double param_freq)
{
   printf("DEBUG : CAntPatternParser::Read(%s,%d,%d,%.2f) , \n",pattern_file,n_phi_samples,n_theta_samples,param_freq);   

   if( m_AutoDetectTheta > 0 ){
      m_SimulMaxTheta = AutoCheckMaxTheta(pattern_file);
   }

	MyIFile in( pattern_file );
	CFreqPattern current_freq_rec(param_freq,n_theta_samples,n_phi_samples);

	int bNewFreqStarted=0;
	int freq_count=0;
	double sum_gain=0.00;
	double min_gain=10000000.00;
	double max_gain=-10000000.00;
	
   int prev_phi_index=-1;
   int prev_theta_index=-1;	
   int prev_line = -1;
   mystring szPrevLine;
   double phi_step = -1.00,theta_step = -1.00;

   if( m_BradleyVersion > 0 ){
      // In this version pattern file header information is not required and is passed in the command line:
      bNewFreqStarted=1;
   }
	
   const char* pLine=NULL;
   int line_counter=0;
   int shown_res_info=0;
   int header_lines=0;
   
   while( pLine = in.GetLine( TRUE ) ){
      MyParser pars=pLine;
      CMyStrTable items;
      pars.GetItems( items );
                  

      if( mystring::get_first_non_white( pLine )=='#' )
      {
         if( strstr(pLine,"Frequency") ){
            if( items.GetCount() >= 2 ){
               current_freq_rec.freq = atof(items[1].c_str())/1e6;
            }else{
               printf("ERROR : could not parse frequency line |%s|, at least 2 items required !!!\n",pLine);            
            }
                     
            bNewFreqStarted=0;
            freq_count=0;
            sum_gain=0.00;
            min_gain=10000000.00;
            max_gain=-10000000.00;               
            continue;
         }
         if( strstr(pLine,"#No. of Theta Samples:") ){
            if( items.GetCount() >= 5 ){
               current_freq_rec.m_ThetaSamples = atol(items[4].c_str());            
               printf("DEBUG : current_freq_rec.m_ThetaSamples := %d\n",current_freq_rec.m_ThetaSamples);
            }else{
               printf("ERROR : could not parse theta samples line |%s|, at least 5 items required !!!\n",pLine);
            }
         }         
         if( strstr(pLine,"#No. of Phi Samples:") ){
            if( items.GetCount() >= 5 ){
               current_freq_rec.m_PhiSamples = atol(items[4].c_str());                                   
               printf("DEBUG : current_freq_rec.m_PhiSamples := %d\n",current_freq_rec.m_PhiSamples);                     
               bNewFreqStarted=1;
               if( gBGPrintfLevel > 1 ){
                  printf("------------------------------------------------------------------------- %.2f [MHz] -------------------------------------------------------------------------\n",current_freq_rec.freq);
                  printf("New frequency %.2f started with %d phi and %d theta samples started\n",current_freq_rec.freq,current_freq_rec.m_PhiSamples,current_freq_rec.m_ThetaSamples);
               }
            }else{
               if( gBGPrintfLevel > 1 ){
                  printf("ERROR : could not parse phi samples line |%s|, at least 5 items required !!!\n",pLine);
               }
            }
         }
         
         header_lines++;
      }

      phi_step = m_SimulMaxPhi / (current_freq_rec.m_PhiSamples-1);
      theta_step = m_SimulMaxTheta / (current_freq_rec.m_ThetaSamples-1);                     
      
      if( shown_res_info <= 0 && bNewFreqStarted ){
         printf("Resolutions of input file are : phi_step = %.4f = %.2f / %d, theta_step = %.4f = %.4f / %d\n",phi_step,m_SimulMaxPhi,(current_freq_rec.m_PhiSamples-1),theta_step,m_SimulMaxTheta,(current_freq_rec.m_ThetaSamples-1));fflush(stdout);
         shown_res_info = 1;
      }      
                                                                                            
      if( items.size() > 0 ){
         // parse line 
         const char* szFirstItem = items[0].c_str();
//       printf("bNewFreqStarted=%d : %s -> %c\n",bNewFreqStarted,szFirstItem,szFirstItem[0]);
         if( bNewFreqStarted > 0 ){
            if( szFirstItem[0] != '#' && szFirstItem[0] != '*' ){
               double theta_deg = atof(items[CAntPatternParser::m_InpitFileFormat_ThetaColIdx].c_str());
               double phi_deg = atof(items[CAntPatternParser::m_InpitFileFormat_PhiColIdx].c_str());
               double gain_db = atof(items[CAntPatternParser::m_InputFileFormat_GainColIdx].c_str()); // was just 8 (not can be defined)
               double gain = db2num(gain_db);
               
//               if( m_GainInLinearScale > 0 ){
                  // gain not in dB (FEKO) but linear - Bradley                   
//                  gain = gain_db;
//                  gain_db = 10.00 * log10( gain );
//               }
               
               if( m_BradleyVersion > 0 ){
                  // Bradley Mayers gives :
                  // ELEVATION AZIMUTH GAIN_LINEAR[NOT DB] :
                  if( m_AnglesInRadians > 0 ){
                     phi_deg = phi_deg * (180.00/M_PI);
                     theta_deg = theta_deg * (180.00/M_PI);
                  }

                  if( m_ElevationInAntPattern > 0 ){
                     theta_deg = (90.00 - theta_deg);
                  }
                  if( theta_deg < 0 ){
                     printf("WARNING : theta_deg = %.20f < 0 -> forced to 0  !\n",theta_deg);
                     theta_deg = 0.00;
                  }
                  
                  phi_deg = ( 90.00 - phi_deg);
                  if( phi_deg < 0 ){
                     phi_deg = phi_deg + 360.00;
                  }
                  if( phi_deg > 360.00 ){
                     phi_deg = phi_deg - 360.00;
                  }
                  
               }

               if( m_GainInLinearScale > 0 ){                  
                  gain = gain_db;
                  if( gain > 0 ){
                     gain_db = 10.00 * log10( gain );
                  }else{
                     gain_db = -1000.00;
                  }
               }

//            double gain  = gain_db; // TEST ONLY 

               phi_step = m_SimulMaxPhi / (current_freq_rec.m_PhiSamples-1);
               theta_step = m_SimulMaxTheta / (current_freq_rec.m_ThetaSamples-1);
               int phi_index = floor(phi_deg/phi_step);
               int theta_index = floor(theta_deg/theta_step);
               
               if( phi_step < 0.9 ){
                  phi_index = round( phi_deg * round(1.00/phi_step) );
               }
               if( theta_step < 0.9 ){
                  theta_index = round( theta_deg * round(1.00/theta_step) );
               }
               
//               if( gDebug>2 ){
//                   printf("Read: (phi_deg,theta_deg) = (%.4f,%.4f) [deg] -> indexes (%d,%d)\n",phi_deg,theta_deg,phi_index,theta_index);
//               }
               
//               if( phi_index<=1630 && phi_index>=1620 && theta_index==8 ){
//                  printf("odo");
//               }
               
               if( phi_index == prev_phi_index && theta_index == prev_theta_index ){
                  printf("ERROR : at (phi,theta) = (%.4f,%.4f) [deg] same indexes (%d,%d) as previously (%d,%d) at line = %d\n",phi_deg,theta_deg,phi_index,theta_index,prev_phi_index,prev_theta_index,prev_line);                  
                  printf("\t\tcurrent line  = %s\n",pLine);
                  printf("\t\tprevious line = %s\n",szPrevLine.c_str()); 
               }
               
//               if( fabs(phi_deg-126.6)<0.10 && fabs(theta_deg-10.5)<0.10 ){
//                  printf("odo");
//               }

               prev_phi_index = phi_index;
               prev_theta_index = theta_index;
               prev_line = line_counter;
               szPrevLine = pLine;

               if( phi_index >= MAX_PHI_COUNT ){
                  printf("ERROR : phi_index=%d (= (phi/step) = %.8f/%.8f ) exceeds size of antenna pattern array = %d\n",phi_index,phi_deg,phi_step,MAX_PHI_COUNT);
                  printf("Line = |%s|\n",pLine);
                  printf("phi_deg = %.8f , item[1] = %s\n",phi_deg,items[1].c_str());
                  exit(-1);
               }
               if( theta_index >= MAX_THETA_COUNT ){
                  printf("ERROR : phi_index=%d (= (phi/step) = %.8f/%.8f ) exceeds size of antenna pattern array = %d\n",theta_index,theta_deg,theta_step,MAX_THETA_COUNT);
                  printf("Line = |%s|\n",pLine);
                  printf("theta_deg = %.8f , item[1] = %s\n",theta_deg,items[0].c_str());
                  exit(-1);
               }
               
               current_freq_rec.ant_pattern[theta_index][phi_index] = gain;


               if( gain > max_gain ){
                  max_gain = gain;
               }
               if( gain < min_gain ){
                  min_gain = gain;
               }
            
               if( CAntPatternParser::m_bEnableInterpolation ){
#ifndef _NO_ROOT_
                   if( !current_freq_rec.pRootGraph ){
                      current_freq_rec.pRootGraph = new TGraph2D();                   
                   }
                   current_freq_rec.pRootGraph->SetPoint( freq_count, phi_deg, theta_deg, gain );
#else
                   printf("ERROR : root interpolation of antenna pattern is disabled by pre-compile directive _NO_ROOT_ . Re-compile the program or disable root interpolation\n");
                   exit(-1);                   
#endif
               }
            
               if( gBGPrintfLevel > 1 ){
                  printf("\tFilled gain of (theta,phi) = (%.2f,%.2f) [deg] at [%d,%d] = %.2f ( %.2f dB )\n",theta_deg,phi_deg,phi_index,theta_index,gain,gain_db);fflush(stdout);
               }
               freq_count++;
               sum_gain += gain;
            }else{
               if( freq_count > 0 ){
                  // new comment found - waiting for new frequency began 
                  m_FreqPatterns.push_back( current_freq_rec );
//                  if( gBGPrintfLevel > 1 ){
                     printf("Freq %.2f [MHz] finished #thetas=%d , #phis = %d, average_gain = %.2f (min=%.2f, max=%.2f)\n",current_freq_rec.freq,current_freq_rec.m_ThetaSamples,current_freq_rec.m_PhiSamples,(sum_gain/freq_count),min_gain,max_gain);
                     printf("-------------------------------------------------------------------------------------------------------------------------\n");
//                  }
               }else{
                  if( gBGPrintfLevel > 1 ){
                     printf("Comment line skipped : |%s|\n",pLine);
                  }
               }
            }
         }
      }else{
         printf("WARNING : empty line skipped !\n");
      }
      
      if( (line_counter % 10000) == 0 ){ // || line_counter>=323960000 ){
         printf("Read %d lines\n",line_counter);fflush(stdout);
      }
      line_counter++;
   }

   if( header_lines <= 0 ){
      printf("ERROR : missing header with FEKO file information -> cannot continue , please verify FEKO file %s and try again\n",pattern_file);
      exit(-1);      
   }
   
   printf("Read %d lines in total\n",line_counter);fflush(stdout);

   int zero_counter=0;
   for(int t=0;t<current_freq_rec.m_ThetaSamples;t++){
      for(int p=0;p<current_freq_rec.m_PhiSamples;p++){
         if( current_freq_rec.ant_pattern[t][p] == DEFAULT_PATTERN_VALUE ){
            // int phi_index = floor(phi_deg/phi_step);
            // int theta_index = floor(theta_deg/theta_step);
            double phi_deg   = p*phi_step;
            double theta_deg = t*theta_step;
                                    
            printf("NOT FILLED PATTERN at theta_index=%d (~%.4f deg) and phi_index=%d (~%.4f deg) (pattern = %d)\n",t,theta_deg,p,phi_deg,DEFAULT_PATTERN_VALUE);
            zero_counter++;
         }
      }
   }
   if(zero_counter>0){
      printf("ERROR !!!!!\n");
      printf("Number of zeros detected = %d of of total %d cells\n",zero_counter,(current_freq_rec.m_ThetaSamples*current_freq_rec.m_PhiSamples));fflush(stdout);
//      sleep(120);
      printf("Program might not be ready to work with your current resultion of Theta,Phi in the antenna pattern file - talk to MS to solve it !!!\n");
      printf("Exiting now !!!\n");
      exit(-1);      
   }

//   if( m_FreqPatterns.size() <= 0 && freq_count > 0 ){
   if( freq_count > 0 ){ // 20181003 fix to add last frequency in FFE file too !!! otherwise the last frequency is skipped !!!
      printf("Only 1 frequency read and added\n");
      m_FreqPatterns.push_back( current_freq_rec );
   }
   printf("Read %d frequnecies from FEKO formatted beam pattern file %s\n",(int)m_FreqPatterns.size(),pattern_file);fflush(stdout);
   
   return m_FreqPatterns.size();
}

int CAntPatternParser::FindFreq( double freq_mhz )
{
   double prev_freq=0.00;

   for(int i=0;i<m_FreqPatterns.size();i++){
      CFreqPattern& freq_info = m_FreqPatterns[i];
      
      if( fabs(freq_mhz-freq_info.freq) < 0.000001 ){
         return i;
      }else{
         if( prev_freq<freq_mhz && freq_mhz<freq_info.freq ){
            int ret = (i-1);
            if( ret < 0){
               ret = 0;
            }
            return ret;               
         }         
      }
      
      prev_freq=freq_info.freq;
   }   
   
//   return m_FreqPatterns.size();
   return -1;
}

double CAntPatternParser::GetGain(double freq,double phi_deg,double theta_deg)
{
   int freq_idx = FindFreq(freq);
   
   return GetGain( freq, freq_idx , phi_deg, theta_deg );
}

double CAntPatternParser::GetGain(double freq,CFreqPattern& pat_f1,CFreqPattern& pat_f2,double phi_deg,double theta_deg)
{
   int theta_samples = pat_f1.m_ThetaSamples;
   double theta_step = (m_SimulMaxTheta/(theta_samples-1));
   int phi_samples   = pat_f1.m_PhiSamples;
   double phi_step = (m_SimulMaxPhi/(phi_samples-1));

//if( pat_f1.m_ThetaSamples != m_FreqPatterns[0].m_ThetaSamples || pat_f1.m_PhiSamples != m_FreqPatterns[0].m_PhiSamples ){
//   printf("ERROR !!!\n");
//   exit(-1);
//}

   int phi_idx = (int)(phi_deg/phi_step);
   int theta_idx = (int)(theta_deg/theta_step);
   
   if( phi_idx <0 || phi_idx>=MAX_PHI_COUNT || theta_idx<0 || theta_idx>=MAX_THETA_COUNT ){
      printf("ERROR in CAntPatternParser::GetGain angle index (%d or %d) out of range (0,%d) and (0,%d) of (%.2f,%.2f) [deg]\n",phi_idx,theta_idx,MAX_PHI_COUNT,MAX_THETA_COUNT,phi_deg,theta_deg);
      exit(-1);
   }

   double gain1 = pat_f1.ant_pattern[theta_idx][phi_idx];
   double gain2 = pat_f2.ant_pattern[theta_idx][phi_idx];

#ifndef _SUPERFAST_OPTIMIZED_
#ifndef _NO_ROOT_
   if( pat_f1.pRootGraph && pat_f2.pRootGraph ){
      gain1 = (pat_f1.pRootGraph)->Interpolate( phi_deg, theta_deg );
      gain2 = (pat_f2.pRootGraph)->Interpolate( phi_deg, theta_deg );
         
      if( gain1 < 0.000001 ){
         gain1 = pat_f1.ant_pattern[theta_idx][phi_idx];
      }
      if( gain2 < 0.000001 ){
         gain2 = pat_f2.ant_pattern[theta_idx][phi_idx];
      }
   }
#endif   
#endif
      
   double gain = interpolate(freq,pat_f1.freq,gain1,pat_f2.freq,gain2);
   return gain;   
}
      
double CAntPatternParser::GetGain(double freq,int freq_idx,double phi_deg,double theta_deg)
{
//   int freq_idx = FindFreq(freq);

   int theta_samples = m_FreqPatterns[0].m_ThetaSamples;
   double theta_step = (m_SimulMaxTheta/(theta_samples-1));
   int phi_samples   = m_FreqPatterns[0].m_PhiSamples;
   double phi_step = (m_SimulMaxPhi/(phi_samples-1));

   int phi_idx = (int)(phi_deg/phi_step);
   int theta_idx = (int)(theta_deg/theta_step);
         
   
   if( freq_idx>=0 && freq_idx<m_FreqPatterns.size() && m_FreqPatterns.size()>=2 ){
      // WARNING : the if is over complicated not to deal with the references ...
      //           I should change it into POINTERS :
      if( (freq_idx+1) < m_FreqPatterns.size() ){
          CFreqPattern& pat_f1 = m_FreqPatterns[freq_idx];
          CFreqPattern& pat_f2 = m_FreqPatterns[freq_idx+1];
      
      
          double gain1 = pat_f1.ant_pattern[theta_idx][phi_idx];
          double gain2 = pat_f2.ant_pattern[theta_idx][phi_idx];

//      double phi_check = phi_deg/phi_step - round(phi_deg/phi_step);
//      double theta_check = theta_deg/theta_step - round(theta_deg/theta_step);
      
//      if( pat_f1.pRootGraph && pat_f2.pRootGraph && theta_deg<179.999 && theta_deg>0.00 && phi_deg<359.999 && phi_deg>0 ){

//      if( pat_f1.pRootGraph && pat_f2.pRootGraph && (phi_check>0.2 || theta_check>0.2) ){ // but keep values from the simulation, only interpolate inbetween them 
#ifndef _NO_ROOT_
          if( pat_f1.pRootGraph && pat_f2.pRootGraph ){
             gain1 = (pat_f1.pRootGraph)->Interpolate( phi_deg, theta_deg );
             gain2 = (pat_f2.pRootGraph)->Interpolate( phi_deg, theta_deg );
         
             if( gain1 < 0.000001 ){
                gain1 = pat_f1.ant_pattern[theta_idx][phi_idx];
             }
             if( gain2 < 0.000001 ){
                gain2 = pat_f2.ant_pattern[theta_idx][phi_idx];
             }
          }
#endif
      
          double gain = interpolate(freq,pat_f1.freq,gain1,pat_f2.freq,gain2);
          return gain;
      }else{
          CFreqPattern& pat_f1 = m_FreqPatterns[freq_idx];


          double gain1 = pat_f1.ant_pattern[theta_idx][phi_idx];

/*#ifndef _NO_ROOT_
          if( pat_f1.pRootGraph && pat_f2.pRootGraph ){
             gain1 = (pat_f1.pRootGraph)->Interpolate( phi_deg, theta_deg );
             gain2 = (pat_f2.pRootGraph)->Interpolate( phi_deg, theta_deg );
         
             if( gain1 < 0.000001 ){
                gain1 = pat_f1.ant_pattern[theta_idx][phi_idx];
             }
             if( gain2 < 0.000001 ){
                gain2 = pat_f2.ant_pattern[theta_idx][phi_idx];
             }
          }
#endif*/     
          double gain = gain1;
          return gain;      
      } 
   }else{
      if( freq_idx < 0 ){
         return m_FreqPatterns[0].ant_pattern[theta_idx][phi_idx];
      }else{
         return m_FreqPatterns[m_FreqPatterns.size()-1].ant_pattern[theta_idx][phi_idx];
      }
   }

   printf("ERROR : no gain value found between %d frequency patterns !!!\n",(int)m_FreqPatterns.size());   
   return -10000000.00;
}

/*double CFreqPattern::GetGain( double phi_deg,double theta_deg)
{
   double theta_step = (CAntPatternParser::m_SimulMaxTheta/(m_ThetaSamples-1));
   double phi_step = (CAntPatternParser::m_SimulMaxPhi/(m_PhiSamples-1));

   int phi_idx = (int)(phi_deg/phi_step);
   int theta_idx = (int)(theta_deg/theta_step);
   
   if( phi_idx <0 || phi_idx>=MAX_PHI_COUNT || theta_idx<0 || theta_idx>=MAX_THETA_COUNT ){
      printf("ERROR in CFreqPattern::GetGain angle index out of range of (%.2f,%.2f) [deg]\n",phi_deg,theta_deg);
      exit(-1);
   }

   double gain = ant_pattern[theta_idx][phi_idx];
   return gain;
}*/

CFreqPattern* CAntPatternParser::GetFreqPattern( int freq_idx )
{
   if( freq_idx>=0 && freq_idx<m_FreqPatterns.size() ){
      return &(m_FreqPatterns[freq_idx]);
   }else{
      if( freq_idx < 0 ){
         return &(m_FreqPatterns[0]);
      }
      if( freq_idx >= m_FreqPatterns.size() ){
         return &(m_FreqPatterns[m_FreqPatterns.size()-1]);
      }
   }
   
   return NULL;
}

CAntPatMap::CAntPatMap()
: map(NULL),m_Count(0)
{
}

CAntPatMap::~CAntPatMap()
{
   if( map ){
      delete [] map;
   }
}

CAntPatMap::CAntPatMap( const CAntPatMap& right )
: map(NULL),m_Count(0)
{
   (*this) = right;
}

CAntPatMap& CAntPatMap::operator=( const CAntPatMap& right )
{
   if( map ){
      if( m_Count != right.m_Count ){
         delete map;
         map = new double[right.m_Count];
      }
   }else{
      map = new double[right.m_Count];
   }   
   
   m_Count = right.m_Count;
   m_szFreqFile = right.m_szFreqFile;
   memcpy(map,right.map,sizeof(double)*m_Count);   
   
   return (*this);
}

int CAntPattern::m_bCacheON=0;
double CAntPattern::m_RotFromAxis=0.00; // for Muresk : 20-30 deg rotation of South direction from the X-axis
                                         // for NS 90 (Caiguna data)
string CAntPattern::m_AntennaPatternFile;
double CAntPattern::gSuppressAntPatternAboveZenAngle=10000;
int CAntPattern::gAntennaPatternFormula=0;
int          CAntPattern::gMWAGridpoint=-1;
CValueVector CAntPattern::m_AntennaPatternIntegralsdOMEGA;
eAntennaPatternType CAntPattern::m_AntennaPatternType=eAntPattIsotropic;
// double (*CAntPattern::m_UserDefinedPatternFormula)( double freq_mhz, double phi_deg, double theta_deg )=NULL;

CAntPattern::CAntPattern()
: m_AntennaPattern(NULL),m_AntennaBeamWidth(DEFAULT_ANT_BEAM_WIDTH)
{

}

CAntPattern::~CAntPattern()
{
   if( m_AntennaPattern ){
      delete m_AntennaPattern;
   }
}


double CAntPattern::antenna_pattern_formula( double freq_mhz, double phi_deg, double theta_deg )
{
   //if( m_UserDefinedPatternFormula && m_AntennaPatternType == eAntPattUserDefined ){
      // user pattern function is defined :
      //return (*m_UserDefinedPatternFormula)( freq_mhz, phi_deg, theta_deg );
   //}
   if( m_AntennaPatternType == eMWA2016 ){
      double ret =  CalcMWABeam( phi_deg, theta_deg, freq_mhz*1e6, 'X', gMWAGridpoint );
      return ret;
   }
   
   if( m_AntennaPatternType == eAntPattUserDefined ){
      return UserDefinedPatternFormula( freq_mhz, phi_deg, theta_deg );
   }

   double ant_patt=(1.00/(4.00*M_PI)); // default isotropic 

   if( m_AntennaPatternType == eAntPattGaussian ){
      double theta0 = m_AntennaBeamWidth;
      ant_patt = 1000.00*exp( -(theta_deg*theta_deg)/(theta0*theta0) );

      if( gSuppressAntPatternAboveZenAngle > 0 ){
         if( theta_deg > gSuppressAntPatternAboveZenAngle ){
            ant_patt = 0.00;
         }
      }
//      printf("Gaussuan beam = 1000*exp( -(theta_deg/%.2f)^2 )\n",theta0);
   }
   
   if( m_AntennaPatternType == eAntPattDipolOverGndScreen ){
      ant_patt = antenna_pattern_dipol( freq_mhz, phi_deg, theta_deg );
   }
               
   if( m_AntennaPatternType == eAntPattIsotropic ){
      ant_patt = (1.00/(4.00*M_PI));
   }
                              
   return ant_patt;   
}                                 
  
// double antenna_pattern( double phi_rad, double theta_rad, double freq_mhz )
double CAntPattern::antenna_pattern_dipol( double freq_mhz, double phi_deg, double theta_deg )
{
  double phi_rad   = deg2rad(phi_deg);
  double theta_rad = deg2rad(theta_deg);
     
  /*  if( theta_rad < M_PI_2 ){ // only above the ground :
      out = cos(theta_rad)*cos(theta_rad);
   }*/
  double freq_hz = freq_mhz * 1000000.00;
  double lambda = (c/freq_hz);
  double k_wave_num = (2*M_PI)/lambda;
  double h = 0.52; // [m] over the ground screen - ZMIERZYC WIECZOREM !!!
                   
  // horizontal infinitesmal dipol over ground screen :
  // double alpha = deg2rad(20.00); // 20-30 deg rotation of South direction from the X-axis   
  // double phi_rad = alpha + (2*M_PI-azim_rad);
  // double theta_rad = zenith_dist_rad;
  //  double ce300e_gain = get_ce300e_gain(freq_mhz);
                           
  // MODELLING 2012-06-27 or ce300e_gain
  // with FREQ factor included 
  //  double f_theta_phi = (1.00-sqr(sin(theta_rad))*sqr(sin(phi_rad)))*sqr(sin(2*M_PI*h/lambda*cos(theta_rad)));
                            
  // with FREQ factor from table :
  //  double f_theta_phi = (1.00-sqr(sin(theta_rad))*sqr(sin(phi_rad)))*get_ce300e_gain(freq_mhz);  
                            
  // without FREQ factor, because it was included in the calibration procedure :
  double f_theta_phi=0.00;
  //  if( theta_rad < M_PI_2 ){ - TEST 2012-10-28
  //     f_theta_phi = (1.00-sqr(sin(theta_rad))*sqr(sin(phi_rad)));
  f_theta_phi = (1.00-sqr(sin(theta_rad))*sqr(sin(phi_rad)))*sqr(sin(2*M_PI*h/lambda*cos(theta_rad)));
  //  }
                                         
  return f_theta_phi; 
}


double CAntPattern::antenna_pattern( double freq_mhz, 
                        CFreqPattern* pPatternFreq , CFreqPattern* pPatternFreqNext,
                        double azim_deg, double zenith_dist_deg )
{   
  // horizontal infinitesmal dipol over ground screen :
  double alpha_deg = m_RotFromAxis;
  double phi_deg = alpha_deg + (360-azim_deg);
  if( phi_deg < 0 ){
     phi_deg = 360 + phi_deg;
  }
  
  double theta_deg = zenith_dist_deg;

/*   if( gIncludeIonosphere ){     
      if( freq_mhz > gPlasmaFrequency ){
         double refraction_index = 1.00;
         refraction_index = sqrt( 1.00 - (gPlasmaFrequency*gPlasmaFrequency)/(freq_mhz*freq_mhz) );
         double sin_theta_ref = (1.00/refraction_index) * sin(theta_deg*M_PI/180.00);
      }
   }*/
           
   if( m_AntennaPattern ){
//      double check_val = antenna_pattern_dipol( freq_mhz, phi_deg, theta_deg );

//      double phi_simul_deg = 270.00 + alpha_deg - azim_deg;

      // This is formula in convention of libnova S=0 deg, W=90deg, N=180 deg, E=270 deg ( on both hemispheres ) - see libnova file :
      // /opt//caastro/ext/dload/libnova/libnova-0.13.0/src/transform.c
      // OR 
      // http://fornax.phys.unm.edu/lwa/doc/lsl-0.6.x/astro.html 
      // as long as all calculations are in the same convention it is OK
      // PHI in the simulation counts from the arm of the antenna it has been confirmed with : awk '{if($1==45){print $0;}}' 40.00.txt
      // 
      // FEKO OUTPUT (according to Shatanu) vs AZIMUTH (in libnova - as above) :
      // EAST(az=270)  : phi_simul_deg =   0 deg 
      // NORTH(az=180) : phi_simul_deg =  90 deg 
      // WEST(az=90)   : phi_simul_deg = 180 deg 
      // SOUTH(az=0)   : phi_simul_deg = 270 deg 
      // vs Randall convention (azimuth from N:az=0 , E:az=90, S:az=180, E:az=270) it transforms to libnova by az' := az - 180 deg 
// OLD(before 2013-05-29) :
//      double phi_simul_deg = azim_deg + (90.00 - alpha_deg);
// NEW - after 2013-05-29 - after checking with Adrian and Shatnanu :
/*      double phi_simul_deg = 270 - azim_deg + alpha_deg;
      if( phi_simul_deg >= 360 ){
         phi_simul_deg = phi_simul_deg - 360;
      }
      if( phi_simul_deg < 0 ){
         phi_simul_deg = 360 + phi_simul_deg;
      }*/
      double phi_simul_deg = azim2phi( azim_deg, alpha_deg );
      double ret = m_AntennaPattern->GetGain(freq_mhz, *pPatternFreq, *pPatternFreqNext, phi_simul_deg,theta_deg);

//      if( fabs(freq_mhz-prev_freq_mhz)>0.1 ){
//         printf("TEST_INFO : %.2f %.2f\n",freq_mhz,ret);
//         prev_freq_mhz = freq_mhz;
//      }
      
//      if( gDebug > 2 ){
//         printf("DEBUG : Gain(%.2f,%.2f) = %.8f vs DIPOL = %.8f \n",phi_deg,theta_deg,ret,check_val);
//      }
      return ret;
   }
   return antenna_pattern_formula( freq_mhz, phi_deg, theta_deg );   
}   

double CAntPattern::calculate_antenna_integral( double freq_mhz )
{
   cValue* pIntegral = m_AntennaPatternIntegralsdOMEGA.find_value( freq_mhz );
   if( pIntegral ){
      if( gBGPrintfLevel>=2 ){printf("CAntPattern::calculate_antenna_integral : returning pre-calculated integral at %.2f MHz = %.8f\n",pIntegral->x,pIntegral->y);}
      return pIntegral->y;
   }

   double dt_phi_deg=0.1; // 0.1 deg 
   double dt_theta_deg=0.1; // 0.1 deg
   double dt_theta_rad = deg2rad( dt_theta_deg );
   double dt_phi_rad = deg2rad( dt_phi_deg );

   double sum=0.00;
   double theta_deg=0;
   while( theta_deg <= 90.00 ){
   
      double phi_deg=0;      
      while( phi_deg <= 360.00 ){
         double theta_rad = deg2rad( theta_deg );
//         double phi_rad = deg2rad( phi_deg );
      
         double dOMEGA = sin(theta_rad)*dt_theta_rad*dt_phi_rad;         
         double gain = 0.00;
                       
         
         if( CAntPattern::gAntennaPatternFormula > 0 ){
              gain = antenna_pattern_formula( freq_mhz, phi_deg, theta_deg );
         }else{
              gain = m_AntennaPattern->GetGain(freq_mhz, phi_deg, theta_deg );
         }                
         sum += gain*dOMEGA;
         
         phi_deg += dt_phi_deg;
      }
      
      theta_deg += dt_theta_deg;
   }   

   if( gBGPrintfLevel>=2 ){printf("CAntPattern::calculate_antenna_integral integral of antenna gain (freq=%.2f MHz) = %.8f\n",freq_mhz,sum);}
   m_AntennaPatternIntegralsdOMEGA.add_value(freq_mhz, sum );
   return sum;
}

double CAntPattern::antenna_gain_from_azh( double freq_mhz, double azim_deg, double zenith_dist_deg, int bNormalised ) {
  // horizontal infinitesmal dipol over ground screen :
  double alpha_deg = m_RotFromAxis;
  double phi_deg = alpha_deg + (360-azim_deg);
  if( phi_deg < 0 ){
     phi_deg = 360 + phi_deg;
  }
  
  double theta_deg = zenith_dist_deg;
  
  double integral_freq = calculate_antenna_integral( freq_mhz );
  printf("DEBUG antenna_gain_from_azh integral_freq(%.2f MHz) = %.4f\n",freq_mhz,integral_freq);

  if( CAntPattern::gAntennaPatternFormula > 0 ){
     double ret =  antenna_pattern_formula( freq_mhz, phi_deg, theta_deg );
     if( gBGPrintfLevel>=2 ){printf("Returning value based on formula %.8f / (%.8f/4PI)\n",ret,integral_freq);}
     ret = ret / (integral_freq/(4.00*M_PI));
     return ret;
  }

   if( m_AntennaPattern ){
      int freq_index = m_AntennaPattern->FindFreq(freq_mhz);  // MS 2014-05-07 - otherwise calibration script was crashing : m_AntennaPattern-> added 
      CFreqPattern* pPatternFreq = m_AntennaPattern->GetFreqPattern( freq_index ); // MS 2014-05-07 - otherwise calibration script was crashing : m_AntennaPattern-> added 
      CFreqPattern* pPatternFreqNext = m_AntennaPattern->GetFreqPattern( freq_index+1 ); // MS 2014-05-07 - otherwise calibration script was crashing - m_AntennaPattern-> added 
   
      double phi_simul_deg = azim2phi( azim_deg, alpha_deg );
      printf("DEBUG : azim = %.2f [deg] -> FEKO_phi = %.2f [deg]\n",azim_deg,phi_simul_deg);
      double ret = m_AntennaPattern->GetGain(freq_mhz, *pPatternFreq, *pPatternFreqNext, phi_simul_deg, theta_deg );

      if( bNormalised ){      
         // added 2016-12-29 :         
         if( gBGPrintfLevel>=2 ){printf("Returning value %.8f / (%.8f/4PI)\n",ret,integral_freq);}
         ret = ret / (integral_freq/(4.00*M_PI));
      }else{
         if( gBGPrintfLevel>=2 ){printf("Returning value %.8f (not normalised)\n",ret);}
      }

//      if( fabs(freq_mhz-prev_freq_mhz)>0.1 ){
//         printf("TEST_INFO : %.2f %.2f\n",freq_mhz,ret);
//         prev_freq_mhz = freq_mhz;
//      }
      
//      if( gDebug > 2 ){
//         printf("DEBUG : Gain(%.2f,%.2f) = %.8f vs DIPOL = %.8f \n",phi_deg,theta_deg,ret,check_val);
//      }
      return ret;
   }
   return antenna_pattern_formula( freq_mhz, phi_deg, theta_deg );   
}   




void CAntPattern::Init( vector<double>& ni_list, int n_phi_samples, int n_theta_samples, double param_freq )
{
   if( strlen(CAntPattern::m_AntennaPatternFile.c_str()) ){
      m_AntennaPattern = new CAntPatternParser();
      if( m_AntennaPattern->Read( CAntPattern::m_AntennaPatternFile.c_str(), n_phi_samples, n_theta_samples, param_freq ) <= 0 ){
         printf("ERROR : could not antenna pattern file %s -> cannot continue\n",CAntPattern::m_AntennaPatternFile.c_str());
         exit(-1);
      }
   }                                    

   if( m_bCacheON  ){
      // COMMENT : remapped to map in AZIM,ALT , origianlly was in PHI,THETA - meaning results might be a little different ...
      printf("Antenna pattern cache is ON -> caching patterns for specfied frequencies\n");
      if( m_AntennaPattern ){   
         m_FreqPatterns.clear();
   
         if( n_phi_samples > MAX_PHI_COUNT || n_theta_samples > MAX_THETA_COUNT ){
            printf("ERROR in code : n_phi_samples=%d exceeds array size = %d or n_theta_samples=%d exceeds array size = %d\n",n_phi_samples,MAX_PHI_COUNT,n_theta_samples,MAX_THETA_COUNT);
            exit(-1);
         }
   
         for(int i=0;i<ni_list.size();i++){         
            double freq_mhz = ni_list[i];
            printf("Caching pattern for %.2f MHz\n",freq_mhz);fflush(stdout);
            
            int freq_index = m_AntennaPattern->FindFreq(freq_mhz);
            CFreqPattern* pPatternFreq = m_AntennaPattern->GetFreqPattern( freq_index );
            CFreqPattern* pPatternFreqNext = m_AntennaPattern->GetFreqPattern( freq_index+1 );
            
            if( !pPatternFreq || !pPatternFreqNext ){
               printf("ERROR in code : could not find antenna pattern for frequency index %d or %d (freq_mhz=%.2f MHz)\n",freq_index,(freq_index+1),freq_mhz);
               exit(-1);
            }    
            
            if( ! ( pPatternFreq->freq<=freq_mhz && freq_mhz<=pPatternFreqNext->freq ) ){
               printf("WARNING : frequency %.2f MHz outside interpolation range %.2f - %.2f MHz (interpolation of antenna pattern probably does not make sens !)\n",freq_mhz,pPatternFreq->freq,pPatternFreqNext->freq);
            }                                   
            
            if( pPatternFreq->m_PhiSamples > n_phi_samples ){
               n_phi_samples = pPatternFreq->m_PhiSamples;
            }
            if( pPatternFreq->m_ThetaSamples > n_theta_samples ){
               n_theta_samples = pPatternFreq->m_ThetaSamples;
            }
         
            double phi_step = m_AntennaPattern->m_SimulMaxPhi / (n_phi_samples-1);
            double theta_step = m_AntennaPattern->m_SimulMaxTheta / (n_theta_samples-1);

            CFreqPattern freq_patt( freq_mhz, n_theta_samples, n_phi_samples );
            for(int i_theta=0;i_theta<n_theta_samples;i_theta++){
               double theta_deg = theta_step * i_theta;
            
               for(int i_phi=0;i_phi<n_phi_samples;i_phi++){
                  double phi_deg = phi_step * i_phi;

                  // RE-MAPPING FROM (PHI,THETA) +/- alfaRot and 270 -> (AZIM,ALT) -> (little change of results) :                  
//                  freq_patt.ant_pattern[i_theta][i_phi] = antenna_pattern( freq_mhz, pPatternFreq, pPatternFreqNext, phi_deg, theta_deg ); 
                  // NO RE-MAPPING FROM (PHI,THETA) -> (AZIM,ALT) :
                  freq_patt.ant_pattern[i_theta][i_phi] = m_AntennaPattern->GetGain(freq_mhz, *pPatternFreq, *pPatternFreqNext, phi_deg, theta_deg );
               }
            }
         
            m_FreqPatterns.push_back( freq_patt );                    
            printf("MEMORY_REPORT : %.2f kB added\n",sizeof(freq_patt)/1000.00);
         }
      }
   }
}

void print_tab_diff( double* tab1, double* tab2, int count )
{
    if( !tab1 || !tab2 ){
       printf("WARNING : one of tables is not allocated tab1=0x%p tab2=0x%p\n",tab1,tab2);
       return;
    }

    double max_diff=0.00;
    int i_max=-1;
    for(int i=0;i<count;i++){
      double diff = fabs(tab1[i] - tab2[i]);
      if( diff > max_diff ){
         max_diff = diff;
         i_max=i;
      }
   }

   if( i_max >= 0 ){
      printf("Max difference between cached and current map = %e at i=%d %e != %e\n",max_diff,i_max,tab1[i_max],tab2[i_max]);
   }else{
      printf("Max difference between cached and current map = 0\n");
   }
}

double* CAntPattern::antenna_pattern( const char* infile, double freq_mhz, vector<cSkyInfo>& sky_intensities )
{
   // initialization of data structure for current antenna map :
   m_CurrentMap.m_szFreqFile = infile;
   if( !m_CurrentMap.map || m_CurrentMap.m_Count!=sky_intensities.size() ){
      if( m_CurrentMap.map ){
         delete m_CurrentMap.map;
      }
        
      m_CurrentMap.m_Count = sky_intensities.size();
      m_CurrentMap.map = new double[m_CurrentMap.m_Count];
   }

   if( m_bCacheON ){
      int freq_idx = FindFreq( freq_mhz );
      if( freq_idx >= 0 ){
         CFreqPattern* pAntPatFreq = GetFreqPattern( freq_idx );
         
         if( pAntPatFreq ){
            int size = sky_intensities.size();
            int k=0;
            for(vector<cSkyInfo>::iterator it=sky_intensities.begin();it!=sky_intensities.end();it++,k++){
               // REMAPPING : double phi_simul_deg = it->azim; when re-mapping (PHI,THETA) -> (AZIM,ALT)
               double phi_simul_deg = azim2phi( it->azim , m_RotFromAxis );
               
//               if( fabs(phi_simul_deg-123.7)<1.00 && fabs((90.00-it->alt)-11.40)<1.00 ){
//                  printf("odo");
//               }
               
               m_CurrentMap.map[k] = pAntPatFreq->GetGain( phi_simul_deg, (90.00-it->alt), CAntPatternParser::m_SimulMaxPhi, CAntPatternParser::m_SimulMaxTheta );
            }
           
           return m_CurrentMap.map;
         }else{
            printf("ERROR : antenna gain for frequency %.2f not found in cache\n",freq_mhz);
         }
      }
   }
   
   // cache not enabled :
   int freq_index = -1;
   CFreqPattern* pPatternFreq = NULL;
   CFreqPattern* pPatternFreqNext = NULL;
       
   if( m_AntennaPattern ){
       freq_index = m_AntennaPattern->FindFreq(freq_mhz);
       pPatternFreq = m_AntennaPattern->GetFreqPattern( freq_index );
       pPatternFreqNext = m_AntennaPattern->GetFreqPattern( freq_index+1 );           

       printf("integrate_galactic_spectrum : Freq index for %.2f MHz = %d -> pointers are 0x%p and 0x%p\n",freq_mhz,freq_index,pPatternFreq,pPatternFreqNext);
            
       if( !pPatternFreq || !pPatternFreqNext ){
          printf("ERROR in code : could not find antenna pattern for frequency index %d or %d (freq_mhz=%.2f MHz)\n",freq_index,(freq_index+1),freq_mhz);
          exit(-1);
       }
   }
       
   int k=0;
   for(vector<cSkyInfo>::iterator it=sky_intensities.begin();it!=sky_intensities.end();it++,k++){
       // strata czasu bo azim i alt sie zmieniaja - nie mozna tego tak zrobic !
       m_CurrentMap.map[k] = antenna_pattern( freq_mhz, pPatternFreq, pPatternFreqNext, it->azim, (90.00-it->alt) );
   }

   return m_CurrentMap.map;
}


double CFreqPattern::GetGain( double phi_deg, double theta_deg, double simul_max_phi, double simul_max_theta )
   {
      double theta_step = (simul_max_theta/(m_ThetaSamples-1));
      double phi_step = (simul_max_phi/(m_PhiSamples-1));
         
      int phi_idx = (int)(phi_deg/phi_step);
      int theta_idx = (int)(theta_deg/theta_step);

      if( theta_deg > CAntPatternParser::m_SimulMaxTheta && CAntPatternParser::m_SimulMaxTheta<180 ){
         // if simulation only up to 90 deg -> gain below horizon is 0.00 , but if we are requesting theta > 180 -> error in code !!!
         return 0.00;
      }else{
         if( phi_idx <0 || phi_idx>=MAX_PHI_COUNT || theta_idx<0 || theta_idx>=MAX_THETA_COUNT ){
            printf("ERROR in CFreqPattern::GetGain angle index out of range of (%.2f,%.2f) [deg]\n",phi_deg,theta_deg);
            printf("phi_idx   = %d vs. allowed range (0-%d)\n",phi_idx,MAX_PHI_COUNT);
            printf("theta_idx = %d vs. allowed range (0-%d)\n",theta_idx,MAX_THETA_COUNT);
            printf("ERROR in code - there should not be request for theta and phi outside range [0,180] and [0,360] deg\n");
            exit(-1);
         }     
      }
                                      
      double gain = ant_pattern[theta_idx][phi_idx];
      return gain;                                             
   }
