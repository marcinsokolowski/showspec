#ifndef _BG_FITS_H__
#define _BG_FITS_H__

#include <fitsio.h>
#include <string>
#include <vector>
#include "bg_array.h"
#include "bg_globals.h"
#include "bg_total_power.h"

// enum eBgFitsAction {eBgAct_None=0,eBgAct_Log10,eBgAct_DivideByConst,eBgAct_Sqrt,eBgAct_AstroRootImage};
enum eCalcFitsAction_T  {eNone=0,eAdd,eSubtract,eMultiply,eDivide,eSubtractLines,eCompare,eGetStat,eDivideConst, eLog10File, eSqrtFile, eAstroRootImage, eDivideLines,eTimesConst, eNormalizeByMedian, eFindValue, eDumpForHisto, eAvgChannels, eSubtractSpectrum, eDivbySpectrum, eDBFile, eAddConst, eSubtractMedian, eFindZeroInt, eLin2DB };
enum eFindValueType_T   {eFindValueExact=0,eFindValueSmaller,eFindValueLarger};


// integration types :
enum eIntType { eIntTypeUndefined=0, eIntTypeANT=1, eIntTypeREF=2, eIntTypeNSON=3 };

// Criteria statistics :
// change to BIT MASK to be able to know everything !!!
enum eRejectionReasonType { eRejectionReason_NOT_REJECTED=0, eRejectionReason_TotalPower, eRejectionReason_ChannelPower, eRejectionReason_ORBCOMM_Power, eRejectionReason_TotalPowerFreq, eRejectionReason_TotalPowerLocalCut };

#define REJECTION_NOT_REJECTED  0
#define REJECTION_TOTAL_POWER   0x01
#define REJECTION_CHANNEL_POWER 0x02
#define REJECTION_ORBCOMM_POWER 0x04
#define REJECTION_TOTAL_POWER_FREQ 0x08
#define REJECTION_TOTAL_POWER_LOCALCUT 0x10

#define BIGHORNS_RFI_VALUE -1000

using namespace std;

// integration range :
struct cIntRange
{ 
   int start_int;
   int end_int;
   int inttype;
   string m_szName; // state name
   int    StateIdx; // index of given type of integration 

   cIntRange(int start=0,int end=0): start_int(start), end_int(end) {}
   static const char* GetIntTypeDesc(eIntType inttype);
   static eIntType ParseIntType(const char* szIntType);
   int IsInRange( int integr );
};

  //! Structure to store header records
  struct HeaderRecord {         
                        std::string Keyword;
                        std::string Value;
                        std::string Comment;
                        char keytype;
                          
                        HeaderRecord(){
                           keytype = '\0';
                        }
                        
                        HeaderRecord& operator=(const HeaderRecord& right){
                           Keyword = right.Keyword;
                           Value   = right.Value;
                           Comment = right.Comment;
                           keytype = right.keytype;
                           
                           return (*this);
                        }

                        };



class CBgFits
{
protected :
  fitsfile* m_fptr;
  int total_counter;
  int m_lines_counter;
  string m_FileName; 
  int m_SizeX;
  int m_SizeY;
  int bitpix;
  int image_type;
  
  // fits DATA :  
  float* data;
  //! Fits header   
  vector<HeaderRecord> _fitsHeaderRecords;
  vector<cIntRange> m_IntegrationRanges;
    
public :
  double inttime; // exposure time 
  time_t dtime_fs;
  int    dtime_fu;
  int    naccum;
  double start_freq;
  double stop_freq;
  double delta_freq;
  vector<string> m_AverList;

  // just for some data problems :
  static int gFitsUnixTimeError;
  static const int m_TypicalBighornsChannels;
  static const int m_TypicalBighornsYSize;

  CBgFits( const char* fits_file=NULL, int xSize=0, int ySize=0 );
  CBgFits( int xSize, int ySize );
  ~CBgFits();
  void Clean();
  
  // memory allocation :
  void Realloc( int sizeX, int sizeY, int bKeepOldData=TRUE );

  // writing fits files piece by piece :
  int Create( const char* fits_file );
  int DumpFitsLine( float* buffer, int size );
  int add_line( float* buffer, int size );
  int add_line( CBgArray& line );
  int WriteKeys();
  int WriteKeys( fitsfile* fptr );
  int Close();
  int get_lines_counter() { return m_lines_counter; }
  int reset_lines_counter(){ m_lines_counter = 0; }
  

  int ReadFits( const char* fits_file=NULL, int bAutoDetect=0 );  
  int WriteFits( const char* fits_file, int bUpdateSizeY=0 );

  // headers :
  void SetKeyValues();
  void SetKeyword(const char *keyword, const char* new_value, char keytype='C', const char* comment=NULL);  
  void SetKeyword(const char *keyword, int new_value );
  void SetKeywordFloat(const char *keyword, float new_value );
  HeaderRecord* GetKeyword(const char *keyword);
  static HeaderRecord* GetKeyword(vector<HeaderRecord>& keys_list, const char *keyword);
  std::vector<HeaderRecord>& GetKeys(){ return  _fitsHeaderRecords; }
  void SetKeys( std::vector<HeaderRecord>& keys ){ _fitsHeaderRecords = keys; }
  double GetIntTime(){ return inttime; }
  void SetIntTime( double _inttime ){ inttime = _inttime; }
  double GetUnixTime();
  double GetUnixTime(int y);
  int SetKeysWithoutStates( std::vector<HeaderRecord>& keys ); // does not remove old keys 
  void ClearKeys();
  void SetIntTimeKeyword( double _inttime );
  void PrepareBigHornsHeader( double ux_start, double _inttime, double freq_start, double delta_freq_mhz );

  // ch -> freq 
  double ch2freq(int ch);
  int    freq2ch(double freq);

  // Flagging : ANTENNA vs REFERENCE integrations 
  int ParseSkyIntegrations();
  int ParseStates( const char* szStatesList );
  vector<cIntRange>& GetIntRanges(){ return m_IntegrationRanges; }
  int GetRangesCount(){ return m_IntegrationRanges.size(); }
  cIntRange* GetRange(int y, int& out_range_idx);
  int IsAntenna(int y,int bDefaultYes=1);
  int IsReference(int y);
  eIntType GetIntType(int y);
  cIntRange* GetRange(int idx, eIntType inttype);
  double GetChannelWidth();
  int GetAntRanges( vector<cIntRange>& ranges );

  // 
  int GetXSize(){ return m_SizeX; }
  int GetYSize(){ return m_SizeY; }
  void SetYSize( int ySize ){ m_SizeY = ySize; }
  
  inline float* get_data(){ return data; }
  const char* GetFileName(){ return m_FileName.c_str(); }  

  CBgFits& operator+=(const CBgFits& right);    
  void Normalize(double norm_factor);
  float setXY( int x, int y, float value );
  int setY( int y, float value );
  float* set_line( int y, float* buffer );
  float* set_line( int y, CBgArray& line );
  float* set_line( int y, vector<cValue>& line );
  void set_ysize( int lines_counter=-1 );

  float valXY( int x, int y );
  inline float getXY( int x, int y ){ return valXY(x,y); }
  float value( int y, int x );
  float* get_line( int y );
  float* get_line( int y, float* buffer );
  float* get_line( int y, CBgArray& buffer );
  
  // char image :
  char valXY_char( int x, int y );
  float valXY_auto( int x, int y );
  
  // arithmetical operations :
  void Divide( CBgFits& right );
  void Multiply( CBgFits& right );
  void Subtract( CBgFits& right );
  void SubtractLines( int y1, int y0, const char* outfile="diff.txt");
  void DivideLines( int y1, int y0, const char* outfile="ratio.txt");
  int Compare( CBgFits& right, float min_diff=0.00001, int verb=0 );
  double GetStat( CBgArray& avg_spectrum, CBgArray& rms_spectrum, int start_int=0, int end_int=-1, const char* szState=NULL,
                  CBgArray* min_spectrum=NULL, CBgArray* max_spectrum=NULL,
                  double min_acceptable_value=-1e20, int* out_number_of_used_integrations=NULL,
                  CBgFits* rfi_flag_fits_file=NULL );
  void Divide( double value );
  int Recalc( eCalcFitsAction_T action, double value=0.00 );
  int GetMedianInt( CBgArray& median_int, CBgArray& rms_iqr_int );
  int GetMedianInt( vector<cIntRange>& int_ranges, CBgArray& median_int, CBgArray& rms_iqr_int );
  void Normalize( CBgArray& median_int );
  int FindValue(double value, double delta=1000, eFindValueType_T type=eFindValueExact);
  void SetValue(double value);
  void AvgChannels(int n_channels, CBgFits& outfits  );
  int SubtractSpectrum( CBgArray& spectrum );
  int DivideBySpectrum( vector<cValue>& spectrum );
  int ZeroWrongValues( double limit=10e20 );
  
  // operations on list of files 
  int CalcMedian( vector<string>& fits_list, CBgFits& out_rms, int bDoAverage=0 );

  // range operations :
  void dump_max_hold( int start_int, int end_int, const char* szOutFile, int bShowFreq=0 );
  void dump_min_hold( int start_int, int end_int, const char* szOutFile, int bShowFreq=0 );
  
  // physical operations :
  // calculate sum of values in given FITS-row (integration <=> Y) :
  double GetTotalPower( int integration, int start_channel=0, int end_channel=-1 );
  double GetMaxPower( int integration, double& max_freq, int start_channel=0, int end_channel=-1 );
  double GetTotalPowerFreq( int integration, double start_freq, double end_freq );
  
  // fit pol N
  double FitPoly( int y, double fit_min_freq, double fit_max_freq, double& A, double& B );
  
  // RFI checking etc :
  CBgFits* m_pRFIMask;
  CBgFits* GetRFIMask();  
  int FlagInt( int y );
  int FlagInt( int y, CBgFits* pRFIMask  );
  int FlagFile();
  static string gInAOFlaggerDir; // directory with ao-flagger masked files 
  int IsFlagged( int integration );
  int IsRFI_OK(  int integration, double& out_total_power, double& max_ch_power_dbm, double& max_freq, double& local_threshold, double& orbcomm_total_power, int& out_rejection_reason );
  int SaveAsByte( const char* outfile );
  
  // managing output files :
  CBgFits* AllocOutFits( const char* fname, int _y_size, int bAddStates=TRUE );
};

#endif
