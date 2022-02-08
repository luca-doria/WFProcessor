#ifndef WaveformProcessor_h
#define WaveformProcessor_h

class TH1F;
class TGraph;
class DataFile;
class TF1;
class Waveform;

struct Pulse{
  double mBaseline;
  double mAbsAmp;
  double mAmp;
  double mTime; 
  double mQ;
  double mWidth;
  double mSPTemplateChi2; 
  double mFitLowLimit;
  double mFitHighLimit;
  double mFitBaseline;
  double mFitTime;
  double mFitAmp;
  double mFitRiseTime;
  double mFitFallTime;
  double mFitTime2Frac;
  double mFitFallTime2;
  double mFitChi2;
  double mFitNDF;
  double mFitQ;
  double mRefitChi2;
  double mRefitNDF;
  int mFirstPulseInGroup;
};

class WaveformProcessor{
 public:
  WaveformProcessor(const char* aSetupFileName, int aRun,  int numPulsesFit=6, int aChannel=1);
  ~WaveformProcessor();
  static WaveformProcessor* instanceForRootMacro(const char* aSetupFileName, int aRun,  int numPulsesFit=6, int aChannel=1);
  
  // ---
  static int secDiffForRootMacros(const char* aFirstFileName,
				  const char* aFileName);

  // ---
  char* getFileName(){return mFileName;}
  int getWaveformCount();
  double getTriggerTime();

  // --- Access  to waveform
  bool readNextWaveform();
  void readWaveform(int aIndex);
  void readWaveform(const char* aHName);
  int getCurrentWaveformIndex();
  Waveform* getCurrentWaveform(){return mWF;}
  
  // --- baseline
  int processBaseline();
  TH1F* getBaselineHisto(){return mHBaseline;}
  double getBaselineMu(){return mBmu;}
  double getBaselineRMS(){return mBRMS;}
  bool PulseFound(){return mPulseFound;}
  
  // --- pulse finding
  int findPulse();
  int getNPulse(){return mNPulse;}
  double getPulseTime(int aPulse){return mPulse[mPulseIndex[aPulse]].mTime;}
  double getPulseAbsAmplitude(int aPulse){return mPulse[mPulseIndex[aPulse]].mAbsAmp;}
  double getPulseAmplitude(int aPulse){return mPulse[mPulseIndex[aPulse]].mAmp;}
  double getPulseBaseline(int aPulse){return mPulse[mPulseIndex[aPulse]].mBaseline;}
  double getPulseCharge(int aPulse){return mPulse[mPulseIndex[aPulse]].mQ;}
  double getPulseWidth(int aPulse){return mPulse[mPulseIndex[aPulse]].mWidth;}
  TGraph* getPulseGraph();
  
  // --- fit Pulse
  void calcSinglePulseTemplateChi2();
  void setFitOption(const char* aOption);
  void fit(int aType);  

  //  
  double getpulsepolarity(){return mPolarity;}
  double getSPTemplateChi2(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mSPTemplateChi2;}
  double getFitBaseline(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mFitBaseline;}
  double getFitRiseTime(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mFitRiseTime;}
  double getFitFallTime(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mFitFallTime;}
  double getFitTime2Frac(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mFitTime2Frac;}
  double getFitFallTime2(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mFitFallTime2;}
  double getFitAmplitude(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mFitAmp;}
  double getFitTime(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mFitTime;}
  double getFitIntegral(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mFitQ;} //INTEGRAL
  double getChi2(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mFitChi2;}
  double getNDF(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mFitNDF;}
  double getChi2Refit(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mRefitChi2;}
  double getNDFRefit(int aPulse=0){return mPulse[mPulseIndex[aPulse]].mRefitNDF;}
  double get_blrms(){return mDiffAmpThresh;}
  TF1* getFitFunction();
  TF1* getFitFunction2();

  TH1F* getPulseHisto();
  void setPulseHisto(int nPulseBins, double preFactor){}
  void setQuietMode(bool mode= true);

  int getChannel(){return mChannel;}

 private:
  static WaveformProcessor* mInstanceForRootMacro;
  static int mRun;
  static int mChannelUnique;
  void init(const char* aSetupFileName, int aRun, int numPulsesFit, int aChannel=1);
  void clear();
  void reFitMulti(); // controlled by type
  
  // >>> parameters
  char mSetupFileName[200];	
  char mFileName[200];
  int mChannel;
  double mWFAmpBining;
  double mAmpMin;
  double mAmpMax;
  double mNoise1Bin;
  double mDiffAmpThresh;
  double mRiseTimeSig;
  double mFallTimeTau;
  double mTime2Frac;
  double mFallTime2Tau;
  int mPolarity;
  double mMinChi2ForRefit;
  double mMaxSPTemplateChi2;

  // >>> Raw data
  DataFile* mDataFile;
  Waveform* mWF;

  // >>> Tools
  TH1F* mHBaseline;
  double mBmu;
  double mBRMS;
  
  // >>> Pulse
  int mNPulse;
  int mNPulseMax;
  Pulse* mPulse;
  int* mPulseIndex;
  double* mPulseTime;//for sorting only
  bool mPulseFound;

  // >>> Fit 
  TF1* mFFit;
  TF1* mFExpGaus;
  TF1* mFExp2Gaus;
  TF1* mFExpGaus2;
  TF1* mFExpRise;
  char mFitOption[10];	
  int mNFitPulseMax;

  TH1F* pulseHisto; // follow convention
  int fitstatus; // follow convention
  bool quietMode;// follow convention
		
};

#endif
