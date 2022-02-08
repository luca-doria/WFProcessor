#ifndef Waveform_h
#define Waveform_h
#include "TH1.h"
#include "TF1.h"
#include "TProfile.h"
#include "TH2.h"

class Waveform : public TH1D{
 public:
  Waveform(int aNBin,  double aMin, double aMax,
	   const char* aName=0, int aRebin=1); // always created in time domain
  Waveform(const TH1D& aWF, int aRebin=1, const char* aName=0);
  virtual ~Waveform();
  void setNegPolarity(){mPolarity=-1;}
  void setFilter(TF1& aFunction);
  
  Waveform* getSubset(int aFirstBin, int aLastBin, const char* aName=0);

  void setFixError(double aVal);
  void setVernierRange(int aMin, int aMax, double aDt){
    mVerMin=aMin; mVerScale = aDt/(aMax-aMin);}
  void setVernierRange(int aVerMean, double aScale){
    mVerMin=aVerMean; mVerScale = aScale;}
  void setSecondOrderVernier(TProfile* aProfile){mVernierProfile = aProfile;mVernierProfileFirstBin= (int) aProfile->GetBinLowEdge(1);}
  void setQuietMode(bool mode= true);

  // return an empty histogram with the same axis
  TH1D* getFrame(const char* aName=0);

  // Fill in the data
  void update(float* aWF);
  void update(int* aWF, 
	      int aVer=0);
  void clear();
  void timeShift(double aShift);
  void rebin(int aRebin);

  double findCFDTime(int aFirstBin, int aLastBin);
  double findTimeGauss(int aFirstBin, int aLastBin, int aShowFit=0);
  double findTimeLTC6400(double aMinRange, double aMaxRange, int aShowFit=0);
  int timeShiftLTC6400(double aMinRange, double aMaxRange, double aTimeMean);

  // --- Baseline/pedestal handling functions
  void initBaseline(int aNBin, double aMin, double aMax, int aSignalSign);
  void subtractBaseline(int aFit=0, int aAv=0); // try to measure baseline event by event
  void updateBaseline(int aAv=0);
  void fitBaseline(int aFit=0);
  void simpleBaselineCalc();
  TH1D* getBaseline() {return mHBaseline;}
  double getMeanBaseline();
  double getSigmaBaseline();
  int validBaseline(){return (mBaselineMean>mBaselineMeanMin && 
			      mBaselineMean<mBaselineMeanMax &&
			      mBaselineRMS>mBaselineRMSMin &&
			      mBaselineRMS<mBaselineRMSMax);}
  void setBaselineLimits(double aMeanMin, double aMeanMax,
			 double aRMSMin, double aRMSMax){
    mBaselineMeanMin=aMeanMin; mBaselineMeanMax=aMeanMax; 
    mBaselineRMSMin=aRMSMin; mBaselineRMSMax=aRMSMax;
  }
  // Function below used for pedestal
  void projectY(TH1D& aH);
  void setPedestal(double aPed){mPed=aPed;}

  // --- Average
  void resetAverage();
  void fillAverage(double aWeight=1., double aBinOffset=0.);
  TH1D* getAverage();
  TH2D* getAverage2D() {return mHAverage2D;}

  // FFT 
  int getNFourierBin() const {return mNFourierBin;}
  void drawPowerSpectrum(const char* opt="");
  void drawRealSpectrum(const char* opt="");
  void drawCompSpectrum(const char* opt="");
  void drawPhaseSpectrum(const char* opt="");
  void drawAmpSpectrum(const char* opt="");
  void fillFrequencySpectra(double aWeight=1.);
  void writeFrequencySpectra();
  void setFrequencySpectra(const char* aFileName);
  void subtractFrequencySpectra();
  void resetFrequencySpectra();
  //Waveform* getAverageFromFrequencySpectra();
  
  void switchToTimeDomain();
  void switchToFrequencyDomain();
  int convoluteWith(Waveform* aWF);
  int deconvoluteWith(Waveform* aWF);
  double* getArrayCopy() const; 
  double* getArray(); // use to speed up calculations
  void filter(int aMin=0, int aMax=0);
  double complexIntegralSquared() const;
  // These functions are not used for lxe
  void write(); // use write only with processing stage
  void setProcessingStage(int aStage); // not used for lxe

  // statistics
  double getIntegral() const;
  double getIntegralTrunc() const;
  double getIntegralPos() const;
  double getIntegralNeg() const;
  double getWidthTrunc() const;
  double getMin() const;
  double getMax() const;
  double getMinPos() const;
  double getMaxPos() const;

 private:
  //Waveform(const char* aName, TH1D* aHReal, TH2D* aHComp);


  void initFourier();
  void fourierTransform();
  void reverseFourierTransform();
  int mWriteId;
  int mProcessingStage;
  int mLastProcessingStage;

  // Average
  TH1D* mHAverage;
  int mNAverage;
  TH2D* mHAverage2D;

  // for FFT
  int mTimeDomain;
  double mCos;
  double mSin;
  int mNFourierBin;
  TH1D* mHPowerSpectrum;
  TH1D* mHRealSpectrum;
  TH1D* mHCompSpectrum;
  TH1D* mHPhaseSpectrum;
  TH1D* mHAmpSpectrum;
  Waveform* mWFFilter;


  // Vernier correction for GHz digitizer
  int mVerMin;
  double mVerScale;
  TProfile* mVernierProfile;
  int mVernierProfileFirstBin;

  //
  TH1D* mHBaseline;
  TF1* mFBaseline;
  mutable double mMin;
  mutable double mMax;
  double mBaselineMean;
  double mBaselineRMS;
  double mBaselineMeanMin;
  double mBaselineMeanMax;
  double mBaselineRMSMin;
  double mBaselineRMSMax;
  int mSignalSign;
  double mPed;
  int mRebin;
  int mPolarity;
  
  // statistics
  void calcStat() const;
  mutable int mStat;
  mutable double mInt;
  mutable double mIntTrunc;
  mutable double mIntP;
  mutable double mIntN;
  mutable double mMinPos;
  mutable double mMaxPos;
  mutable double mWidthTrunc;

  TF1* mFLTC6400;

  bool quietMode;

};

inline void Waveform::clear(){
  Reset("ICE");
}
inline void Waveform::setProcessingStage(int aStage){
  mProcessingStage = aStage;
}
inline void Waveform::switchToTimeDomain(){
  if(!mTimeDomain) reverseFourierTransform();  
}
inline void Waveform::switchToFrequencyDomain(){
  if(mTimeDomain) fourierTransform();
}
inline double* Waveform::getArray(){
  return fArray+1;
}
inline double Waveform::getIntegralPos() const{
  if(!mStat) calcStat();
  return mIntP;
}
inline double Waveform::getIntegralNeg() const{
  if(!mStat) calcStat();
  return mIntN;
}
inline double Waveform::getIntegral() const{
  if(!mStat) calcStat();
  return mInt;
}
inline double Waveform::getIntegralTrunc() const{
  if(!mStat) calcStat();
  return mIntTrunc;
}
inline double Waveform::getWidthTrunc() const{
  if(!mStat) calcStat();
  return mWidthTrunc;
}
inline double Waveform::getMin() const{
  if(!mStat) calcStat();
  return mMin;
}
inline double Waveform::getMax() const{
  if(!mStat) calcStat();
  return mMax;
}
inline double Waveform::getMinPos() const{
  if(!mStat) calcStat();
  return mMinPos;
}
inline double Waveform::getMaxPos() const{
  if(!mStat) calcStat();
  return mMaxPos;
}


#endif
