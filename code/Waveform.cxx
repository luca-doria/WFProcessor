#include "Waveform.h"
#include <iostream>
#include <string>
#include <cstring>
#include "TMath.h"
#include "fft2f.c"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH2.h"

Waveform::Waveform(int aNBin, double aMin, double aMax,
		   const char* aName, int aRebin):
  TH1D("HWF","",aNBin, aMin, aMax),
  mWriteId(0),
  mLastProcessingStage(0),
  mTimeDomain(1),
  mHBaseline(0),
  mBaselineMeanMin(-1e6),
  mBaselineMeanMax(1e6),
  mBaselineRMSMin(0.),
  mBaselineRMSMax(1e6),
  mFBaseline(0),
  mPed(0.),
  mStat(0),
  mRebin(aRebin),
  mVerMin(0),mVerScale(1),
  mVernierProfile(0),
  mHAverage(0),
  mNAverage(0),
  mHAverage2D(0),
  mPolarity(1),
  mHPowerSpectrum(0),
  mWFFilter(0),
  mFLTC6400(0)
{
  if(aName){
    SetName(aName);
  }
  else{
    int ti=0;
    char tName[50];
    sprintf(tName,"HWF%i",ti);
    while(gDirectory->Get(tName)){
      ti++;
      sprintf(tName,"HWF%i",ti);
    }
    SetName(tName);
  }
  //if(quietMode == false) std::cout << "Instantiating waveform " << GetName() << std::endl;
  initFourier();
}

Waveform::Waveform(const TH1D& aHisto, int aRebin, const char* aName):
  TH1D(aHisto),
  mWriteId(0),
  mLastProcessingStage(0),
  mTimeDomain(1),
  mHBaseline(0),
  mBaselineMeanMin(-1e6),
  mBaselineMeanMax(1e6),
  mBaselineRMSMin(0.),
  mBaselineRMSMax(1e6),
  mFBaseline(0),
  mPed(0.),
  mStat(0),
  mRebin(aRebin),
  mVerMin(0),mVerScale(1),
  mVernierProfile(0),
  mHAverage(0),
  mNAverage(0),
  mHAverage2D(0),
  mPolarity(1),
  mHPowerSpectrum(0),
  mWFFilter(0),
  mFLTC6400(0)
{
  if(aName){
    SetName(aName);
  }
  else{
    char tName[50];
    sprintf(tName,"%sWF",aHisto.GetName());
    SetName(tName);
  }
  //if(quietMode == false) std::cout << "Instantiating waveform " << GetName() << std::endl;
  initFourier();
}


void Waveform::initFourier(){
  int tNBin = GetNbinsX(); 
  mNFourierBin = (int) pow(2,floor(log(tNBin)/log(2.)));
  mCos = cos(TMath::Pi()/mNFourierBin);
  mSin = sin(TMath::Pi()/mNFourierBin);
 
  char ttName[100];
  if(mHPowerSpectrum){
    mHPowerSpectrum->Reset("ICE");
    mHRealSpectrum->Reset("ICE");
    mHCompSpectrum->Reset("ICE");
    mHPhaseSpectrum->Reset("ICE");
  }
  else{
    sprintf(ttName,"%sPS",GetName());
    //if(quietMode == false) std::cout << "Init Fourier:  NBins = " << mNFourierBin << " | ptr " <<mHPowerSpectrum << " | Name " << ttName <<  std::endl;
    mHPowerSpectrum = new TH1D(ttName,"",mNFourierBin/2,0.,mNFourierBin/2.);
    sprintf(ttName,"%sReal",GetName());
    mHRealSpectrum = new TH1D(ttName,"",mNFourierBin/2,0.,mNFourierBin/2.);
    sprintf(ttName,"%sComp",GetName());
    mHCompSpectrum = new TH1D(ttName,"",mNFourierBin/2,0.,mNFourierBin/2.);
    sprintf(ttName,"%sAmp",GetName());
    mHAmpSpectrum = new TH1D(ttName,"",mNFourierBin/2,0.,mNFourierBin/2.);
    sprintf(ttName,"%sPhase",GetName());
    mHPhaseSpectrum = new TH1D(ttName,"",mNFourierBin/2,0.,mNFourierBin/2.);
  }
}

void Waveform::setFilter(TF1& aFunction){
  if(mWFFilter) delete mWFFilter;
  char tName[50];
  sprintf(tName,"%sFilter",GetName());
  mWFFilter = new Waveform(GetNbinsX(),
			   GetXaxis()->GetXmin(),
			   GetXaxis()->GetXmax(),
			   tName);
  for(int ti=1; ti<=GetNbinsX(); ti++){
    mWFFilter->SetBinContent(ti,aFunction.Eval(mWFFilter->GetBinCenter(ti)));
  }   
}


Waveform::~Waveform(){
  if(mHBaseline) {
    mHBaseline->Delete();
    mHBaseline=0;
  }
  if(mFBaseline) {
    mFBaseline->Delete();
    mFBaseline=0;
  }
  mHPowerSpectrum->Delete();
  mHPowerSpectrum=0;
  mHRealSpectrum->Delete();
  mHCompSpectrum->Delete();
  mHPhaseSpectrum->Delete();
  mHAmpSpectrum->Delete();
  if(mWFFilter) delete mWFFilter;
}

void Waveform::setFixError(double aVal){
  for(int ti=1; ti<=GetNbinsX(); ti++){
    SetBinError(ti,aVal);
  }
}

void Waveform::setQuietMode(bool mode){
  quietMode = mode;
}

void Waveform::update(float* aWF){
  mStat=0;
  mTimeDomain=1;
  for(int ti=0; ti<GetNbinsX(); ti++){
    fArray[ti+1] = aWF[ti]*mPolarity;
    fArray[ti+1] -= mPed;;
  }
  if(mWFFilter) convoluteWith(mWFFilter);
}

void Waveform::update(int* aWF, int aVer){
  mStat=0;
  mTimeDomain=1;
  double tVal;

  for(int ti=0; ti<GetNbinsX(); ti++){
    tVal=0;
    for(int tj=0; tj<mRebin; tj++){
      tVal += (double) aWF[ti*mRebin+tj];
      tVal -= mPed;
    }
    fArray[ti+1] = tVal*mPolarity;
  }
  if(aVer) {
    double tTimeShift = (aVer-mVerMin)*mVerScale;
    if(mVernierProfile){
      tTimeShift += mVernierProfile->GetBinContent(aVer-mVernierProfileFirstBin+1);
    }
    timeShift(tTimeShift);
  }
  if(mWFFilter) convoluteWith(mWFFilter);
  //}
}

void Waveform::rebin(int aRebin){
  Rebin(aRebin);
  initFourier();
}

void Waveform::projectY(TH1D& aH){
  for(int ti=0; ti<GetNbinsX(); ti++){
    aH.Fill(fArray[ti+1]);
  }
}

void Waveform::initBaseline(int aNBin, double aMin, double aMax, int aSignalSign){
  if(mHBaseline) mHBaseline->Delete();
  char tBaselineName[50];
  sprintf(tBaselineName,"HBaseline%s",GetName());
  mHBaseline = new TH1D(tBaselineName,tBaselineName,aNBin,aMin,aMax);
  sprintf(tBaselineName,"FBaseline%s",GetName());
  mFBaseline = new TF1(tBaselineName,"gaus");
  mFBaseline->SetParameter(2,1.);
  mMin=aMin;
  mMax=aMax;
  mSignalSign=aSignalSign;
}
void Waveform::updateBaseline(int aAv){
  if(!mHBaseline){    
    mMin=1e6;
    mMax=-1e6;
    for(int ti=1; ti<=GetNbinsX(); ti++){
      if(mMin>fArray[ti]) mMin = fArray[ti];
      if(mMax<fArray[ti]) mMax = fArray[ti];
    }
    //double tLowerLim = floor(6.*mMin-5.*mMax)-0.5;
    //double tUpperLim = floor(6.*mMax-5.*mMin)+0.5;
    initBaseline((int) (floor(mMax)+1.-floor(mMin)),
		 floor(mMin),
		 floor(mMax),1);

//(int)(tUpperLim-tLowerLim),tLowerLim,tUpperLim,1);    
  }
  else{
    if(!aAv) mHBaseline->Reset("ICE");
  }
  double tSumX;
  for(int ti=1; ti<=GetNbinsX(); ti++){
    mHBaseline->Fill(GetBinContent(ti));
  }
  


}
void Waveform::fitBaseline(int aFit){
  int ti=1;
  double tVal=mHBaseline->GetBinContent(ti);
  double tEntries = mHBaseline->Integral();
  while(tVal/tEntries<0.18) {
    ti++;
    tVal+=mHBaseline->GetBinContent(ti);
  }
  mBaselineRMS =mHBaseline->GetBinLowEdge(ti);
  while(tVal/tEntries<0.5) {
    ti++;
    tVal+=mHBaseline->GetBinContent(ti);
  }
  mBaselineMean = mHBaseline->GetBinLowEdge(ti);
  mBaselineRMS = mBaselineMean-mBaselineRMS;
  double tMinBaselineRMS = 0.2887;//1./sqrt(12.)
  if(mBaselineRMS<tMinBaselineRMS) mBaselineRMS = tMinBaselineRMS;
  //cout << mBaselineMean << " " << mBaselineRMS << " " << aFit << endl;
  if(aFit){
    mFBaseline->SetRange(mBaselineMean-5.*mBaselineRMS,
			 mBaselineMean+3.*mBaselineRMS);
    mFBaseline->SetParameter(1,mBaselineMean);
    mFBaseline->SetParameter(2,mBaselineRMS);
    mHBaseline->Fit(mFBaseline,"QR0");
    mBaselineMean = mFBaseline->GetParameter(1);
    mBaselineRMS = mFBaseline->GetParameter(2);
    if(mBaselineRMS<tMinBaselineRMS) mBaselineRMS = tMinBaselineRMS;
  }
  //cout << mBaselineMean << " " << mBaselineRMS << endl;
}

void Waveform::subtractBaseline(int aFit, int aAv){
  //updateBaseline(aAv);
  //fitBaseline(aFit);
  simpleBaselineCalc();
  for(int ti=1; ti<=GetNbinsX(); ti++){    
    fArray[ti]-=mBaselineMean;
  }
  //SetMinimum(mMin-mBaselineMean);
  // SetMaximum(mMax-mBaselineMean);
}
double Waveform::getMeanBaseline(){
  return mBaselineMean;
}
double Waveform::getSigmaBaseline(){
  return mBaselineRMS;
}

void Waveform::simpleBaselineCalc(){
  if(!mHBaseline){
    char tBaselineName[50];
    sprintf(tBaselineName,"HBaseline%s",GetName());
    mHBaseline = new TH1D(tBaselineName,"",4000,-2000.,2000.);
  }
  else{
    mHBaseline->Reset("ICE");
  }
  double tMin=1e12;
  //double tMax=-1e12;
  double tAv=0.;
  //double tRMS;
  int tN =GetNbinsX(); 
  for(int ti=1; ti<=tN; ti++){    
    mHBaseline->Fill(fArray[ti]);
    tAv+=fArray[ti];
    //tRMS+=fArray[ti]*fArray[ti];
    if(tMin>fArray[ti]) tMin=fArray[ti];
    //if(tMax<fArray[ti]) tMax=fArray[ti];
  }
  tAv=mHBaseline->GetXaxis()->GetBinCenter(mHBaseline->GetMaximumBin());
  mBaselineMean=0.;
  mBaselineRMS=0.;
  tN=0;
  for(int ti=1; ti<=GetNbinsX(); ti++){ 
    //cout << fArray[ti] << endl;
    if(fArray[ti]<(tAv+(tAv-tMin))){
      mBaselineMean+=fArray[ti];
      mBaselineRMS+=(fArray[ti]*fArray[ti]);
      tN++;
    }
  }
  mBaselineMean/=tN;
  mBaselineRMS/=tN;
  mBaselineRMS-=(mBaselineMean*mBaselineMean);
  if(mBaselineRMS>1/12.){
    mBaselineRMS = sqrt(mBaselineRMS);
  }
  else{
    mBaselineRMS=1./sqrt(12.);
  }
}
  
TH1D* Waveform::getFrame(const char* aName){
  char tName[100];
  if(aName){
    std::strcpy(tName,aName);
  }
  else{
    sprintf(tName,"$sFrame",GetName());
  }
  TH1D* tH = (TH1D*) gROOT->FindObjectAny(tName);
  if(tH) tH->Delete();
  tH = new TH1D(tName,"",GetNbinsX(),GetXaxis()->GetXmin(),GetXaxis()->GetXmax());
  return tH;
}



int Waveform::convoluteWith(Waveform* aWaveform){
  aWaveform->switchToFrequencyDomain();
  switchToFrequencyDomain();
  double* tFFT1 = getArray();
  double* tFFT2 = aWaveform->getArray();
  tFFT1[0] *= tFFT2[0];
  tFFT1[1] *= tFFT2[1];
  double tR1;
  double tI1;
  for(int tj=1; tj<(mNFourierBin/2); tj++){
    tR1 = tFFT1[2*tj];
    tI1 = tFFT1[2*tj+1];
    tFFT1[2*tj]   = tR1*tFFT2[2*tj]  -tI1*tFFT2[2*tj+1];
    tFFT1[2*tj+1] = tR1*tFFT2[2*tj+1]+tI1*tFFT2[2*tj];
  }
  switchToTimeDomain();
  return 0;
}

int Waveform::deconvoluteWith(Waveform* aWaveform){
  aWaveform->switchToFrequencyDomain();
  switchToFrequencyDomain();
  double* tFFT1 = getArray();
  double* tFFT2 = aWaveform->getArray();
  tFFT1[0] /= tFFT2[0];
  tFFT1[1] /= tFFT2[1];
  double tR1;
  double tI1;
  double tNorm2;
  for(int tj=1; tj<(mNFourierBin/2); tj++){
    tR1 = tFFT1[2*tj];
    tI1 = tFFT1[2*tj+1];
    tNorm2 = tFFT2[2*tj]*tFFT2[2*tj]+tFFT2[2*tj+1]*tFFT2[2*tj+1];
    tFFT1[2*tj]   = (tR1*tFFT2[2*tj]+tI1*tFFT2[2*tj+1])/tNorm2;
    tFFT1[2*tj+1] = (tI1*tFFT2[2*tj]-tR1*tFFT2[2*tj+1])/tNorm2;
  }
  switchToTimeDomain();
  return 0;
}

void Waveform::fourierTransform(){
  // Use directly the array of the histo to save time
  // fArray[1] is the firs bin.
  //cout << GetName() << " " << mTimeDomain << "T" << endl;
  rdft(mNFourierBin, mCos, mSin, fArray+1);
  mTimeDomain=0;
} 

void Waveform::reverseFourierTransform(){
  // Use directly the array of the histo to save time
  // fArray[1] is the firs bin.
  //cout << GetName() << " " << mTimeDomain << "R" << endl;
  rdft(mNFourierBin, mCos, -1.*mSin, fArray+1);
  mTimeDomain=1;
  double* tFFT1 = getArray();
  for(int tj=0; tj<mNFourierBin; tj++){
    tFFT1[tj]/=mNFourierBin;
    tFFT1[tj]*=2.;
  }
}

double* Waveform::getArrayCopy() const{
  double* tArray = new double[GetNbinsX()];
  for(int ti=0; ti<GetNbinsX(); ti++){
    tArray[ti] = fArray[ti+1];
  }
  return tArray;
}

void Waveform::write(){
  if(mProcessingStage<=mLastProcessingStage) mWriteId++;
  mLastProcessingStage = mProcessingStage;
  char tName[100];
  sprintf(tName,"WF_E%iS%i",mWriteId,mProcessingStage);
  SetName(tName);
  Write();
}

double Waveform::complexIntegralSquared() const{
  double tSum=0.;
  for(int ti=1; ti<=mNFourierBin/2; ti++){
    tSum += fArray[2*ti-1]*fArray[2*ti-1]+fArray[2*ti]*fArray[2*ti];
  }
  return sqrt(tSum/GetBinWidth(1)/1e-9);
}


void Waveform::calcStat() const{
  mStat=1;
  mMin = 1e6;
  mMax = -1e6;
  mInt=0.;
  mIntP=0.;
  mIntN=0.;
  double tVal;
  for(int ti=0; ti<GetNbinsX(); ti++){
    tVal = fArray[ti+1];
    mInt+= tVal;
    if(tVal>0.){
      mIntP+=tVal;
    }
    else{
      mIntN-=tVal;
    }
    if(mMin>tVal) {
      mMin = tVal; 
      mMinPos = ti;
    }
    if(mMax<tVal){
      mMax = tVal;
      mMaxPos = ti;
    }
  }
  mIntTrunc=0.;
  double tN=0.;
  double tMinTrunc=-1;  
  mWidthTrunc=-1.;
  for(int ti=0; ti<GetNbinsX(); ti++){
    tVal = fArray[ti+1];
    if(tVal/mMin>0.1) {
      mIntTrunc-=tVal;
      tN++;
      if(tMinTrunc<0.) tMinTrunc=ti;
    }
    else{
      if(tMinTrunc>0.) mWidthTrunc = ti-tMinTrunc;
    }
  }
  mIntTrunc/=tN;
    
  // smooth max/min over 5 bins
  //int tLow = mMinPos>2? mMinPos-2:1;
  //int tUp = mMinPos<(GetNbinsX()-2)? mMinPos+2 : GetNbinsX();
  //for(int ti= tLow; ti<tUp; ti++){
}

void Waveform::drawRealSpectrum(const char* opt){
  fillFrequencySpectra();
  mHRealSpectrum->SetLineColor(GetLineColor());
  mHRealSpectrum->DrawCopy(opt);
}

void Waveform::drawCompSpectrum(const char* opt){
  fillFrequencySpectra();
  mHCompSpectrum->SetLineColor(GetLineColor());
  mHCompSpectrum->DrawCopy(opt);
}

void Waveform::drawPhaseSpectrum(const char* opt){
  fillFrequencySpectra();
  mHPhaseSpectrum->SetLineColor(GetLineColor());
  mHPhaseSpectrum->DrawCopy(opt);
}

void Waveform::drawAmpSpectrum(const char* opt){
  fillFrequencySpectra();
  mHAmpSpectrum->SetLineColor(GetLineColor());
  mHAmpSpectrum->DrawCopy(opt);
}

void Waveform::drawPowerSpectrum(const char* opt){
  fillFrequencySpectra();
  mHPowerSpectrum->SetLineColor(GetLineColor());
  mHPowerSpectrum->DrawCopy(opt);
}

void Waveform::resetFrequencySpectra(){
  for(int ti=1; ti<=mNFourierBin/2.; ti++){
    mHPowerSpectrum->SetBinContent(ti,0);
    mHRealSpectrum->SetBinContent(ti,0);
    mHCompSpectrum->SetBinContent(ti,0);
    mHPhaseSpectrum->SetBinContent(ti,0);
  }
}

void Waveform::fillFrequencySpectra(double aWeight){
  if(mTimeDomain) fourierTransform();
  for(int ti=1; ti<=mNFourierBin/2.; ti++){
    mHPowerSpectrum->Fill(ti-1,sqrt(fArray[2*ti]*fArray[2*ti]+
				  fArray[2*ti-1]*fArray[2*ti-1])*aWeight);
    mHRealSpectrum->Fill(ti-1,fArray[2*ti-1]*aWeight);
    mHCompSpectrum->Fill(ti-1,fArray[2*ti]*aWeight);
    mHPhaseSpectrum->Fill(ti-1,atan(fArray[2*ti]/fArray[2*ti-1])*aWeight);
    mHAmpSpectrum->Fill(ti-1,sqrt(fArray[2*ti]*fArray[2*ti]+fArray[2*ti-1]*fArray[2*ti-1])*aWeight);
  }
  switchToTimeDomain();
}

void Waveform::writeFrequencySpectra(){
  mHPowerSpectrum->Write();
  mHRealSpectrum->Write();
  mHCompSpectrum->Write();
}


void Waveform::setFrequencySpectra(const char* aFileName){
  TFile tFIn(aFileName);
  //  tFIn.ls();
  char ttName[100];
  sprintf(ttName,"%sPS",GetName());	  
  TH1D* tHPS = (TH1D*) tFIn.Get(ttName);
  sprintf(ttName,"%sReal",GetName());
  TH1D* tHReal = (TH1D*) tFIn.Get(ttName);
  sprintf(ttName,"%sComp",GetName());
  TH1D* tHComp = (TH1D*) tFIn.Get(ttName);
  for(int ti=1; ti<=tHPS->GetNbinsX(); ti++){
    mHPowerSpectrum->SetBinContent(ti,tHPS->GetBinContent(ti));
    mHRealSpectrum->SetBinContent(ti,tHReal->GetBinContent(ti));
    mHCompSpectrum->SetBinContent(ti,tHComp->GetBinContent(ti));
  }
}


void Waveform::subtractFrequencySpectra(){
  if(mTimeDomain) fourierTransform();
  for(int ti=1; ti<=mNFourierBin/2.; ti++){
    fArray[2*ti-1] -= mHRealSpectrum->GetBinContent(ti);
    fArray[2*ti]   -= mHCompSpectrum->GetBinContent(ti);
    mHPowerSpectrum->Fill(ti,sqrt(fArray[2*ti]*fArray[2*ti]+
				  fArray[2*ti-1]*fArray[2*ti-1]));
  }
  reverseFourierTransform();
}


void Waveform::filter(int aMin, int aMax){
  if(mTimeDomain) fourierTransform();
  //fArray[2*982-1]=0.;
  //fArray[2*982]=0.;
  //for(int ti=950; ti<=1024; ti++){
  //fArray[2*ti-1] =0.;//2;
  //fArray[2*ti] =0.;//2;
  //}
  //for(int ti=aMax; ti<=mNFourierBin/2; ti++){
  //fArray[2*ti-1] =0.;
  //fArray[2*ti] = 0.;
  //}
  int nMasked=8;
  int iMasked[8] = {206,308,410,411,413, 615, 616, 820};
  //for(int ti=200; ti<300; ti++){
  for(int ti=0; ti<nMasked; ti++){
    fArray[2*iMasked[ti]-1] =0.;
    fArray[2*iMasked[ti]] = 0.;
  }

  reverseFourierTransform();
}

void Waveform::resetAverage(){
  if(!mHAverage){
    char tName[50];
    sprintf(tName,"%sAv",GetName());
    mHAverage = (TH1D*) gROOT->FindObjectAny(tName);
    if(mHAverage) mHAverage->Delete();
    mHAverage = new TH1D(*((TH1D*) this));
    mHAverage->SetName(tName);
  }
  mHAverage->Reset();
  mNAverage=0;
  if(!mHAverage2D){
    char tName[50];
    sprintf(tName,"%sAv2D",GetName());
    mHAverage2D = (TH2D*) gROOT->FindObjectAny(tName);
    if(mHAverage2D) mHAverage2D->Delete();
    mHAverage2D = new TH2D(tName,tName,
			   GetNbinsX(),GetXaxis()->GetXmin(),GetXaxis()->GetXmax(),
			   4096*4,-2048.,2048.);
  }
  mHAverage2D->Reset();
}

void Waveform::timeShift(double aShift){
  switchToTimeDomain();
  int tBinShift = (int) floor(aShift/GetBinWidth(1));
  double tDiff = aShift/GetBinWidth(1)-tBinShift;
  double* tArray = new double[GetNbinsX()];
  int tiBin;
  for(int ti=0; ti<GetNbinsX(); ti++){
    tiBin = (ti+1)+tBinShift;
    if(tiBin>0 && tiBin<GetNbinsX()){
      tArray[ti] = GetBinContent(tiBin+1)*tDiff+ 
	GetBinContent(tiBin)*(1-tDiff);
    }
    else{
      tArray[ti] = 0.;
    }
  }
  for(int ti=0; ti<GetNbinsX(); ti++){
    fArray[ti+1] = tArray[ti];
  }
  delete[] tArray;
}

void Waveform::fillAverage(double aWeight, double aBinOffset){
  if(aBinOffset) timeShift(aBinOffset);
  //int tiBinOffset = (int) floor(aBinOffset);
  //double tDiff = aBinOffset-tiBinOffset;
  double tVal;
  mNAverage++;
  for(int ti=1; ti<=GetNbinsX(); ti++){
    tVal = GetBinContent(ti);
    mHAverage->SetBinContent(ti,mHAverage->GetBinContent(ti)+tVal);
    mHAverage->SetBinError(ti,mHAverage->GetBinError(ti)+tVal*tVal);
    mHAverage2D->Fill(mHAverage->GetBinLowEdge(ti),tVal);
    //tVal = GetBinContent(ti+1)*aWeight*tDiff;
    //mHAverage->SetBinContent(ti+tiBinOffset, 
    //		     mHAverage->GetBinContent(ti+tiBinOffset)+tVal);
    //mHAverage->SetBinError(ti+tiBinOffset, 
    //		   mHAverage->GetBinError(ti+tiBinOffset)+tVal*tVal);
    //tVal = GetBinContent(ti)*aWeight*(1.-tDiff);
    //mHAverage->SetBinContent(ti+tiBinOffset, 
    //		     mHAverage->GetBinContent(ti+tiBinOffset)+tVal);
    //mHAverage->SetBinError(ti+tiBinOffset, 
    //		   mHAverage->GetBinError(ti+tiBinOffset)+tVal*tVal);
    //if(mNAverage==1 && ti>200 && ti<700) cout << ti << " " <<  tVal << " " << tiBinOffset  << endl;
  }
}

TH1D* Waveform::getAverage(){
  if(mNAverage){
    for(int ti=1; ti<=GetNbinsX(); ti++){
      mHAverage->SetBinContent(ti,mHAverage->GetBinContent(ti)/mNAverage);
      mHAverage->SetBinError(ti,sqrt((mHAverage->GetBinError(ti)/mNAverage-
				      mHAverage->GetBinContent(ti)*mHAverage->GetBinContent(ti))/mNAverage));
    }
    mNAverage=0;
  }
  return mHAverage;
}


double Waveform::findCFDTime(int aFirstBin, int aLastBin){
  double tLocBaseline=0.;
  int tiBin;
  for(tiBin=0; tiBin<5; tiBin++) tLocBaseline += GetBinContent(tiBin);
  tLocBaseline/=5.;
  double tMax=-1000.;
  int tiBinMax;
  for(tiBin=aFirstBin; tiBin<=aLastBin; tiBin++){
    if(tMax<GetBinContent(tiBin)) {
      tMax = GetBinContent(tiBin);
      tiBinMax = tiBin;
    }
  }
  tMax=0.;
  for(tiBin=tiBinMax-1; tiBin<tiBinMax+2; tiBin++) tMax+=GetBinContent(tiBin);
  tMax/=3.;
  //tMax-=tLocBaseline;
  tiBin = tiBinMax;  
  while(GetBinContent(tiBin)/tMax>0.5) tiBin--;  
  double tDiff = (tMax/2.-GetBinContent(tiBin))/(GetBinContent(tiBin+1)-GetBinContent(tiBin));
  //std::cout << tLocBaseline << " " << tMax << " " << " " << tiBinMax << " " << tiBin << " " << tDiff << std::endl;
  return tDiff*(tiBin+1)+(1.-tDiff)*tiBin;

}

double Waveform::findTimeGauss(int aFirstBin, int aLastBin, int aShowFit){
  double tMax=-1000.;
  int tiBinMax;
  for(int tiBin=aFirstBin; tiBin<=aLastBin; tiBin++){
    if(tMax<GetBinContent(tiBin)) {
      tMax = GetBinContent(tiBin);
      tiBinMax = tiBin;
    }
  }
  TF1* myGauss = (TF1*) gROOT->FindObjectAny("myGauss");
  if(!myGauss) myGauss = new TF1("myGauss","gaus",GetBinLowEdge(aFirstBin-20),GetBinLowEdge(aLastBin+1));
  myGauss->SetParameters(tMax,GetBinLowEdge(tiBinMax),1.);
  myGauss->SetRange(GetBinLowEdge(aFirstBin-20),GetBinLowEdge(tiBinMax+2));
  if(aShowFit){
    this->Fit("myGauss","RQ");
    this->DrawCopy();
  }
  else{
    this->Fit("myGauss","RQ0");
  }
  return myGauss->GetParameter(1);
}


double FuncLTC6400(double* x, double* par){
  return par[1]+par[2]/2./par[3]*exp((par[4]*par[4]/2./par[3]-(x[0]-par[0]))/par[3])*
    TMath::Erfc(1./sqrt(2.)*(par[4]/par[3]-(x[0]-par[0])/par[4]));
}

double Waveform::findTimeLTC6400(double aMinRange, double aMaxRange, 
				 int aShowFit){
  //TF1* tFLTC6400 = (TF1*) gROOT->FindObjectAny("FLTC6400");
  if(!mFLTC6400) {
    mFLTC6400 = new TF1("FLTC6400",FuncLTC6400, aMinRange, aMaxRange,5);
    mFLTC6400->FixParameter(1,0.);
    mFLTC6400->FixParameter(3,22.32);
    mFLTC6400->FixParameter(4,1.25);
    mFLTC6400->SetParLimits(0,510.,600.);
  }
  mFLTC6400->SetParameter(0,537.);
  mFLTC6400->SetParameter(1,0.);
  mFLTC6400->SetParameter(2,850.);
  if(aShowFit){
    this->Fit("FLTC6400","RQ");
    this->DrawCopy();
  }
  else{
    this->Fit("FLTC6400","RQ0");
  }
  if(mFLTC6400->GetChisquare()<300.){
    return mFLTC6400->GetParameter(0);
  }
  else{
    return -1.;
  }
}


int Waveform::timeShiftLTC6400(double aMinRange, double aMaxRange,
				double aTimeMean){
  switchToTimeDomain();  
  double tTimeShift = findTimeLTC6400(aMinRange,aMaxRange,0);
  if(tTimeShift<0.) return 0;
  tTimeShift-=aTimeMean;


  int tBinShift = (int) floor(tTimeShift/GetBinWidth(1));
  double tDiff = tTimeShift/GetBinWidth(1)-tBinShift;
  double* tArray = new double[GetNbinsX()];
  int tiBin;
  for(int ti=0; ti<GetNbinsX(); ti++){
    tiBin = (ti+1)+tBinShift;
    if(tiBin>0 && tiBin<GetNbinsX()){
      
      tArray[ti] = 
	(GetBinContent(tiBin+1)-
	 mFLTC6400->Eval(GetBinLowEdge(tiBin+1)))*tDiff+ 
	(GetBinContent(tiBin)-
	 mFLTC6400->Eval(GetBinLowEdge(tiBin)))*(1-tDiff)+
	 mFLTC6400->Eval(GetBinLowEdge(tiBin)+tDiff);
    }
    else{
      tArray[ti] = 0.;
    }
  }
  for(int ti=0; ti<GetNbinsX(); ti++){
    fArray[ti+1] = tArray[ti];
  }
  delete[] tArray;

  return 1;
}


Waveform* Waveform::getSubset(int aFirstBin, int aLastBin, 
			      const char* aName){
  char tName[50];
  if(aName){    
    strcpy(tName,aName);
  }
  else{
    sprintf(tName,"%sSub%iTo%i",GetName(),aFirstBin,aLastBin);
  }  
  Waveform* tWaveform = (Waveform*) gROOT->FindObjectAny(tName);
  if(tWaveform) {
    //if(quietMode == false) std::cout << "Delete " << tName << std::endl;
    delete tWaveform;
  }
  tWaveform = new Waveform(aLastBin-aFirstBin+1,
			 GetBinLowEdge(aFirstBin),
			 GetBinLowEdge(aLastBin+1),
			 tName,
			 0);
  for(int tiBin=aFirstBin; tiBin<=aLastBin; tiBin++){
    tWaveform->SetBinContent(tiBin-aFirstBin+1,GetBinContent(tiBin));
  }
  return tWaveform;
}
