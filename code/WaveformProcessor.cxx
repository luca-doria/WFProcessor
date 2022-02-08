/*******************************************************************************
*  WaveformProcessor.cxx 
* 
*  Main TRIUMF code for SiPM pulse fitting
*
*  Code for waveform analysis. Right now the code is optimized for Hamamatsu Devices
*  Some work must be done in order to make it suitable for FBK devices
*
*  Description:
*  Reads a wfm file and produces an ntuple with the distribution of the pulses
*  for waveforms in that wfm file.
*
*  Comments:
*
* 1)  I have changed the way in which the Pulsewidth is computed. The previous condition (line 403) 
*     required that the Pulse End Amplitude was 5 times the noise, If the noise level it is high 
*     Then the Pulse end amplitude was too small and the final condition fails. The only way is to
*     adjust this value in function of the noise level
*
* 2)  Changing the frequency of pulsefinding line 581 is useful when you need to find pulses! (For example for PDE) but 
*     change also the fit behaviour in fact it tends to add pulses also where they shouldn't be. For this reason the most 
*     safiest thing is to not change the code here and to use the new macro Threshold to cunt the pulses in a certain windoe
*     for the PDE. If you change and you don't come back before fitting, fitting will be much worste
*
*  3) In order to compile properly code and macro you need to use the make file in particular 
*     type: make PF (comple pulsefinding) and make lib(compile macro)
*
*  4) Hard coded selection of the fit function 
*
*     giacomo@triumf.ca


 History:

 - 26/02/2019 Line 1100: Added condition to add pulses accordingly to pulse polarity
              Line 732 : Extended Left Integration Limit for pulse fitting
  -27/02/2019 Line 959 : Added sorting for the .1 and .2 fitting. Before there was only the sorting for .101 fitting  


 Known Bugs

 - 27/02/2019 In the refitting: If the added pulse produce an increase in the ChiSquare and no pulses are added to the waveform, then the sorting not happen. 
              (Ex. you start with 2 pulses and you end with 2 pulses since the third you try to add produce an increase in ChiSquare, then the sorting of the pulses  
              don't happen).  
*
*******************************************************************************/

using namespace std;

#include <sys/stat.h>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "DataFile.h"
#include "LecroyFile.h"
#include "V1730File.h"
#include "TextFile.h"
#include "WaveformProcessor.h"
#include "Waveform.h"


#include "TFile.h"
#include "TF1.h"
#include "TROOT.h"
#include "TH1.h"
#include "TGraph.h"
#include "TMath.h"

#include <Math/MinimizerOptions.h>


//___________________________________________________________________
// --- Fit function
double FuncExpGausMulti(double* x, double*p){
  //std::cout<<"Try fitting FuncExpGausMulti..."<<std::endl;
  // p[0]: baseline
  // p[1]: gaussian sig
  // p[2]: exponential decay constant
  // p[3]: number of pulses
  // p[4+2*i]: pulse[i] amplitude
  // p[5+2*i]: pulse[i] time
  // convolution of an exponential and a gaussian, formula found in a
  // chromatography paper on exponentially modified gaussians.
  // http://www.springerlink.com/content/yx7554182g164612/
  //https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution	
	double val=p[0];
	for(int iPulse=0; iPulse<p[3]; iPulse++){
		double time = x[0]-p[5+2*iPulse];
		if(((p[1]*p[1]/p[2]/2.-time)/p[2])<700){ //otherwise exponential explodes
		  val+=p[4+2*iPulse]/2./p[2]*	
			  exp((p[1]*p[1]/p[2]/2.-time)/p[2])*
			  TMath::Erfc((p[1]*p[1]/p[2]-time)/sqrt(2)/p[1]);
		}
		//exp(1/2*p[1]*p[1]/p[2]/p[2]-time/p[2])*
		//(TMath::Erf(1/sqrt(2)*(p[5+2*iPulse]/p[1]+p[1]/p[2]))+
		// TMath::Erf(1/sqrt(2)*(time/p[1]-p[1]/p[2])));
	}
	return val;
}

double FuncExp2GausMulti(double* x, double*p){
  //std::cout<<"FuncExpGausMult..."<<std::endl;
  // p[0]: baseline
  // p[1]: gaussian sig
  // p[2]: first  exponential decay constant
  // p[3]: second exponential fraction
  // p[4]: second exponential decay constant
  // p[5]: number of pulses
  // p[6+2*i]: pulse[i] amplitude
  // p[7+2*i]: pulse[i] time

                                                                         
                                                                     

    
  double val=0.;
  for(int iPulse=0; iPulse<p[5]; iPulse++){
    double time = x[0]-p[7+2*iPulse];
    if(((p[1]*p[1]/p[2]/2.-time)/p[2])<700){ //otherwise exponential explodes
      val+=(p[6+2*iPulse]/2./p[2]*exp((p[1]*p[1]/p[2]/2.-time)/p[2])*
	    TMath::Erfc((p[1]*p[1]/p[2]-time)/sqrt(2)/p[1]))*(1-fabs(p[3]));
      if(p[3]>0.){
	val+=(p[3]* p[6+2*iPulse]/2./p[2]*exp((p[1]*p[1]/p[4]/2.-time)/p[4])*
	      TMath::Erfc((p[1]*p[1]/p[4]-time)/sqrt(2)/p[1]));
      }
      if(p[3]<0.){
	val+=(fabs(p[3])* p[6+2*iPulse]/2./p[2]*exp((p[4]*p[4]/p[2]/2.-time)/p[2])*
	      TMath::Erfc((p[4]*p[4]/p[2]-time)/sqrt(2)/p[4]));
      }
    }
  }

  

  val+=p[0];
  return val;
}

double FuncExpRiseMulti(double* x, double*p){
  //std::cout<<"FuncExpRiseMulti.."<<std::endl;
  // p[0]: baseline
  // p[1]: rise time
  // p[2]: fast fall time - rise time
  // p[3]: slow fraction
  // p[4]: slow fall time - rise time - fast fall time
  // p[5]: number of pulses
  // p[6+2*i]: pulse[i] charge
  // p[7+2*i]: pulse[i] time
  // subtract rise to avoid dividing by zero by setting limit >0. on time constants and ensure time constants remain in order
  double val=0.;
  for(int iPulse=0; iPulse<p[5]; iPulse++){
    double time = x[0]-p[7+2*iPulse];
    if(time>0.){
      val+=p[6+2*iPulse]*((1.-p[3])*(exp(-time/(p[2]+p[1]))-exp(-time/p[1]))/p[2]+p[3]*(exp(-time/(p[4]+p[1]+p[2]))-exp(-time/p[1]))/(p[4]+p[2]));
      //val+=p[6+2*iPulse]*((1.-p[3])*(exp(-time/p[2])-exp(-time/p[1]))/(p[2]-p[1])+p[3]*(exp(-time/p[4])-exp(-time/p[1]))/(p[4]-p[1]));
    }
  }
  val+=p[0];
  //cout << x[0] << " " << val <<endl;
  return val;
}


//___________________________________________________________________
// --- Constuctor
WaveformProcessor::WaveformProcessor(const char* aSetupFileName, int aRun,  int numPulsesFit,int aChannel){
  init(aSetupFileName,aRun,numPulsesFit,aChannel);
}
// ---
void WaveformProcessor::init(const char* aSetupFileName, int aRun,  int numPulsesFit, int aChannel){

  //  cout<<"sono in init"<<endl;

  mHBaseline=0;
  mFExpGaus=0;
  mFExpRise=0;
  mPulse=0;
  mPulseIndex=0;
  mPulseTime=0;
  mChannel=aChannel;
  int tNBuf;
  if(aChannel==-1){// Setup file contains channel #
    tNBuf=13;
  }
  else{
    tNBuf=12;
  }
  
  strcpy(mSetupFileName,aSetupFileName);
  // --- Open run info file
  ifstream parFile(aSetupFileName);
  char buf[200];
  for(int iBuf=0; iBuf<tNBuf; iBuf++) parFile >> buf; 
	//mWFAmpBining >> mAmpMin >> mAmpMax not currently used in this code. See WaveformProcessor::processBaseline()
  std::cout << buf << std::endl;
  parFile >> mRun; 
  while(!parFile.eof() && mRun!=aRun){
    for(int iBuf=0; iBuf<tNBuf-1; iBuf++) parFile >> buf; 
    parFile >> mRun; 
  }
  if(parFile.eof()){
    std::cerr << "Could not find run requested " << aRun << std::endl;
    exit(0);
  }
  if(mChannel==-1) parFile >> mChannel;
  parFile >> mFileName >>  mWFAmpBining >> mAmpMin >> mAmpMax 
	  >> mNoise1Bin >> mRiseTimeSig >> mFallTimeTau >>  mTime2Frac >>  mFallTime2Tau >> mPolarity >> mMinChi2ForRefit;

  //std::cout << mFallTime2Tau << std::endl;


  mDiffAmpThresh = 5*mNoise1Bin; // 5 times the noise baseline sigma should be sufficient to suppress bg efficient
  mMaxSPTemplateChi2 = 1e9;
  
  // --- Open input file
  //>>> check file extension
  int iChar=0;
  while(iChar<strlen(mFileName) && strncmp(mFileName+iChar,".",1)!=0) iChar++;
  //std::cout << mFileName << " " << mFileName+iChar+1 << std::endl;
  if(strncmp(mFileName+iChar+1,"root",4)==0){
    mDataFile = new LecroyFile(mFileName);
  }
  else{
    
    if(strncmp(mFileName+iChar+1,"mid.lz4",6)==0 || strncmp(mFileName+iChar+1,"mid.gz",6)==0){
      cout << "Opening V1730 file" << endl;
      mDataFile = new V1730File(mFileName, mChannel);
    }
    if(strncmp(mFileName+iChar+1,"mid.lz4",6)==0 || strncmp(mFileName+iChar+1,"txt",3)==0){
      cout << "Opening text file" << endl;
      mDataFile = new TextFile(mFileName, mChannel);
    }
    else{
      std::cerr << "cannot tell file type. Aborting " << std::endl;
      exit(1);
    };
    
  }	
  //if(quietMode == false) std::cout << "Waveform found : " << mDataFile->getWaveformCount() << std::endl;

  cout << "Prepare Vector of Pulses.." << endl;
  
  // --- prepare vector of pulses
  mNFitPulseMax = numPulsesFit;
  mNPulse=0;
  mNPulseMax= mNFitPulseMax*2;
  mPulse = new Pulse[mNPulseMax];
  mPulseIndex = new int[mNPulseMax];
  mPulseTime = new double[mNPulseMax];
  
  //cout<<"mNPulseMax "<<mNPulseMax<<endl;

  cout << "Prepare functions for pulse fitting.." << endl;

  // --- Function for pulse fitting
  //std::cout << "Find object :" << std::endl;
  mFExpGaus = (TF1*) gROOT->FindObject("FExpGaus");
  //mFExpGaus = (TF1*) gROOT->FindObjectAny("FExpGaus");
  //std::cout << "Deleting " << std::endl;
  if(mFExpGaus) mFExpGaus->Delete();
  //std::cout << "Deleted " << std::endl;
  mFExpGaus = new TF1("FExpGaus",FuncExpGausMulti,
		      0.,0.5e-6,4+mNFitPulseMax*2); // range reset later
  mFExpGaus->SetParLimits(1,0.,1e-6);// sigma
  mFExpGaus->SetParLimits(2,0.,1e-6);// tau
  strcpy(mFitOption,"QR0");
  
  mFExp2Gaus = (TF1*) gROOT->FindObject("FExp2Gaus");
  if(mFExp2Gaus  && gROOT->FindObjectAny(mFExp2Gaus->GetName())) mFExp2Gaus->Delete();
  mFExp2Gaus = new TF1("FExp2Gaus",FuncExp2GausMulti,
		       0.,0.5e-6,6+mNFitPulseMax*2); // range reset later
   mFExp2Gaus->SetParLimits(1,0.,1e-6);// sigma
   mFExp2Gaus->SetParLimits(2,0.,1e-6);// tau

  mFExpRise = (TF1*) gROOT->FindObject("FExpRise");
  if(mFExpRise  && gROOT->FindObjectAny(mFExpRise->GetName())) mFExpRise->Delete();
  mFExpRise = new TF1("FExpRise",FuncExpRiseMulti,
		   0.,0.5e-6,6+mNFitPulseMax*2); // range reset later
  mFExpRise->SetParLimits(1,1e-9,1e6);// rise time
  mFExpRise->SetParLimits(2,1e-9,1e6);// fast fall time
  mFExpRise->SetParLimits(4,1e-9,1e6);// slow fall time

 mFFit=mFExp2Gaus; //FBK LF

 // mFFit = mFExpRise;//VUV4 // hard coded selection of function
  
  // mFFit=mFExpGaus;

 cout << "WaveformProcessor Init done." << endl;
 
}

// ---
WaveformProcessor* WaveformProcessor::mInstanceForRootMacro=0;
int WaveformProcessor::mRun=-1;
int WaveformProcessor::mChannelUnique=-1;
WaveformProcessor* WaveformProcessor::instanceForRootMacro(const char* aSetupFileName, int aRun,  int numPulsesFit, int aChannel){
  if(!mInstanceForRootMacro){ //!!!! Must remain a singleton. Force init if needed somewhere else
    mInstanceForRootMacro = new WaveformProcessor(aSetupFileName, aRun, numPulsesFit, aChannel);
    mChannelUnique=aChannel;
  }
  else{
    if(mRun!=aRun || mChannelUnique!=aChannel){     
      mInstanceForRootMacro->clear();
      mInstanceForRootMacro->init(aSetupFileName, aRun, numPulsesFit, aChannel);
      mChannelUnique=aChannel;
    }
  }
  return mInstanceForRootMacro;
}

// ---
int WaveformProcessor::getWaveformCount(){
  return mDataFile->getWaveformCount();
}

//___________________________________________________________________
// ---
#include <time.h>
#include <sys/stat.h>
int WaveformProcessor::secDiffForRootMacros(const char* aFirstFileName,
					    const char* aFileName){
  struct stat buf;
  stat(aFirstFileName,&buf);
  time_t firstTime = buf.st_mtime;
  stat(aFileName,&buf);
  return  difftime(buf.st_mtime,firstTime);
}


//___________________________________________________________________
// ---
WaveformProcessor::~WaveformProcessor(){
  //std::cout << "WP Destructor... " << std::endl;
  clear();
  //std::cout << "WP Destructed... " << std::endl;
}
void WaveformProcessor::clear(){
  delete mDataFile;
  //std::cout << "Finished delet mDataFile " << std::endl;
  if(mHBaseline && gROOT->FindObjectAny(mHBaseline->GetName())) 
    mHBaseline->Delete();
  if(mFExpGaus && gROOT->FindObjectAny(mFExpGaus->GetName())) 
    mFExpGaus->Delete();
  if(mFExp2Gaus && gROOT->FindObjectAny(mFExp2Gaus->GetName())) 
    mFExp2Gaus->Delete();
  if(mFExpRise && gROOT->FindObjectAny(mFExpRise->GetName())) 
    mFExpRise->Delete();
  if(mPulse) delete[] mPulse;
  if(mPulseIndex) delete[] mPulseIndex;
  if(mPulseTime) delete[] mPulseTime;
  //std::cout << "Finished WP clear " << std::endl;
}

//___________________________________________________________________
// --- Access the waveform
bool WaveformProcessor::readNextWaveform(){
  mWF = mDataFile->getNextWaveform(); // index is not used!

 
  if(mWF){
    mWF->setQuietMode(quietMode);
    return true;
  }
  else{
    return false;
  }
}
void WaveformProcessor::readWaveform(int aIndex){

  std::cout << "Start Reading Waveform.." << std::endl;
  
  mWF = mDataFile->getWaveform(aIndex);

  std::cout << "mDataFile->getWaveform done." << std::endl;

  //mWF->setQuietMode(quietMode);  
  //std::cout << "Quietmode set." << std::endl;

  //Read event from file and fill the WF

  


  
/*
 
  //Idea di base creo un istogramma temporaneo con la media e opportuno ribinnaggio e poi lo associo a mWF

  //Valori originali
    
    double BinXmax = (mWF->GetNbinsX())*(mWF->GetBinWidth(1));
    double BinXmin= 0;
    int nBins = (mDataFile->getWFYNBins());
    
    nBins=500;
    std::cout<<"nBins:="<<nBins<<std::endl;
    
    //int number_bin_temp=500;
    
    std::cout<<"Preparing temporary histogram ..."<<std::endl;

    TH1F *temporary;
    temporary = new TH1F("TEMP","TEMP",nBins,BinXmin,BinXmax);
    
    float sum=0;

   
    //Corre sui bin di mWF e ogni 5 salva il valore
    for(int iBin=1; iBin <= (mWF->GetNbinsX()); iBin++){
        
        sum=sum+(mWF->GetBinContent(iBin));
        
        if(iBin%5==0){
        temporary->Fill(sum/5.0);
        std::cout<< "Mean"<<(sum/5.0)<<std::endl;
        sum=0;
        }
    }

  //temporary->Draw();
 */
    
}
void WaveformProcessor::readWaveform(const char* aHName){
  mWF = mDataFile->getWaveform(aHName);	
  mWF->setQuietMode(quietMode);										 
}
int WaveformProcessor::getCurrentWaveformIndex(){
  char ciwf[50];
  strcpy(ciwf,mWF->GetName()+7);
  return atoi(ciwf);
}

//___________________________________________________________________	
// --- baseline
int WaveformProcessor::processBaseline(){
  // --- Create histogram for baseline calculation
  mHBaseline = (TH1F*) gROOT->FindObjectAny("HBaseline");
  if(mHBaseline)  mHBaseline->Reset("C");//mHBaseline->Delete();
  else {
    int nBinsBaseline = (mAmpMax-mAmpMin)/mWFAmpBining; 
    //cout<<"mAmpMAx:="<<mAmpMax<<endl;
    //cout<<"mAmpMin:="<<mAmpMin<<endl;
    //cout<<"mWFAmpBining"<<mWFAmpBining<<endl;
if(mDataFile->getWFYLimitsKnown()){
      mAmpMax = mDataFile->getWFYMax();
      mAmpMin = mDataFile->getWFYMin();
      nBinsBaseline = mDataFile->getWFYNBins();
      // cout<<"mAmpMAX:="<<mAmpMax<<endl;
      //cout<<"mAmpMin:="<<mAmpMin<<endl;
      //cout<<"mWFAmBining:="<<nBinsBaseline<<endl;
 }
    mHBaseline = new TH1F("HBaseline","HBaseline",
			  nBinsBaseline,mAmpMin,mAmpMax);
  }


  //mAmpMin = mWF->GetMinimum();
  //mAmpMax = mWF->GetMaximum();
  //int nBinsBaseline = 500; !!! Why using this?
 

  for(int iBin=1; iBin <= mWF->GetNbinsX(); iBin++){
    mHBaseline->Fill(mWF->GetBinContent(iBin));
  }
  int blMostProbBin=0;
  int blMostProbCont=0;
  for(int iBin=1; iBin <= mHBaseline->GetNbinsX(); iBin++){
    if(mHBaseline->GetBinContent(iBin) > blMostProbCont){
      blMostProbCont = mHBaseline->GetBinContent(iBin);
      blMostProbBin = iBin;
    }				
  }
  mBmu=mHBaseline->GetBinLowEdge(blMostProbBin);

 /****************

13/02/2017

Giacomo Gallina

giacomo@triumf.ca

RMS chose for faster computation insted of fit gaussian peak. RMS consider the whole waveform so it overstimate the baseline RMS. Used only as reference
MU actually used in the following but should be reliable even withou fit

In general the rms and the std is not the same. The rms is commonly used as an estimator for the std because its often the only readily availalable quantity (it can be easily calculated from the data). The rms is an unbiased estimator of the std. If the underlying distribution is gaussian and systematic errors such as ADC quantization errors are small compared too the signal amplitude, but it can lead too terribly wrong results if these two requiremnts are not satisfied.

*/

  mBRMS=mHBaseline->GetRMS();
  
  //std::cout << "The baseline mu is: " << mBmu << std::endl;
  //std::cout << "The number of baseline bins is: " << nBinsBaseline << std::endl;
  //std::cout << mBmu << " " << mBRMS << std::endl;
  
  return mBRMS<25*mNoise1Bin; // Otherwise noisy. i.e. pick up
}

//___________________________________________________________________	
// --- pulse finding
int WaveformProcessor::findPulse(){
  mPulseFound = false;
  int nBinRiseTime = round(mRiseTimeSig*3/mWF->GetBinWidth(1));
  int nBinFallTime = round((mFallTimeTau*3+mFallTime2Tau)/mWF->GetBinWidth(1));

  //std::cout<<"findPulse()"<<std::endl;
  //std::cout<<"nBinFallTime_(2_width): "<<nBinFallTime<<std::endl;
  //std::cout<<"nBinRiseTime_(2_width): "<<nBinRiseTime<<std::endl;
  //std::cout<<"BinWidth: "<<(mWF->GetBinWidth(1))<<std::endl;
  //std::cout<<"mNPulseMax:"<<mNPulseMax<<std::endl;
  // std::cout<<std::endl;
  mNPulse=0;
  //std::cout << mBmu << " " << nBinRiseTime << " " << nBinFallTime << " " << mNoise1Bin << std::endl;
  for(int iBin=nBinRiseTime+1; iBin<=mWF->GetNbinsX(); iBin++){
    


     //>>Find Pulse

    // first  condition: pulse is more than 3 sigma away from baseline
    // second condition: pulse is negative and has a certain height
    
//if( (abs(mWF->GetBinContent(iBin-nBinRiseTime)-mBmu)>mDiffAmpThresh) &&

//    if((mPolarity*(mWF->GetBinContent(iBin)-mBmu)>6*mNoise1Bin) && 
//     (mPolarity*(mWF->GetBinContent(iBin)-
//    mWF->GetBinContent(iBin-nBinRiseTime))>5*mNoise1Bin) 

//    std::cout<<" "<<(mPolarity*(mWF->GetBinContent(iBin)-mBmu))<<" "<<6*mNoise1Bin<<" "<<mPolarity*(mWF->GetBinContent(iBin)-
//												    mWF->GetBinContent(iBin-nBinRiseTime))<<" "<<5*mNoise1Bin<<std::endl;

    if((mPolarity*(mWF->GetBinContent(iBin)-mBmu)>6.*mNoise1Bin) && // cannot cut too hard on this
       (mPolarity*(mWF->GetBinContent(iBin)-
		   mWF->GetBinContent(iBin-nBinRiseTime))>5.*mNoise1Bin) // 
       //(mPolarity*(mWF->GetBinContent(iBin-nBinRiseTime)-mBmu)>2.*mNoise1Bin) 
       ){
      mPulseFound = true;
      // std::cout<<""<<std::endl;
      //           std::cout<<"New pulse Found ....!"<<std::endl;
      //std::cout<<"First two condition fulfilled!"<<std::endl;
      

// >>> Find maximum amplitude without baseline subtraction

      double pulseAbsAmp=mPolarity*mWF->GetBinContent(iBin);
      //std::cout<<"pulseAbsAmp"<<pulseAbsAmp<<std::endl;
      iBin++;
      while(iBin <= mWF->GetNbinsX() && pulseAbsAmp <= mPolarity*mWF->GetBinContent(iBin)){
	// std:cout<<"eccomi!!"<<std::endl;	
	pulseAbsAmp=mPolarity*mWF->GetBinContent(iBin);
	iBin++;
      }
      //      pulseAbsAmp+=(mBmu*mPolarity);
      pulseAbsAmp*=mPolarity;
      // std::cout<<"iBin"<<iBin<<std::endl;
      //std::cout<<"pulseAbsAmp"<<pulseAbsAmp<<std::endl;

// >>> Calculate the baseline before the pulse

      double pulseBaseline=0.;
      double prePulseMin=0.;
      for(int iBinBase=iBin-nBinRiseTime-15; 
	  iBinBase<(iBin-nBinRiseTime-5); 
	  iBinBase++){
	pulseBaseline+=mWF->GetBinContent(iBinBase);
	if(prePulseMin>mWF->GetBinContent(iBinBase)*mPolarity){
	  prePulseMin=mWF->GetBinContent(iBinBase)*mPolarity;
	}
      }
      pulseBaseline/=10.;
      double pulseAmp=(pulseAbsAmp-pulseBaseline)*mPolarity;
      // std::cout << mWF->GetBinCenter(iBin) << " " << pulseBaseline << " " << pulseAbsAmp << " " <<  pulseAmp << std::endl;
      
// >>> Calculate the charge (for removing noise)

      double tCharge(0);
      for(int iBinCharge=iBin-nBinRiseTime-2;
	  iBinCharge<iBin+nBinFallTime;
	  iBinCharge++){
	tCharge+=mWF->GetBinContent(iBinCharge);
	tCharge-=mBmu;//pulseBaseline;
      }
      tCharge*=mPolarity;
      tCharge*=mWF->GetBinWidth(1);


      int iBinInt=iBin-nBinRiseTime-2;
      //std::cout <<" iBinInt: " <<iBinInt<<" Content: "<< mWF->GetBinContent(iBinInt) << " Baseline : " << pulseBaseline <<std::endl;

      //Prima anche qui c'era un 5, metto un 2 !!!!
        
        
      while(iBinInt<iBin+nBinFallTime && 
	    (mPolarity*(mWF->GetBinContent(iBinInt)-pulseBaseline))<1*mNoise1Bin) iBinInt++;
      int iPulseStart=iBinInt;      
      //std::cout << iPulseStart << " " << mWF->GetBinContent(iBinInt) << std::endl;

      //Qui ci metto un 2  invece che un 5 la condizione credo sia troppo forte
      iBinInt=iBin;
      while(iBinInt<iBin+nBinFallTime && 
	    (mPolarity*(mWF->GetBinContent(iBinInt)-pulseBaseline))>=2*mNoise1Bin){ 
	//std::cout<<"roba"<<(mWF->GetBinContent(iBinInt)-pulseBaseline)<<std::endl;
	//std::cout<<"ibinInit"<<iBinInt<<std::endl;	
	iBinInt++;
      }
      int iPulseEnd=iBinInt-1;      
      double pulseWidth=iPulseEnd-iPulseStart+1;

      /*
      std::cout <<"iPulseStart: "<<iPulseStart << std::endl;
      std::cout <<"iPulseEnd: "<<iPulseEnd << std::endl;
                  
      std::cout <<"GetBinCenter(iBin): "<<mWF->GetBinCenter(iBin) << std::endl;
      std::cout <<"GetBinCenter(iPulseStart): "<< mWF->GetBinCenter(iPulseStart) <<std::endl;
      std::cout<<"GetBinCenter(iPulseEnd): "<< mWF->GetBinCenter(iPulseEnd) << std::endl;
      std::cout << "pulseWidth: "<<pulseWidth << std::endl;
      std::cout<< "mBmu: "<<mBmu << std::endl;
      std::cout<<"tCharge: "<<tCharge<<std::endl;
      std::cout<<"pulseBaseline: "<< pulseBaseline << std::endl;
      std::cout<<"pulseAbsAmp: "<< pulseAbsAmp << std::endl;
      std::cout<<"pulseAmp: "<<pulseAmp << std::endl;
      std::cout<<"tCharge: "<<tCharge << std::endl;
      std::cout<<"3*sqrt(nBinRiseTime+nBinFallTime)*mNoise1Bin: "<<3*sqrt(nBinRiseTime+nBinFallTime)*mNoise1Bin << std::endl;
      std::cout<<"pulseWidth: "<<pulseWidth << std::endl;
      std::cout<<"(nBinRiseTime+nBinFallTime)/5: "<<(nBinRiseTime+nBinFallTime)/5. << std::endl;
      std::cout << "2nd check " <<  mWF->GetBinCenter(iBin) << " " 
      		<< pulseBaseline << " " 
      		<< tCharge << " " << 3*sqrt(nBinRiseTime+nBinFallTime)*mNoise1Bin << " " 
      	<< pulseWidth << "  " << (nBinRiseTime+nBinFallTime)/5. << " " 
		<< prePulseMin << " " << mBmu << " "<< (-prePulseMin-mBmu)<<" "<<mNoise1Bin<<std::endl;
      */
      

      
      if(iBin < mWF->GetNbinsX() && //This condition it is used to not keep pulses that are at the edge of the waveform
	 tCharge>3.*sqrt(nBinRiseTime+nBinFallTime)*mNoise1Bin && //FBK gain is small --> Reduced to 1* (was 3*)
	 (-prePulseMin-mBmu)<mNoise1Bin*4.// &&  //This condition fails for PMT usually ---> Reduced for FBK (was *4) --> pulse smaller
	 // pulseWidth>(nBinRiseTime+nBinFallTime)/13. // bit risky... --> Incresed for FBK (was /5.)
	 ){ 
	iBin--;
	//std::cout<<"saving .."<<std::endl;

	if(mNPulse<mNPulseMax){
	  //std::cout<<"Saving ...."<<std::endl;
	  mPulse[mNPulse].mTime=mWF->GetBinCenter(iBin);
	  mPulse[mNPulse].mAbsAmp=pulseAbsAmp;
	  mPulse[mNPulse].mAmp=pulseAmp;
	  mPulse[mNPulse].mBaseline=pulseBaseline;
	  mPulse[mNPulse].mQ=tCharge;
	  mPulse[mNPulse].mWidth=pulseWidth;
	  mPulseIndex[mNPulse]=mNPulse;	  
	  mNPulse++;
	}
      }
      iBin+=nBinFallTime;
    }
  }

  //questa ritorna sempre 1 se trova un numero di pulses minori di mNPulseMax
  return mNPulse<mNPulseMax;
}
// ---
TGraph* WaveformProcessor::getPulseGraph(){
  double* tTime = new double[mNPulse];
  double* tAmp = new double[mNPulse];
  for(int iPulse=0; iPulse<mNPulse; iPulse++){
    tTime[iPulse]=mPulse[iPulse].mTime;
    tAmp[iPulse]=mPulse[iPulse].mAbsAmp;
    //std::cout << iPulse << " " << tTime[iPulse] << " " << tAmp[iPulse] << std::endl;
  }
  TGraph* tG = new TGraph(mNPulse,tTime,tAmp); //TGraph(number of points, xArray, yArray)
  delete[] tTime;
  delete[] tAmp;
  return tG;
}

//
// --- Single PE checker
void WaveformProcessor:: calcSinglePulseTemplateChi2(){
  // >>> set the function parameter
  mFFit->SetParameter(0,mBmu);
  mFFit->SetParameter(1,mRiseTimeSig);
  mFFit->SetParameter(2,mFallTimeTau);
  mFFit->SetParameter(3,mTime2Frac);
  mFFit->SetParameter(4,mFallTime2Tau);
  mFFit->SetParameter(5,1);
  for(int iPulse=1; iPulse<mNFitPulseMax; iPulse++){
    mFFit->FixParameter(6+2*iPulse,0);
    mFFit->FixParameter(7+2*iPulse,0);
  }
  for(int iPulse=0; iPulse<mNPulse; iPulse++){
    // >>>
    mFFit->SetParameter(6,mPulse[iPulse].mQ);
    mFFit->SetParameter(7,mPulse[iPulse].mTime-mRiseTimeSig*3); 
    // >>>
    int iFirstBin = mWF->GetXaxis()->FindBin(mPulse[iPulse].mTime-mRiseTimeSig*5);
    if(iFirstBin<1) iFirstBin=1;
    int iLastBin = mWF->GetXaxis()->FindBin(mPulse[iPulse].mTime+mFallTimeTau*3)+mFallTime2Tau;
    if(iLastBin>mWF->GetNbinsX()) iLastBin=mWF->GetNbinsX();
    
    // >>>
    mPulse[iPulse].mSPTemplateChi2=0;
    for(int iBin=iFirstBin; iBin<=iLastBin; iBin++){
      double val = (mWF->GetBinContent(iBin)-mFFit->Eval(mWF->GetBinCenter(iBin)))/mNoise1Bin;
      mPulse[iPulse].mSPTemplateChi2+= (val*val);
    }
    mPulse[iPulse].mSPTemplateChi2/=(iLastBin-iFirstBin);
  }
}
  
//___________________________________________________________________
// --- fit Pulse
void WaveformProcessor::setFitOption(const char* aOption){	
  strcpy(mFitOption,aOption);
}
void WaveformProcessor::fit(int aType){

  //std::cout << mFallTime2Tau << std::endl;

  // aType = multiFit + 10* templateCheck . multiFit=2 for refiting 1 otherwise
  int tFitMulti = aType/100;
  int tTemplateCheck = (aType%100)/10;
  int tFitType = (aType%100)%10;
 
  //std::cout << tFitMulti << " " << tTemplateCheck << " " << tFitType << std::endl;

  // >>> Template calculation 
  // does not really work, so phase out
  if(tTemplateCheck )calcSinglePulseTemplateChi2();

  // >>> Set Histo error
  for(int iBin=1; iBin<=mWF->GetNbinsX(); iBin++){
    mWF->SetBinError(iBin,mNoise1Bin);
  }

  if(aType>0){
    int iPulse=0;
    while(iPulse<mNPulse){   

      //      std::cout<<" Setting .."<<iPulse<<" of .."<<mNPulse<<std::endl;
      
      // --- Skip fit if template comparison was good
 
      /*
      if(tTemplateCheck && mPulse[iPulse].mSPTemplateChi2<mMaxSPTemplateChi2){
	std::cout<<"Ciccio"<<std::endl;
	mPulse[iPulse].mFitChi2 = mPulse[iPulse].mSPTemplateChi2;
	mPulse[iPulse].mFitNDF = -1;
	iPulse++;
      }
      */
      if(0){

      }else{
	//std::cout<<"Fitting..."<<std::endl;
	// >>> Calculate fit limits
	// >>> Low limit given by first pulse
	
	//	mPulse[iPulse].mFitLowLimit = mPulse[iPulse].mTime-mRiseTimeSig*5-5*mWF->GetBinWidth(1); // fit range (low limit). Add 5 bins for baseline
	//     mPulse[iPulse].mFitLowLimit = mPulse[iPulse].mTime-mRiseTimeSig*10-5*mWF->GetBinWidth(1);// was 10 increase to 20

	mPulse[iPulse].mFitLowLimit = mPulse[iPulse].mTime-mRiseTimeSig*10-5*mWF->GetBinWidth(1); 

	if(mPulse[iPulse].mFitLowLimit<mWF->GetXaxis()->GetXmin()) 
	  mPulse[iPulse].mFitLowLimit = mWF->GetXaxis()->GetXmin();
	// >>> Pulses may be merged in a group to be fitted together. The limit is adjusted accordingly
	
	
	//Prima qui avevo mFallTimeTau*3 ora ci metto mFallTimeTau*5
	mPulse[iPulse].mFitHighLimit = 
	  mPulse[iPulse].mTime+mFallTimeTau*3.+mFallTime2Tau+mRiseTimeSig*5.+5*mWF->GetBinWidth(1); //Try to reduce mFallTimeTau*1 --> FBK has a long tail
	int tNPulseInGroup=1; // total number of pulses fitted together


	//Prima qui avevi mFallTimeTay*3 ora metto mFallTimeTau*5
	while(iPulse+tNPulseInGroup<mNPulse && 
	      mPulse[iPulse].mFitHighLimit>mPulse[iPulse+tNPulseInGroup].mTime){	
	  mPulse[iPulse].mFitHighLimit = mPulse[iPulse+tNPulseInGroup].mTime+mFallTimeTau*3+mFallTime2Tau+mRiseTimeSig*5+5*mWF->GetBinWidth(1);
	  tNPulseInGroup++;
	}
	if(mPulse[iPulse].mFitHighLimit>mWF->GetBinLowEdge(mWF->GetNbinsX())) 
	  mPulse[iPulse].mFitHighLimit = mWF->GetBinLowEdge(mWF->GetNbinsX());	
	for(int iFitPulse=iPulse; iFitPulse<iPulse+tNPulseInGroup; iFitPulse++){
	  mPulse[iFitPulse].mFirstPulseInGroup = iPulse;
	  mPulse[iFitPulse].mFitLowLimit=mPulse[iPulse].mFitLowLimit;
	  mPulse[iFitPulse].mFitHighLimit=mPulse[iPulse].mFitHighLimit;
	}
        
	//	std::cout<<"tNPulseInGroup: "<<tNPulseInGroup<<std::endl;
	// --- Set fit paremeters 
	mFFit->FixParameter(5,tNPulseInGroup);
	switch(tFitType){
	case 1: // fixed time constants but free baseline
	  
	  mFFit->ReleaseParameter(0);
	  mFFit->SetParameter(0,mPulse[iPulse].mBaseline);
	  //mFFit->FixParameter(0,mPulse[iPulse].mBaseline);
	  if(strcmp(mFFit->GetName(),"FExpRise")==0){
	    mFFit->FixParameter(1,mRiseTimeSig);
	    mFFit->FixParameter(2,mFallTimeTau-mRiseTimeSig);
	    mFFit->FixParameter(3,mTime2Frac);

	    //Changing 23/02/2018

            /*

              This changing allow for exp rise  to set the second expotential to zero in order to analyze cold data

              giacomo@triumf.ca

              Giacomo Gallina

	    */

            if(mFallTime2Tau!=0){

	      mFFit->FixParameter(4,mFallTime2Tau-mRiseTimeSig-mFallTimeTau);  

            }else{

              //Set paramenter 4 to zero
              mFFit->FixParameter(4,0);

            }




	  }
	  else{

	    mFFit->FixParameter(1,mRiseTimeSig);
	    mFFit->FixParameter(2,mFallTimeTau);
	    mFFit->FixParameter(3,mTime2Frac);
	    mFFit->FixParameter(4,mFallTime2Tau);

	  }
	  break;
	case 2: // free time constant but fixed baseline
	  mFFit->FixParameter(0,mPulse[iPulse].mBaseline);
	  mFFit->ReleaseParameter(1);
	  mFFit->SetParameter(1,mRiseTimeSig);
	  mFFit->ReleaseParameter(2);
	  mFFit->SetParameter(2,mFallTimeTau);
	  if(strcmp(mFFit->GetName(),"FExpRise")==0){
	    mFFit->SetParLimits(1,1e-9,1e6);
	    mFFit->SetParameter(2,mFallTimeTau-mRiseTimeSig);
	    mFFit->SetParLimits(2,1e-9,1e6);
	  }
	  if(mTime2Frac!=0.){

	    


	    mFFit->ReleaseParameter(3);
	    mFFit->SetParameter(3,mTime2Frac);
	    mFFit->SetParLimits(3,0.,1.);	    
	    mFFit->ReleaseParameter(4);
	    mFFit->SetParameter(4,mFallTime2Tau); // don't subtract because will be fitted anyway
	    mFFit->SetParLimits(4,1e-9,1e6); 
	    
    
	    mPulse[iPulse].mFitHighLimit = mPulse[iPulse].mTime+mFallTimeTau*3+mFallTime2Tau*3.+mRiseTimeSig*5+5*mWF->GetBinWidth(1);


	  }
	  else{


            /*

	      22/02/2018

              At low temperatures the fit function basically becomes a 1 exp instead of a 2 exp. Parameter 3 goes to 0 or to 1
              Basically only 1 exp is onvolved. So to analyze cold data I fixed p3 and p4 to 0
              This must be changed at higher temperatures where the shape is described with  2 exp

              Giacomo Gallina

              giacomo@triumf.ca



	    */


	    mFFit->FixParameter(3,0.);
	    mFFit->FixParameter(4,0.);

	    //Is it fine this limit here ? I think so ...

	    mPulse[iPulse].mFitHighLimit = mPulse[iPulse].mTime+mFallTimeTau*3+mFallTime2Tau*3.+mRiseTimeSig*5+5*mWF->GetBinWidth(1);

	  }
	  break;
	case 3: // free time constant but with tight limits
	  mFFit->ReleaseParameter(0);
	  mFFit->SetParameter(0,mPulse[iPulse].mBaseline);
	  // mFFit->FixParameter(0,mPulse[iPulse].mBaseline);
	  mFFit->ReleaseParameter(1);
	  mFFit->SetParameter(1,mRiseTimeSig);
	  mFFit->SetParLimits(1,mRiseTimeSig*0.1, mRiseTimeSig*10);
	  mFFit->ReleaseParameter(2);
	  mFFit->SetParameter(2,mFallTimeTau);
	  mFFit->SetParLimits(2,mFallTimeTau*0.1, mFallTimeTau*10);
	  if(mTime2Frac!=0.){
	    mFFit->ReleaseParameter(3);
	    mFFit->SetParameter(3,mTime2Frac);
	    mFFit->SetParLimits(3,0.,1.);
	    mFFit->ReleaseParameter(4);
	    mFFit->SetParameter(4,mFallTime2Tau);
	    mFFit->SetParLimits(4,mFallTime2Tau*0.1, mFallTime2Tau*10);
	  }
	  else{
	    mFFit->FixParameter(3,0.);
	    mFFit->FixParameter(4,0.);
	  }
	  break;
	}       
	for(int iFitPulse=0; iFitPulse<tNPulseInGroup; iFitPulse++){
	  mFFit->ReleaseParameter(6+2*iFitPulse);
	  mFFit->SetParameter(6+2*iFitPulse,(mPulse[iFitPulse+iPulse].mQ*mPolarity));
	  

	  if(mPolarity>0)  mFFit->SetParLimits(6+2*iFitPulse,0,1e10);

	  if(mPolarity<0) mFFit->SetParLimits(6+2*iFitPulse,-1e10,0);




	  /*std::cout << "limit " << iFitPulse 
	  	    << " " << mPulse[iFitPulse+iPulse].mQ 
	  	    << "<" << mNoise1Bin		   
	  	    << " " << mPulse[iPulse].mFitLowLimit
	  	    << "<" << mPulse[iFitPulse+iPulse].mTime
	  	    << "<" << mPulse[iPulse].mFitHighLimit << std::endl;
	  */
	  
	  //Prima era commentata
	  //mFFit->SetParLimits(6+2*iFitPulse,-20000.,0);
	 
	  mFFit->ReleaseParameter(7+2*iFitPulse);
	  mFFit->SetParameter(7+2*iFitPulse,mPulse[iFitPulse+iPulse].mTime-mRiseTimeSig*3.);
	  mFFit->SetParLimits(7+2*iFitPulse,mPulse[iPulse].mFitLowLimit,mPulse[iPulse].mFitHighLimit); // bad things happen when pulses are allowed in front. 
	}
	for(int iFitPulse=tNPulseInGroup; iFitPulse<mNFitPulseMax; iFitPulse++){
	  mFFit->FixParameter(6+2*iFitPulse,0);
	  mFFit->FixParameter(7+2*iFitPulse,0);
	}
	mFFit->SetRange(mPulse[iPulse].mFitLowLimit,mPulse[iPulse].mFitHighLimit);
   
	//Increase number of cicles to fit
	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);

	mWF->Fit(mFFit->GetName(),mFitOption);
	
	//>>> Store the fit information
	for(int iFitPulse=0; iFitPulse<tNPulseInGroup; iFitPulse++){
	  mPulse[iFitPulse+iPulse].mFitAmp = mFFit->GetParameter(6+2*iFitPulse);
	  mPulse[iFitPulse+iPulse].mFitTime = mFFit->GetParameter(7+2*iFitPulse);

	  //std::cout<<"Value"<<mFFit->Eval( mPulse[iFitPulse+iPulse].mFitTime)<<std::endl;

	  mPulse[iFitPulse+iPulse].mFitBaseline = mFFit->GetParameter(0);
	  mPulse[iFitPulse+iPulse].mFitRiseTime = mFFit->GetParameter(1);
	  mPulse[iFitPulse+iPulse].mFitFallTime = mFFit->GetParameter(2);
	  mPulse[iFitPulse+iPulse].mFitTime2Frac = mFFit->GetParameter(3);
	  mPulse[iFitPulse+iPulse].mFitFallTime2 = mFFit->GetParameter(4);
	  mPulse[iFitPulse+iPulse].mFitChi2 = mFFit->GetChisquare();
	  mPulse[iFitPulse+iPulse].mFitNDF = mFFit->GetNDF();
      mPulse[iFitPulse+iPulse].mFitQ = mFFit->Integral(mPulse[iFitPulse+iPulse].mFitLowLimit,mPulse[iFitPulse+iPulse].mFitHighLimit); //COMPUTE integral between define limits
	  mPulse[iFitPulse+iPulse].mRefitChi2 = -1;
	  mPulse[iFitPulse+iPulse].mRefitNDF = -1;
	}
	iPulse+=tNPulseInGroup;
      }
    }
  }

  //std::cout<<"tfitMult: "<<tFitMulti<<std::endl;
  if(tFitMulti){

    //NOT REORDER PULSES IN TIME, YOU WILL DO IT LATER
    reFitMulti();
  
}else{

  //SINCE YOU DON'T REFITTED, SORT PULSES                                  
  if(mNPulse>1){
                                                                                                                                                                                                       
    for(int iPulse=0; iPulse<mNPulse; iPulse++){
      mPulseIndex[iPulse]=iPulse;
      mPulseTime[iPulse]=mPulse[iPulse].mFitTime;
    }
    TMath::Sort(mNPulse, mPulseTime, mPulseIndex, false);
  }

}

}
	


// ---
void WaveformProcessor::reFitMulti(){
  //std::cout << "WaveformProcessor::reFitMulti" << std::endl;
  //int iPulse=0;
  //std::cout<<" refit ... cazzo"<<std::endl;
  int nPulseBeforeRefit=mNPulse;

  // std::cout<<" nPulseBeforeRefit"<<nPulseBeforeRefit<<std::endl;

  int newPulseAdded=0;
  for(int iPulse=0; iPulse<nPulseBeforeRefit; iPulse++){
    //    std::cout<<"Changing Main Pulse ..."<<std::endl;
    
    //std::cout<<"Number Fit Parameter"<<mFFit->GetNpar()<<std::endl;
 

    //Once you change the main pulse set fit parameter to 0

        
    /*
    for (int l=6;l<mFFit->GetNpar();l++){

      mFFit->ReleaseParameter(l);

    }

    */

    //std::cout <<" "<<mPulse[iPulse].mFirstPulseInGroup<<"Main Pulse "<<iPulse<<" "<<mPulse[iPulse].mFitChi2 << " " << mPulse[iPulse].mFitNDF << " " << mMinChi2ForRefit << std::endl;
    if(mPulse[iPulse].mFirstPulseInGroup==iPulse &&
       mPulse[iPulse].mFitNDF>0 &&
       //(mPulse[iPulse].mFitChi2/mPulse[iPulse].mFitNDF/mPulse[iPulse].mAmp*mMinChi2ForRefit)>3.5){ // refit this group
       (mPulse[iPulse].mFitChi2/mPulse[iPulse].mFitNDF)>mMinChi2ForRefit) {
      // std::cout<<"iPulse: "<<iPulse<<std::endl;   


      // >>> Set function parameters that will not change
      mFFit->SetRange(mPulse[iPulse].mFitLowLimit,mPulse[iPulse].mFitHighLimit);
      mFFit->ReleaseParameter(0);
      mFFit->SetParameter(0,mBmu);
      mFFit->FixParameter(1,mRiseTimeSig);
      mFFit->FixParameter(2,mFallTimeTau);
      mFFit->FixParameter(3,mTime2Frac);
      mFFit->FixParameter(4,mFallTime2Tau);
      if(strcmp(mFFit->GetName(),"FExpRise")==0){
	mFFit->FixParameter(2,mFallTimeTau-mRiseTimeSig);
	
	//mFFit->FixParameter(4,mFallTime2Tau-mRiseTimeSig-mFallTimeTau);


	if(mFallTime2Tau!=0){

	  mFFit->FixParameter(4,mFallTime2Tau-mRiseTimeSig-mFallTimeTau);

	}else{

	  //Set paramenter 4 to zero
	  mFFit->FixParameter(4,0);

	}



      } 


    
      // >>> Copy chi2 to refit
      mPulse[iPulse].mRefitChi2 = mPulse[iPulse].mFitChi2;
      mPulse[iPulse].mRefitNDF = mPulse[iPulse].mFitNDF;
      
      // >>> add pulses as long as the chi2 is too large and that the added oulse is not too close to a another pulse

      double tMinTimeDiff=1.; // minimum time between two pulses
      //      std::cout << "here " << std::endl;


      int tNFitPulse=1;
      int validRefit=1;// stop refitting if not much progress is being made
      while(tNFitPulse<(mNFitPulseMax-1) && //can add one more pulse due to n pulse fit limit
	    mNPulse<(mNPulseMax-1) &&  //can add one more pulse due to total n pulse limit
	    mPulse[iPulse].mRefitNDF>0 &&
	    //(mPulse[iPulse].mRefitChi2/mPulse[iPulse].mRefitNDF/mPulse[iPulse].mAmp*mMinChi2ForRefit)>3.5 && 
	    (mPulse[iPulse].mRefitChi2/mPulse[iPulse].mRefitNDF)>mMinChi2ForRefit &&
	    validRefit){
	
	//cout<<"Valid refitting !!!"<<endl;

	// >>> Set fit function starting parameters without additional pulse
	// Needed to figure out where to add the next pulse
	tNFitPulse=0;
	for(int iFitPulse=iPulse; iFitPulse<mNPulse; iFitPulse++){
	  if(mPulse[iFitPulse].mFirstPulseInGroup==iPulse){
	    mFFit->ReleaseParameter(6+2*tNFitPulse);
	    mFFit->SetParameter(6+2*tNFitPulse,mPulse[iFitPulse].mFitAmp);

	    if(mPolarity>0)  mFFit->SetParLimits(6+2*tNFitPulse,0,1e10);

	    if(mPolarity<0) mFFit->SetParLimits(6+2*tNFitPulse,-1e10,0);
	   
	    
	    mFFit->ReleaseParameter(7+2*tNFitPulse);
	    mFFit->SetParameter(7+2*tNFitPulse,
				     mPulse[iFitPulse].mFitTime);
	    mFFit->SetParLimits(7+2*tNFitPulse,
				     mPulse[iPulse].mFitLowLimit,
				     mPulse[iPulse].mFitHighLimit);
	    /*
	    std::cout <<"Refit pulse " << tNFitPulse << " "
	          << mPulse[iFitPulse].mFitAmp << " "
	          << mPulse[iFitPulse].mFitTime << " " 
	          <<  mPulse[iPulse].mFitLowLimit << " "
	          <<  mPulse[iPulse].mFitHighLimit << std::endl;
	    */
		  tNFitPulse++;
	  }
	}
	mFFit->FixParameter(5,tNFitPulse);
	
	// >>> Look for most likely position of next pulse
	int tLastBin = mWF->GetXaxis()->FindBin(mPulse[iPulse].mFitHighLimit)-1;
	int iMaxDiff;
	double maxDiff=0.;
	double newPulseTime=0.;
	double newPulseAmp=0;
	double tDiff;
	for(int iBin=mWF->GetXaxis()->FindBin(mPulse[iPulse].mFitLowLimit)+1;
	    iBin<=tLastBin; iBin++){
	  tDiff=mPolarity*(mWF->GetBinContent(iBin)-
			   mFFit->Eval(mWF->GetBinCenter(iBin)));
	  if(maxDiff<tDiff){
	    maxDiff=tDiff;
	    newPulseTime =  mWF->GetBinCenter(iBin);
	    newPulseAmp = mWF->GetBinContent(iBin);

	    //std::cout<<"newPulseTime: "<<newPulseTime<<std::endl;
	    //std::cout<<"newPulseAmp: "<<newPulseAmp<<std::endl;

	  }
	}
	
	// >>> Abort if too close to an existing pulse
	tMinTimeDiff=1.;
	for(int iFitPulse=0; iFitPulse<mNPulse; iFitPulse++){
	  if(mPulse[iFitPulse].mFirstPulseInGroup==iPulse){
	    tDiff=fabs(mPulse[iFitPulse].mFitTime-newPulseTime);
	    if(tMinTimeDiff>tDiff) tMinTimeDiff=tDiff;
	  }
	}

	validRefit=0;

	//previously here was 0.5e-9

	if(tMinTimeDiff>0.5e-9){// for Hamamatsu the pulse shape has a double peak. Had to reduce it for VUV3 from 4 to 1
	  // >>> Add one more pulse
	  mFFit->ReleaseParameter(6+2*tNFitPulse);
	  mFFit->SetParameter(6+2*tNFitPulse,newPulseAmp-mFFit->GetParameter(0));
	  //Da commentare la riga seguente nel caso era -mNoise1Bin e ora 0 culo
	  
	  if(mPolarity>0)  mFFit->SetParLimits(6+2*tNFitPulse,0,1e10);

	  if(mPolarity<0) mFFit->SetParLimits(6+2*tNFitPulse,-1e10,0);
	  
	  mFFit->ReleaseParameter(7+2*tNFitPulse);
	  mFFit->SetParameter(7+2*tNFitPulse,newPulseTime);
	  mFFit->SetParLimits(7+2*tNFitPulse,mPulse[iPulse].mFitLowLimit,mPulse[iPulse].mFitHighLimit);

	  /*
	  	     std::cout <<"Refit pulse second read.: " << tNFitPulse << " "
	      << mFFit->GetParameter(6+2*tNFitPulse) << " "
	       << mFFit->GetParameter(7+2*tNFitPulse) << " "
	        <<  mPulse[iPulse].mFitLowLimit << " "
	        <<  mPulse[iPulse].mFitHighLimit << std::endl;
	  */

	  tNFitPulse++;
	  mFFit->FixParameter(5,tNFitPulse);
	  mPulse[mNPulse].mFirstPulseInGroup=iPulse;
	  mNPulse++;
	  
	  // >>> Fit
	  mWF->Fit(mFFit->GetName(),mFitOption);	 	  

	  
	  if((mPulse[iPulse].mRefitChi2-mFFit->GetChisquare())/mPulse[iPulse].mRefitChi2>0.02){//chi2 has decreased by more than 2% (0.02)
	    newPulseAdded=1;
	    validRefit=1;
	    // >>> Copy fit information
	    int iFitParameter=0;
	    for(int iFitPulse=iPulse; iFitPulse<mNPulse; iFitPulse++){
	      if(mPulse[iFitPulse].mFirstPulseInGroup==iPulse){	      
		mPulse[iFitPulse].mFitAmp = mFFit->GetParameter(6+2*iFitParameter);
		mPulse[iFitPulse].mFitTime = mFFit->GetParameter(7+2*iFitParameter);
		mPulse[iFitPulse].mFitBaseline = mFFit->GetParameter(0);
		mPulse[iFitPulse].mFitRiseTime = mFFit->GetParameter(1);
		mPulse[iFitPulse].mFitFallTime = mFFit->GetParameter(2);
		mPulse[iFitPulse].mFitTime2Frac = mFFit->GetParameter(3);
		mPulse[iFitPulse].mFitFallTime2 = mFFit->GetParameter(4);
		mPulse[iFitPulse].mFitChi2 = mPulse[iPulse].mFitChi2;
		mPulse[iFitPulse].mFitNDF = mPulse[iPulse].mFitNDF;
		mPulse[iFitPulse].mRefitChi2 = mFFit->GetChisquare();
		mPulse[iFitPulse].mRefitNDF = mFFit->GetNDF();
        mPulse[iFitPulse].mFitQ = mFFit->Integral(mPulse[iFitPulse].mFitLowLimit,mPulse[iFitPulse].mFitHighLimit); //COMPUTE integral between define limits
		iFitParameter++;
	      }
	    }
	  }
	  else{

	    //std::cout<<"cazzo"<<std::endl;
	    mNPulse--;
	  }
	}
      }
    }   
   
    //Prima questa non era commentata
    // iPulse++;
  }
  // >>> Sort
  if(mNPulse>1 && newPulseAdded){

    //    cout<<"culo .."<<endl;

    for(int iPulse=0; iPulse<mNPulse; iPulse++){
      mPulseIndex[iPulse]=iPulse;
      mPulseTime[iPulse]=mPulse[iPulse].mFitTime;
    }
    TMath::Sort(mNPulse, mPulseTime, mPulseIndex, false);
  }
  
}

double WaveformProcessor::getTriggerTime(){
  return mDataFile->getTriggerTime();
}

TF1* WaveformProcessor::getFitFunction() {
  // include all pulses and full range of the waveform
  mFFit->SetLineColor(2);
  mFFit->SetRange(mWF->GetXaxis()->GetXmin(),
		      mWF->GetXaxis()->GetXmax());
  int tNFitPulse;
  if(mNPulse>mNFitPulseMax){
    std::cout << "WARNING, number of pulses exceed max allowed for function " << mNPulse << " vs " 
	      << mNFitPulseMax << " allowed" << std::endl;
    tNFitPulse = mNFitPulseMax;
  }
  else{
    tNFitPulse = mNPulse;
  }

  //  cout<<"tNFitPulse: "<<tNFitPulse<<endl;

  mFFit->SetParameter(5,tNFitPulse);
  for(int iPulse=0; iPulse<tNFitPulse; iPulse++){
    mFFit->SetParameter(6+2*iPulse,mPulse[iPulse].mFitAmp);
    mFFit->SetParameter(7+2*iPulse,mPulse[iPulse].mFitTime);
    std::cout<< "getFitFunction: " << mPulseIndex[iPulse] << " " 
    << mPulse[mPulseIndex[iPulse]].mFitTime << " " << mPulse[mPulseIndex[iPulse]].mFitAmp << std::endl;
  }
  for(int iPulse=tNFitPulse; iPulse<mNFitPulseMax; iPulse++){
    mFFit->FixParameter(6+2*iPulse,0);
    mFFit->FixParameter(7+2*iPulse,0);
  }
  return mFFit;
}
TF1* WaveformProcessor::getFitFunction2() {
	return mFFit;
}
TH1F* WaveformProcessor::getPulseHisto(){
	return pulseHisto;
}
void WaveformProcessor::setQuietMode(bool mode){
	quietMode = mode;
}


