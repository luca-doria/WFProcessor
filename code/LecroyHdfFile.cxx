#include <iostream>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "sys/stat.h"

#include "LecroyHdfFile.h"
#include "Waveform.h"


#ifdef __CINT__
#include "H5Cpp.h"
#endif

LecroyHdfFile::LecroyHdfFile(const char* aFileName, int aChannel):DataFile(0){
  struct stat buf;
  stat(aFileName,&buf);
  mFileCreationTime = buf.st_mtime;
 
  mFile = new H5File(aFileName,H5F_ACC_RDONLY );
  char channelName[100];
  sprintf(channelName,"channel%i",aChannel);
  mWFDataSet = mFile->openDataSet(channelName);
  
  // --- Read attribute from header. All attribute shown in lecroy manual are available
  // >>> Attribute were written in the channelx "bank"
  Attribute tAttr;
  tAttr = mWFDataSet.openAttribute("horiz_interval");
  float tDeltaTime;
  tAttr.read(PredType::NATIVE_FLOAT,&tDeltaTime);
  tAttr = mWFDataSet.openAttribute("horiz_offset");
  float tStartTime;
  tAttr.read(PredType::NATIVE_FLOAT,&tStartTime);
  tAttr = mWFDataSet.openAttribute("vertical_gain");
  tAttr.read(PredType::NATIVE_FLOAT,&mDeltaV);
  tAttr = mWFDataSet.openAttribute("vertical_offset");
  tAttr.read(PredType::NATIVE_FLOAT,&mOffsetV);
 
  
  // --- Read the waveform data from the bank
  // >>> Check the bank size
  mWFDataSpace = mWFDataSet.getSpace();
  hsize_t dataDim[2];
  int nDataDim = mWFDataSpace.getSimpleExtentDims(dataDim, NULL);
  mWaveformCount = dataDim[0];
  mNBins = dataDim[1];
  std::cout << "Opening " << aFileName << std::endl;
  std::cout << "channel " << aChannel << "  bank. nEvents= " << mWaveformCount 
	    << " | mNBins= " << mNBins 
	    << " | vertical interval = " << mDeltaV 
	    << " | vertical offset = " << mOffsetV 
	    << std::endl;  
  mWFYMin = mOffsetV-mDeltaV*pow(2,7);
  mWFYMax = mOffsetV+mDeltaV*pow(2,7);
  mWFYNBins = (int) pow(2,8);

  // >>> select memory for 1 event.
  mWFDataOffset[0] = 0; // set to event zero for first event
  mWFDataOffset[1] = 0;
  mWFDataCount[0] = 1; // read events one by one
  mWFDataCount[1] = mNBins;	
  // >>> create memory space for event
  hsize_t nDim=1;
  hsize_t dimOut[1];
  dimOut[0] = mNBins;
  mWFMemSpace = new DataSpace( nDim, dimOut); // not sure what is the difference between this
  mWFDataBuffer = new int[mNBins];             // and this. Why twice
  
  // --- create waveform
  mWF = new Waveform(mNBins,tStartTime,tStartTime+mNBins*tDeltaTime);

  // --- Read time data set if present
  mTrigTimeDataSetAvailable=1;
  try{
    mTimeDataSet = mFile->openDataSet("trigtime");
  }
  // catch failure caused by the H5File operations
  catch( FileIException error ){
    std::cout << std::endl << std::endl << std::endl 
	      << "No time data set. Must be an older file. You can ignore the earlier error messages but the trigger timing information will not be available" << std::endl;
    mTrigTimeDataSetAvailable=0;
  }
  if(mTrigTimeDataSetAvailable){
    mTimeDataSpace = mTimeDataSet.getSpace(); // space in file
    nDataDim = mTimeDataSpace.getSimpleExtentDims(dataDim, NULL);
    std::cout << "time bank. n Dim=2? " << nDataDim << " nEvents = " << dataDim[0] << " dim check =2? " << dataDim[1] << std::endl;
    //
    mTimeDataCount[0]=1; // read events one by one
    mTimeDataCount[1]=2;
    // >>> space in memory  
    nDim=1;
    dimOut[0] = 2;
    mTimeMemSpace = new DataSpace(nDim, dimOut); // only want one value 
  }
}


Waveform* LecroyHdfFile::getNextWaveform(){
  getWaveform(mWFDataOffset[0]);
  mWFDataOffset[0]++;
  return mWF;
}

Waveform* LecroyHdfFile::getWaveform(int index){
  //Attribute tAttr;
  //tAttr = mWFDataSet.openAttribute("trigtime_array");
  //char tTime[8];
  //std::cout << "Attribute "
    //	    << tAttr.getInMemDataSize() << std::endl;
  //tAttr.read(PredType::NATIVE_ULONG,&tTime);
  //std::cout << atof(tTime) << std::endl;
  
  mWFDataOffset[0]=index;
  mWFDataSpace.selectHyperslab( H5S_SELECT_SET, mWFDataCount, mWFDataOffset );
  mWFDataSet.read(mWFDataBuffer, 
		  PredType::NATIVE_INT, 
		  *mWFMemSpace, 
		  mWFDataSpace );
  for(int iBin=0; iBin<mNBins; iBin++){
    mWF->SetBinContent(iBin+1,mOffsetV+mWFDataBuffer[iBin]*mDeltaV);
  }

  if(mTrigTimeDataSetAvailable){
    mTimeDataSpace.selectHyperslab( H5S_SELECT_SET, mTimeDataCount, mWFDataOffset);
    mTimeDataSet.read(mTimeDataBuffer, 
		      PredType::NATIVE_DOUBLE, 
		      *mTimeMemSpace, 
		      mTimeDataSpace);
    mTriggerTime = mTimeDataBuffer[0];//+mTimeDataBuffer[1];
    //std::cout << mTimeDataBuffer[0] << " " <<  mTimeDataBuffer[1] <<  std::endl;
  }
  return mWF;
}
	

LecroyHdfFile::~LecroyHdfFile(){
	delete mFile;
	delete mWFMemSpace;
	delete[] mWFDataBuffer;
	if(mTrigTimeDataSetAvailable) delete mTimeMemSpace;
	delete mWF;
}
