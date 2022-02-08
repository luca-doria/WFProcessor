
/*******************************************************

*
0 1000 0
*********************************************************/
#include <cstring>
#include <inttypes.h>
#include <iostream>
#include <inttypes.h>
#include <stdlib.h>

#include "TV1730RawData.hxx"
#include "V1730File.h"
#include "Waveform.h"

#include "TFile.h"
#include "TKey.h"
#include "TROOT.h"
#include "TList.h"

#include "TDataContainer.hxx"

V1730File::V1730File(const char* aFileName, int aChannel):DataFile(0),mChannel(aChannel),mWF(0){

}

V1730File::~V1730File() {

}


Waveform* V1730File::getNextWaveform(){
  return 0;
  //return mWF;  
}
  

Waveform* V1730File::getWaveform(int index){

  return 0;
  //..return mWF;
}
