

#include <cstring>
#include <inttypes.h>
#include <inttypes.h>
#include <stdlib.h>

#include "TextFile.h"
#include "Waveform.h"

#include "TFile.h"
#include "TKey.h"
#include "TROOT.h"
#include "TList.h"

#include <sstream>
#include <iomanip>

#include "TDataContainer.hxx"


TextFile::TextFile(const char* aFileName, int aChannel):DataFile(0),mChannel(aChannel),mWF(0){

  mFile.open(aFileName);
  std::cout << "(TextFile) opening file " << aFileName << std::endl;

  mWF = new Waveform(1300,0,1300);
  
}

TextFile::~TextFile() {
  mFile.close();
}


Waveform* TextFile::getNextWaveform(){
  return 0;
  //return mWF;  
}
  

Waveform* TextFile::getWaveform(int index){

  std::cout << "(TextFile) Start getWaveform.." << std::endl;
  
  //int Run, WFAmpBining ampMin ampMax noise riseTimeSig fallTimeTau Time2Frac fallTime2Tau Pol MinChi2ForRefit
  std::string FileName;
  std::string line;
  
  //std::size_t found = str.find(str2);

  int i=0,j=0;
  int evt = -1;
  bool read = false;
  while(getline(mFile, line)){
    std::size_t found = line.find("Event n.");
    if (found!=std::string::npos){
      //std::cout << "found at: " << found << '\n';
      evt = stoi(line.substr(8,12));
      //std::cout << evt << std::endl;
      if (evt == index){
	read = true;
      }
      else read = false;
    }
    if (read){
      //std::cout << "Found Event " << evt << std::endl;
      //std::cout << line << std::endl;

      size_t pos = 0;
      std::string col;
      std::string delimiter = "\t";
      double data[3];
      i=0;
      while ((pos = line.find(delimiter)) != std::string::npos) {
	col = line.substr(0, pos);	
	data[i] = atof(col.c_str());

	//std::cout << col << std::endl;
	
	i++;
	line.erase(0, pos + delimiter.length());
      }

      std::cout << data[0] << " " << data[1] << " " << data[2] << std::endl;
      
      //fill the WF
      mWF->SetBinContent(j,data[1]);
      j++;
    }
    
  }
  
  return mWF;

}
