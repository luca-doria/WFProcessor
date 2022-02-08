/*******************************************************
* Tektronix WFM File reader
*
* History:
* v1	2011/08/17	Initial support for WFM003 mFiles (Kyle Boone)
* v2	2011/11/25	Integration with DataFile class (Kyle Boone)
*
* TODO:
* add support for WFM001 and WFM002 mFiles
* handle endianness
*********************************************************/

#include "LecroyFile.h"
#include "Waveform.h"
#include "TFile.h"
#include "TKey.h"
#include "TROOT.h"
#include "TList.h"

//#include <cstring>
//#include <inttypes.h>
#include <iostream>



// --- class functions

LecroyFile::LecroyFile(const char* aFileName):DataFile(0){
	// used with macros. Irrelevant with compiled code
	mFile = (TFile*) gROOT->GetListOfFiles()->FindObject(aFileName);
	if(!mFile){
		std::cout << "opening file " << aFileName << std::endl;
		mFile = new TFile(aFileName);
	}
	mKeyIter = new TIter(mFile->GetListOfKeys());
	//mKeyLink = mFile->GetListOfKeys()->FirstLink();
	mCurIndex=0;
	mWaveformCount=mFile->GetNkeys();	
 
}



LecroyFile::~LecroyFile() {

}

Waveform* LecroyFile::getNextWaveform(){
	Waveform* tWF = (Waveform*) ((TKey*) mKeyIter->Next())->ReadObj();
	mCurIndex++;
	//std::cout << tWF->GetName() << " " << mCurIndex << std::endl;
	return tWF;
}

//Waveform* LecroyFile::getWaveform(const char* aHName){
//return (Waveform*) mFile->Get(aHName);
//}


Waveform* LecroyFile::getWaveform(int index) {



  return 0;//(Waveform*) mFile->Get(wfName);
}




