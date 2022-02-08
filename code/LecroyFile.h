/*******************************************************
* Tektronix WFM File reader
*
* History:
* v1	2011/08/17	Initial support for WFM003 mFiles (Kyle Boone)
* v2	2011/11/25	Integration with DataFile class (Kyle Boone)
*
* TODO:
* add support for LECROY001 and LECROY002 mFiles
* handle endianness
*********************************************************/

#ifndef LECROY_FILE_H 
#define LECROY_FILE_H

#include <fstream>

#include "DataFile.h"

class Waveform;
class TFile;
//class TObjLink;
class TIter;

class LecroyFile : public DataFile {
private:
	TFile* mFile;
	//TObjLink* mKeyLink;
	TIter* mKeyIter;
	int mCurIndex;
  // waveform data
  //Waveform* mWaveform;

public:
  LecroyFile(const char* aFileName);
	Waveform* getNextWaveform();
	//Waveform* getWaveform(const char* aHName);
  Waveform* getWaveform(int index);
  virtual ~LecroyFile();

private:
  // internal helper functions
  //bool openFile(const char* filename, int index);

}; 


#endif


