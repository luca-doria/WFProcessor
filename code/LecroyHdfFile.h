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

#ifndef LECROYHDF_FILE_H 
#define LECROYHDF_FILE_H

#include "DataFile.h"

#ifndef __CINT__
#include "H5Cpp.h"
using namespace H5;
#else
class H5File;
class DataSet;
class DataSpace;
class hsize_t;
#endif /* __CINT __ */

class Waveform;

class LecroyHdfFile : public DataFile {
 public:
  LecroyHdfFile(const char* aFileName, int aChannel=1);
  Waveform* getNextWaveform();
  Waveform* getWaveform(int index);
  virtual ~LecroyHdfFile();
  
private:
  H5File* mFile;
  DataSet mWFDataSet;
  DataSpace mWFDataSpace;
  hsize_t mNBins;
  hsize_t mWFDataOffset[2];
  hsize_t mWFDataCount[2];
  DataSpace* mWFMemSpace;
  int* mWFDataBuffer;
  
  int mTrigTimeDataSetAvailable;
  DataSet mTimeDataSet;
  DataSpace mTimeDataSpace;
  hsize_t mTimeDataCount[2];
  DataSpace* mTimeMemSpace;
  double mTimeDataBuffer[2];
  
  float mOffsetV;
  float mDeltaV;
  
  Waveform* mWF;
  
}; 


#endif


