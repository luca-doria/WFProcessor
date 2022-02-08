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

#ifndef V1730_FILE_H 
#define V1730_FILE_H

//#include <vector>
//#include <stdint.h>

#include "DataFile.h"

class Waveform;
class TDataContainer;

class V1730File : public DataFile {
 public:
  V1730File(const char* aFileName, int aChannel);
  Waveform* getNextWaveform();
  Waveform* getWaveform(int index);

  virtual ~V1730File();
  
private:

  TDataContainer *mDataContainer;
  Waveform* mWF;
  int mChannel;
  int mIndex;

}; 


#endif


