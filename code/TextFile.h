
#ifndef TEXT_FILE_H 
#define TEXT_FILE_H

#include "DataFile.h"
#include <fstream>

class Waveform;
class TDataContainer;

class TextFile : public DataFile {
 public:
  TextFile(const char* aFileName, int aChannel);
  Waveform* getNextWaveform();
  Waveform* getWaveform(int index);

  virtual ~TextFile();
  
private:

  TDataContainer *mDataContainer;
  Waveform* mWF;
  int mChannel;
  int mIndex;
  std::ifstream mFile;
  
}; 


#endif


