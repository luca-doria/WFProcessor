/*******************************************************
* File reader interface
*
* History:
* v1	2011/11/25	Initial version (Kyle Boone)
*********************************************************/

#ifndef DATA_FILE_H
#define DATA_FILE_H

#include <time.h>
#include <sys/stat.h>

class Waveform;


class DataFile {
  protected:
    int mWaveformCount;
    bool mGood;
    int mRun;
    double mTriggerTime;

    int mWFYNBins;
    double mWFYMax;
    double mWFYMin;

    time_t mFileCreationTime;
   

  public:
    // General functions
    
 DataFile(int run) : mGood(true), mWaveformCount(0), mRun(run), mTriggerTime(0), mWFYNBins(-1) {}

   virtual  int getWaveformCount() { return mWaveformCount; }
    int getRun() { return mRun; }
    double getTriggerTime(){return mTriggerTime;}

    // If the file is not opened properly, or if there are errors reading it
    // this will return false.
    bool isGood() { return mGood; }

    virtual ~DataFile() {};

    // Interface to be implemented by subclasses

    // Retrieve the waveform at a given index
    virtual Waveform* getWaveform(int index){return 0;}
    virtual Waveform* getWaveform(const char* aName){return 0;}
    virtual Waveform* getNextWaveform()=0;
		
    // Set the format for the names of waveforms that will be retrieved.
    virtual void setWaveformNameFormat(const char* format) {};

    // when the information about the Y scale of the waveform is stored in the data
    int getWFYLimitsKnown(){return mWFYNBins>0;}
    double getWFYMax(){return mWFYMax;}
    double getWFYMin(){return mWFYMin;}
    int getWFYNBins(){return mWFYNBins;}

    time_t getFileCreationTime(){return mFileCreationTime;}

    static int timeDifferenceBetweenFiles();//(const char* aFirstFileName,
    //const char* aFileName);
};




#endif // include guard
