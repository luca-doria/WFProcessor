#include "DataFile.h"

// going around macro limitations
int DataFile::timeDifferenceBetweenFiles(){//const char* aFirstFileName,
//const char* aFileName){
  //struct stat buf;
  //stat(aFirstFileName,&buf);
  //time_t firstTime = buf.st_mtime;
  //stat(aFileName,&buf);
  //return  difftime(buf.st_mtime,firstTime);

  return 1;
}
