
ROOTANADIR = $(ROOTANASYS)/include

GFLAGS = -g -O0

CFLAGS = -c -fPIC $(GFLAGS) $(shell root-config --cflags) -I. -I$(ROOTANASYS)/include 
LFLAGS = $(GFLAGS) $(shell root-config --libs) -lMinuit -lz -lm $(ROOTANALIBS)

all: PF lib

PF : bin/PulseFinding.exe

lib : lib/libSipmAnalysis.so

lib/libSipmAnalysis.so : obj/Waveform.o obj/WaveformProcessor.o  obj/LecroyFile.o obj/V1730File.o obj/TextFile.o obj/DataFile.o obj/RootDict.o 
	g++ -shared -fPIC -o lib/libSipmAnalysis.so $^ $(LFLAGS)

bin/PulseFinding.exe : obj/PulseFinding.o  obj/LecroyFile.o obj/V1730File.o obj/TextFile.o obj/WaveformProcessor.o obj/Waveform.o obj/DataFile.o 
	g++ -o bin/PulseFinding.exe $^ $(LFLAGS)

obj/RootDict.o : code/Waveform.h code/LecroyFile.h code/V1730File.h code/TextFile.h code/WaveformProcessor.h code/DataFile.h code/Linkdef.h
	rootcint -f code/RootDict.cxx -c -I$(ROOTANADIR) $^
	g++ -c -o $@ code/RootDict.cxx $(CFLAGS)

obj/%.o: code/%.cxx
	g++ -c -o $@ $< $(CFLAGS)

clean:
	rm -f obj/* bin/* lib/*.so
