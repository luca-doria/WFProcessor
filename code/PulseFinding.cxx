

/*
 
 TRIUMF PULSE FINDING AND FITTING ROUTINE
 
 
 Le tre seguenti variabili devono essere date al fine di calcolare correttamente la time distribution e quindi creare l'istogramma dei tempi
 
 float sigma_1pe_pulse;
 float threshold;


Example Usage:

The parameters are:

Order Number in RunInfo.txt - Type Fit - Number of Waveforms to Fit - Directory of the RunInfo.txt File

./bin/PulseFinding.exe 8933 101 100000 0 VUV4PDEApr2018


Note:

-Usually the visualization macro event number is the Pulseifinder number -1. This can however change if some events are discarded.
 Check and compare always the trigger time. 

History:

-26/02/2019 Added check to see if pulse polarity and pulse time is consistent
-27/02/2019 Added a check to not save bad waveforms
 
 */



#include <sys/stat.h>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "LecroyFile.h"
#include "V1730File.h"
#include "WaveformProcessor.h"
#include "Waveform.h"

#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH1.h"

// Global variables


int main(int argc, char** argv){
    int aRun = argc>1 ? atoi(argv[1]) : 0;
    int aFitType = argc>2 ? atoi(argv[2]) : 0;
    int aNEventMax = argc>3 ? atoi(argv[3]) : 100000000;
    int aChannel = argc>4 ? atoi(argv[4]) : 1;
    char aRunInfo[100];
    if(argc>5) sprintf(aRunInfo,"ntp/%s/RunInfo.txt",argv[5],argv[5]);
    else strcpy(aRunInfo,"RunInfo.txt");
    
    //Time distribution Analysis
    int number_pulse;
    float **pulses_values;
    float time;
    float pulse_fit_amp;
    float minimum;
    int pulse_found;
    
    //Variable of amplitude, input
    float sigma_1pe_pulse=50;
    float threshold=6650;
    
    //Third parameter number of pulses allowed. Usally 5 

    std::cout << "Instantiating WF processor.." << std::endl;

    WaveformProcessor wfProc(aRunInfo,aRun,15,aChannel);

    std::cout << "WF processor instantiated.." << std::endl;
    
    // --- Open output file
    //>>> Strip directory name from file name
    int slashIndex=0;
    for(int index=0; index<strlen(wfProc.getFileName()); index++){
        if(strncmp(wfProc.getFileName()+index,"/",1)==0) slashIndex=index;
    }
    char outFileName[200];
    if(argc>5){
        sprintf(outFileName,"ntp/%s/%s.Ch%iFit%i",argv[5],wfProc.getFileName()+slashIndex,wfProc.getChannel(), aFitType);
    }
    else{
        sprintf(outFileName,"ntp/%s.Ch%iFit%i",wfProc.getFileName()+slashIndex,wfProc.getChannel(), aFitType);
    }
    //sprintf(outFileName,"/home/huth/Desktop/nEXO/testdata%s.fanat%i",wfProc.getFileName()+slashIndex,aFitType);
    std::cout << outFileName <<  " " << wfProc.getBaselineRMS() << slashIndex << std::endl;
    
    //Old version
    /*
    TFile outFile(outFileName ,"RECREATE");
    TNtuple ntp("ntp","ntp",
                "evt:tt:blmu:blRMS:np:pbl:paa:pa:pt:pq:pw:tchi2:fa:ft:frt:fft:ff2:fft2:fblmu:fchi2:ndf:frchi2:frndf:fcharge");
    TNtuple ntpE("ntpE","ntpE","evt:tt:blmu:blRMS:Amin:tmin:Amax:tmax:np:qL:qD");
    TNtuple ntptime("ntptime","ntptime","evt:td");
    */
    
    TDirectory* procDir = gDirectory;

    //New version
    TFile *outFile = new TFile(outFileName,"RECREATE");
                                                                                                                               
    TNtuple *ntp=new TNtuple("ntp","ntp","evt:tt:blmu:blRMS:np:pbl:paa:pa:pt:pq:pw:tchi2:fa:ft:frt:fft:ff2:fft2:fblmu:fchi2:ndf:frchi2:frndf:fcharge");                                            
    TNtuple *ntpE=new TNtuple("ntpE","ntpE","evt:tt:blmu:blRMS:Amin:tmin:Amax:tmax:np:qL:qD");                                                                                         
    TNtuple *ntptime=new TNtuple("ntptime","ntptime","evt:td");    

    float ntpCont[100];
    int skippedcount = 0;

    procDir->cd();

    // ---
    int nEvent = wfProc.getWaveformCount();
    int number_waveform_discarded=0;


    std::cout<< "Number of Event: "<<nEvent<<std::endl;
    
    if(nEvent>aNEventMax || nEvent==-1) nEvent=aNEventMax;
    int iEvent=0;
    while(iEvent<nEvent && wfProc.readNextWaveform()){
      //for(int iEvent=0; iEvent<nEvent; iEvent++){
      //wfProc.readNextWaveform();
      iEvent++;
      int skipFit=0;
      if(wfProc.processBaseline() && wfProc.findPulse()){

	wfProc.fit(aFitType);
      }
      else{
	skippedcount++;
	skipFit=1;
      }

      Waveform* curWF = wfProc.getCurrentWaveform();
      ntpE->Fill(iEvent,wfProc.getTriggerTime(),
		 wfProc.getBaselineMu(),wfProc.getBaselineRMS(),
		 curWF->GetMinimum(),curWF->GetBinLowEdge(curWF->GetMinimumBin()),
		 curWF->GetMaximum(),curWF->GetBinLowEdge(curWF->GetMaximumBin()),
		 wfProc.getNPulse(),
		 curWF->Integral(2400,2700),curWF->Integral(1000,1300));
      //std::cout << "Event: " << iEvent <<"\t nPulses: "
      //          << wfProc.getNPulse() << std::endl;
      
      // if(iEvent==9173) std::cout<<"Number pulses ..."<<number_pulse<<std::endl;
      
      //Trovo numero di Pulses
      number_pulse=wfProc.getNPulse();
      

      if(iEvent==0) std::cout<<"Number pulses ..."<<number_pulse<<std::endl;    

      //Alloco un array del tipo pulses[Number_of_pulses][Information_per_pulses]
      //solo se il numero di pulses è diverso da zero
      if(number_pulse!=0){
	
	
	pulses_values = (float**)malloc(number_pulse *sizeof(float*));
	for(int i=0; i<number_pulse; i++){
	  *(pulses_values+i) = (float*)malloc(23 *sizeof(float));
	}
	
        
      }
        
      //std::cout<<"changing event ..."<<std::endl;

      
      int troubles=0;
      for(int iPulse=0; iPulse<wfProc.getNPulse(); iPulse++){


        // CHECK IF THE PULSES ARE NON ORDERED IN TIME                                                                                                                                                                                                                                                               
        if(iPulse!=0){

          if(wfProc.getFitTime(iPulse)<wfProc.getFitTime(iPulse-1)){

	    std::cout<<"Pulses are not ordered in time!!! "<<std::endl;

	    std::cout<<"iEvent: "<<iEvent<<std::endl;

	    // std::cout<<"Fit Amplitude "<<wfProc.getFitAmplitude(iPulse)<<std::endl;

	    number_waveform_discarded++;

	    troubles=1;

	    break;  

          }

        }


	//CHECK IF THE FIT PULSE POLARITY IS OK                                                                                                                                                                                                                                                                      
        if(wfProc.getpulsepolarity()*wfProc.getFitAmplitude(iPulse)<0){

	  //std::cout<<"Fit Amplitude "<<wfProc.getFitAmplitude(iPulse)<<std::endl;

	  std::cout<<"Fit Pulse Polarity is wrong !!! "<<std::endl;

	  std::cout<<"iEvent: "<<iEvent<<std::endl;


	  number_waveform_discarded++;

	  troubles=1;

          break;

        }



      }


      if(troubles==1){

	std::cout<<"Skipping current waveform ..."<<std::endl;

	continue;

      }

      for(int iPulse=0; iPulse<wfProc.getNPulse(); iPulse++){
      

	if(iPulse!=0){

	//ADDITIONAL USELESS CHECK
	if(wfProc.getFitTime(iPulse)<wfProc.getFitTime(iPulse-1)){

	  std::cout<<"Pulses are not ordered in time!!! "<<std::endl;

	  std::cout<<"iEvent: "<<iEvent<<std::endl;

	  // std::cout<<"Fit Amplitude "<<wfProc.getFitAmplitude(iPulse)<<std::endl;                                                                                                                                                                                                                               
	  return 0;

	}

      }


      //CHECK IF THE FIT PULSE POLARITY IS OK                                                                                                                                                                                                                                                                     \
                                                                                                                                                                                                                                                                                                                     
      if(wfProc.getpulsepolarity()*wfProc.getFitAmplitude(iPulse)<0){

	//std::cout<<"Fit Amplitude "<<wfProc.getFitAmplitude(iPulse)<<std::endl;                                                                                                                                                                                                                                  

	std::cout<<"Fit Pulse Polarity is wrong !!! "<<std::endl;

	std::cout<<"iEvent: "<<iEvent<<std::endl;

	return 0;

      }


	//Salvo informazioni nell'ntupla                                                                                                                                                   
        ntpCont[0]=iEvent;
        ntpCont[1]=wfProc.getTriggerTime();
	ntpCont[2]=wfProc.getBaselineMu();
	ntpCont[3]=wfProc.getBaselineRMS();
	ntpCont[4]=wfProc.getNPulse();
	ntpCont[5]=wfProc.getPulseBaseline(iPulse);
	ntpCont[6]=wfProc.getPulseAbsAmplitude(iPulse);
	ntpCont[7]=wfProc.getPulseAmplitude(iPulse);
	ntpCont[8]=wfProc.getPulseTime(iPulse);
	ntpCont[9]=wfProc.getPulseCharge(iPulse);
	ntpCont[10]=wfProc.getPulseWidth(iPulse);
	ntpCont[11]=wfProc.getSPTemplateChi2(iPulse);
	ntpCont[12]=wfProc.getFitAmplitude(iPulse);
	ntpCont[13]=wfProc.getFitTime(iPulse);
	ntpCont[14]=wfProc.getFitRiseTime(iPulse);
	ntpCont[15]=wfProc.getFitFallTime(iPulse);
	ntpCont[16]=wfProc.getFitTime2Frac(iPulse);
	ntpCont[17]=wfProc.getFitFallTime2(iPulse);
	ntpCont[18]=wfProc.getFitBaseline(iPulse);
	ntpCont[19]=wfProc.getChi2(iPulse);
	ntpCont[20]=skipFit? -1 : wfProc.getNDF(iPulse);
	ntpCont[21]=wfProc.getChi2Refit(iPulse);
	ntpCont[22]=wfProc.getNDFRefit(iPulse);
	ntpCont[23]=wfProc.getFitIntegral(iPulse);
	//std::cout<<"ntpCont[23]: "<<ntpCont[23]<<std::endl;
	ntp->Fill(ntpCont);
        
	//Salvo informazione nell'array doppio
        
        if(number_pulse!=0){
	  
	  pulses_values[iPulse][0]=iEvent;
	  pulses_values[iPulse][1]=wfProc.getTriggerTime();
	  pulses_values[iPulse][2]=wfProc.getBaselineMu();
	  pulses_values[iPulse][3]=wfProc.getBaselineRMS();
	  pulses_values[iPulse][4]=wfProc.getNPulse();
	  pulses_values[iPulse][5]=wfProc.getPulseBaseline(iPulse);
	  pulses_values[iPulse][6]=wfProc.getPulseAbsAmplitude(iPulse);
	  pulses_values[iPulse][7]=wfProc.getPulseAmplitude(iPulse);
	  pulses_values[iPulse][8]=wfProc.getPulseTime(iPulse);
	  pulses_values[iPulse][9]=wfProc.getPulseCharge(iPulse);
	  pulses_values[iPulse][10]=wfProc.getPulseWidth(iPulse);
	  pulses_values[iPulse][11]=wfProc.getSPTemplateChi2(iPulse);
	  pulses_values[iPulse][12]=wfProc.getFitAmplitude(iPulse);
	  pulses_values[iPulse][13]=wfProc.getFitTime(iPulse);
	  pulses_values[iPulse][14]=wfProc.getFitRiseTime(iPulse);
	  pulses_values[iPulse][15]=wfProc.getFitFallTime(iPulse);
	  pulses_values[iPulse][16]=wfProc.getFitTime2Frac(iPulse);
	  pulses_values[iPulse][17]=wfProc.getFitFallTime2(iPulse);
	  pulses_values[iPulse][18]=wfProc.getFitBaseline(iPulse);
	  pulses_values[iPulse][19]=wfProc.getChi2(iPulse);
	  pulses_values[iPulse][20]=skipFit? -1 : wfProc.getNDF(iPulse);
	  pulses_values[iPulse][21]=wfProc.getChi2Refit(iPulse);
	  pulses_values[iPulse][22]=wfProc.getNDFRefit(iPulse);
        }
	
        
        
      }
        
      if(iEvent==0) wfProc.getBaselineHisto()->Write();
        
      if(floor(iEvent*200./nEvent)==(iEvent*200./nEvent)){
	std::cout << iEvent << "/" << nEvent << std::endl; //"\r";
	//std::cout.flush();
      }
        
        //  Idee alla base del codice:
        //
        //  1)Cerco il primo pulse: Il primo pulse che supera la soglia è quello che ha dato origine alla tua acquisizione
        //    abbastanza facile da fare basta controllare le ampiezze dei pulses del fit
        //
        //  2)Il pulse è 1p.e event: Prima di fare la misura della distribuzione dei tempi faccio un istrogramma delle ampiezze (o area) in
        //    questo modo riesci a capire qual'è l'ampiezza che corrisponde ad 1 p.e e cerco di avere un valore di ampiezza e una sigma
        //    per quanto riguarda il 1 pe event. Questa condizione la uso per selezionare solo i pulse buoni i.e. pulse che stanno
        //    in questo range come primo pulse
        //
        //  3)Cerco il pulse immediatamente successivo al primo (che sia 1 p.e./ 2.p.e non mi interessa lo cerco e guardo il tempo a cui
        //   è avvenuto salvandolo
        
        
        
        //Ora l'array con le informazioni sui Pulses è pieno
        //Se il numero di pulse nella Waveform è inferiore a 2 i.e. 0 o 1 non fare un
        //cazzo
    
      if(number_pulse>1){
        
        int k=0;
        pulse_found=0;
        
        //Cerco il primo pulse che supera la soglia facendo un for sui pulse di quella waveform
        for(k=0;k<number_pulse;k++){
	  
	  //std::cout<<"Pulse number :"<<k<<std::endl;
          
	  //Usare pulseabsmap oppure ampiezza del fit, per il momento pulseamp
	  //Salvo ampiezza del fit in MODULO
        //    if(pulses_values[k][12]<0){
        //        pulse_fit_amp=-pulses_values[k][12];
        //    }else{
        //
        //     pulse_fit_amp=pulses_values[k][12];
        // }
        
	  //std::cout<<"pulse_fit_amp :"<<pulses_values[k][6]<<std::endl;
        //In entrambi i Casi sottraggo al valore della baseline il valore dell'ampiezza e questo è il valore che poi vado a confrontare con
        //la threshold
        
        //minimum=pulses_values[k][18]-pulse_fit_amp;
        //std::cout<<"minimum :"<<minimum<<std::endl;
        
            //Negative signal, first found means first that has produced the trigger
	  if((pulses_values[k][6]<threshold)){
             
                //Controllo ampiezza del pulse, deve sare tra T e T-sigma
	    if( (pulses_values[k][6]>=(threshold-sigma_1pe_pulse))){
                  
                //Pulse 1 p.e found !
	      
	      pulse_found=1;
                    
	      break;
		    
                    
	    }
                
	  }
            
            
        }
        
        
        if((k<(number_pulse-1))&& pulse_found==1){
	  
	  //Ho un altro pulse che segue quello in analisi
	  //Cerco il pulse successivo e guardo il tempo
          
	  time=pulses_values[k+1][8]-pulses_values[k][8];
	  //std::cout<<"Event: "<<iEvent<<" Time difference between first two pulses: "<<time<<std::endl;
	  
            //Fill ntupla with time information
	  ntptime->Fill(iEvent,time);
          
        }
	
        
        
        
      }
    } 
    outFile->cd();
    ntpE->Write();
    ntp->Write();
    ntptime->Write();
    //std:: cout << "LAST FILE = " << gROOT->GetListOfFiles()->Last()->GetName() << std::endl;
    
    outFile->Close();
    procDir->cd();
    //delete outFile;

   
    std::cout<<"Number Waveform Discarded : "<<number_waveform_discarded<<std::endl;

    
    //if(gROOT->FindObjectAny(gROOT->GetListOfFiles()->Last()->GetName())){
    //std::cout<<"Closing ..."<<std::endl;

    //outFile->Close();
    //delete outFile;


    // }else{
    //std::cout<<"Object not found!"<<std::endl;

    //}

   
}

