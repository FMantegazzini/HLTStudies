/*
C++ code to draw occupancy plots.
This program read spectra_output_total.root, that is the union of the .root files created by plugins/MakeAlCaPhiSymSpectra.cc

To compile: c++ -o finalSpectra `root-config --cflags --glibs` finalSpectra.cpp
To run: ./finalSpectra.cpp scenariosList.txt
*/

using namespace std;

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <map>
#include <vector>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TTree.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TPaveStats.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLine.h"

string uintToString (unsigned int);

string Double_tToString (Double_t);

int main (int argc, char** argv) {

  std::string PU;
  std::string bx;
  std::string yesno;
  std::string multifit;
  
  char* inputLIST = argv[1];
  std::string inputList = std::string(inputLIST);
  
  std::cout << "inputList = " << inputList << std::endl;
  std::vector<std::string> inputFiles;

  char dumps[500];
  FILE *f_dumps;

  f_dumps = fopen(inputList.c_str(),"r"); //read list of files .root
  while(fscanf(f_dumps,"%s \n", dumps) !=EOF ){
    //std::cout << "Reading file from the list: " << dumps << std::endl;
    std::string DUMPS = std::string(dumps);
    if(DUMPS.find("#") != std::string::npos) continue;
    inputFiles.push_back(DUMPS);
  }  

  for (unsigned int ii = 0; ii < inputFiles.size(); ii++) { //loop over the files list
    
    std::cout << std::endl << ">>>> Reading inputFile: " << inputFiles.at(ii) << std::endl;

    /*
    cout << "PU?" << endl;
    cin >> PU;
    cout << "bx?" << endl;
    cin >> bx;
    cout << "multifit? (yes/no)" << endl;
    cin >> yesno;
    if (yesno == "yes")
      multifit = "multifit";
    else if (yesno == "no")
      multifit = "no_multifit";
    */

    //for the default scenarioList.txt use this:
    if (ii==0) {
      PU = "20";
      bx = "25";
      multifit = "multifit";
      cout << "PU = " << PU << "; bx = " << bx << "; " << multifit << endl;
    }
    if (ii==1) {
      PU = "20";
      bx = "25";
      multifit = "no_multifit";
      cout << "PU = " << PU << "; bx = " << bx << "; " << multifit << endl;
    }
    if (ii==2) {
      PU = "20";
      bx = "50";
      multifit = "multifit";
      cout << "PU = " << PU << "; bx = " << bx << "; " << multifit << endl;
    }
    if (ii==3) {
      PU = "20";
      bx = "50";
      multifit = "no_multifit";
      cout << "PU = " << PU << "; bx = " << bx << "; " << multifit << endl;
    }
    if (ii==4) {
      PU = "40";
      bx = "25";
      multifit = "multifit";
      cout << "PU = " << PU << "; bx = " << bx << "; " << multifit << endl;
    }
    if (ii==5) {
      PU = "40";
      bx = "25";
      multifit = "no_multifit";
      cout << "PU = " << PU << "; bx = " << bx << "; " << multifit << endl;
    }
    if (ii==6) {
      PU = "40";
      bx = "50";
      multifit = "multifit";
      cout << "PU = " << PU << "; bx = " << bx << "; " << multifit << endl;
    }
    if (ii==7) {
      PU = "40";
      bx = "50";
      multifit = "no_multifit";
      cout << "PU = " << PU << "; bx = " << bx << "; " << multifit << endl;
    }
    
    //std::string str = uintToString(ii+1);
    //TFile *outputFile = new TFile (("spectra_" + str + ".root").c_str(),"RECREATE"); 
    //std::cout << "New file created: " << "spectra_" << str << ".root" << std::endl; 

    static const int EB_rings = 85;
    static const int EE_rings = 39;

    TFile* f = new TFile(inputFiles.at(ii).c_str());
    f->cd("makeAlCaPhiSymSpectra");
  
    TH1F* h_events = (TH1F*)f->Get("h_nEvents");
    
    ostringstream t;

    //---EB spectra
    TCanvas* e_EBM = new TCanvas("e_EBM","e_EBM");
    TCanvas* e_EBP = new TCanvas("e_EBP","e_EBP");
  
    TCanvas* et_EBM = new TCanvas("et_EBM","et_EBM");
    TCanvas* et_EBP = new TCanvas("et_EBP","et_EBP");

    TCanvas* cal_EBM = new TCanvas("cal_EBM","cal_EBM");
    TCanvas* cal_EBP = new TCanvas("cal_EBP","cal_EBP");
     
    for (int i=0; i<EB_rings; i++) {
     
      t << "EBM_eSpectrum_" << i+1;
      TH1F* h_eSpectrum_EBM = (TH1F*)f->Get(t.str().c_str()); 
      t.str("");     
      t << "EBP_eSpectrum_" << i+1;
      TH1F* h_eSpectrum_EBP = (TH1F*)f->Get(t.str().c_str()); 
      t.str("");

      t << "EBM_etSpectrum_" << i+1;
      TH1F* h_etSpectrum_EBM = (TH1F*)f->Get(t.str().c_str()); 
      t.str("");
      t << "EBP_etSpectrum_" << i+1;
      TH1F* h_etSpectrum_EBP = (TH1F*)f->Get(t.str().c_str()); 
      t.str("");     

      t << "EBM_calibration_" << i+1;
      TH1F* h_calSpectrum_EBM = (TH1F*)f->Get(t.str().c_str()); 
      t.str("");
      t << "EBP_calibration_" << i+1;
      TH1F* h_calSpectrum_EBP = (TH1F*)f->Get(t.str().c_str()); 
      t.str("");     

      e_EBM -> cd();
      h_eSpectrum_EBM -> Draw();
      e_EBP -> cd();
      h_eSpectrum_EBP -> Draw();
     
      et_EBM -> cd();
      h_etSpectrum_EBM -> Draw();
      et_EBP -> cd();
      h_etSpectrum_EBP -> Draw();
      
      cal_EBM -> cd();
      h_calSpectrum_EBM -> Draw();
      cal_EBP -> cd();
      h_calSpectrum_EBP -> Draw();
     
    }
    
    e_EBM -> Print("EBM_eSpectra.png", "png");
    e_EBM -> Print("EBM_eSpectra.pdf", "pdf"); 
    et_EBP -> Print("EBP_etSpectra.png", "png");
    et_EBP -> Print("EBP_etSpectra.pdf", "pdf"); 
    cal_EBP -> Print("EBP_calibSpectra.png", "png");
    cal_EBP -> Print("EBP_calibSpectra.pdf", "pdf"); 
      
    delete e_EBM;
    delete e_EBP;
    delete et_EBM;
    delete et_EBP;
    delete cal_EBM;
    delete cal_EBP;

    //---EE spectra
    TCanvas* e_EEM = new TCanvas("e_EEM","e_EEM");
    TCanvas* e_EEP = new TCanvas("e_EEP","e_EEP");
  
    TCanvas* et_EEM = new TCanvas("et_EEM","et_EEM");
    TCanvas* et_EEP = new TCanvas("et_EEP","et_EEP");

    TCanvas* cal_EEM = new TCanvas("cal_EEM","cal_EEM");
    TCanvas* cal_EEP = new TCanvas("cal_EEP","cal_EEP");
     
    for (int i=0; i<EE_rings; i++) {
     
      t << "EEM_eSpectrum_" << i+1;
      TH1F* h_eSpectrum_EEM = (TH1F*)f->Get(t.str().c_str()); 
      t.str("");     
      t << "EEP_eSpectrum_" << i+1;
      TH1F* h_eSpectrum_EEP = (TH1F*)f->Get(t.str().c_str()); 
      t.str("");

      t << "EEM_etSpectrum_" << i+1;
      TH1F* h_etSpectrum_EEM = (TH1F*)f->Get(t.str().c_str()); 
      t.str("");
      t << "EEP_etSpectrum_" << i+1;
      TH1F* h_etSpectrum_EEP = (TH1F*)f->Get(t.str().c_str()); 
      t.str("");     

      t << "EEM_calibration_" << i+1;
      TH1F* h_calSpectrum_EEM = (TH1F*)f->Get(t.str().c_str()); 
      t.str("");
      t << "EEP_calibration_" << i+1;
      TH1F* h_calSpectrum_EEP = (TH1F*)f->Get(t.str().c_str()); 
      t.str("");     

      e_EEM -> cd();
      h_eSpectrum_EEM -> Draw();
      e_EEP -> cd();
      h_eSpectrum_EEP -> Draw();
     
      et_EEM -> cd();
      h_etSpectrum_EEM -> Draw();
      et_EEP -> cd();
      h_etSpectrum_EEP -> Draw();
      
      cal_EEM -> cd();
      h_calSpectrum_EEM -> Draw();
      cal_EEP -> cd();
      h_calSpectrum_EEP -> Draw();
     
    }
   
    e_EEM -> Print("EEM_eSpectra.png", "png");
    e_EEM -> Print("EEM_eSpectra.pdf", "pdf"); 
    et_EEP -> Print("EEP_etSpectra.png", "png");
    et_EEP -> Print("EEP_etSpectra.pdf", "pdf"); 
    cal_EEP -> Print("EEP_calibSpectra.png", "png");
    cal_EEP -> Print("EEP_calibSpectra.pdf", "pdf"); 
      
    delete e_EEM;
    delete e_EEP;
    delete et_EEM;
    delete et_EEP;
    delete cal_EEM;
    delete cal_EEP;

  } //files list
  
} //main


// function to convert unsigned int into string
string uintToString(unsigned int val) {
  char buff[500];
  sprintf(buff, "%u", val);
  string str = buff;
  return(str);
}

// function to convert Double_t into string
string Double_tToString (Double_t number) {
  std::ostringstream buff;
  buff<<number;
  return buff.str();   
}
