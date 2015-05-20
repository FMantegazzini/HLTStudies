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

void DrawSpectra (TFile* file, Double_t n_events, std::string EBEE, std::string h_name, std::string title, std::string xTitle, Double_t xMax, Double_t yMax);

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
    
    TFile* f = new TFile(inputFiles.at(ii).c_str());

    TH1F* h_events = (TH1F*)f->Get("makeAlCaPhiSymSpectra/h_nEvents");
    Double_t num_events = h_events->GetBinContent(h_events->FindBin(0));
    cout << "Total number of events = " << num_events << endl;

    //EBM
    DrawSpectra(f,num_events,std::string("EBM"),std::string("eSpectrum"),std::string("EBM_eSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("E"),10,1);
    if (PU == "40" && bx == "50")  
      DrawSpectra(f,num_events,std::string("EBM"),std::string("etSpectrum"),std::string("EBM_etSpectrum_PU"+PU+"_"+bx+"ns_"+ multifit),std::string("E_{T}"),10,1);  
    else  
      DrawSpectra(f,num_events,std::string("EBM"),std::string("etSpectrum"),std::string("EBM_etSpectrum_PU"+PU+"_"+bx+"ns_"+ multifit),std::string("E_{T}"),10,1); 
    DrawSpectra(f,num_events,std::string("EBM"),std::string("calibration"),std::string("EBM_calibSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("Calibration"),0.5,1);

    //EBP
    DrawSpectra(f,num_events,std::string("EBP"),std::string("eSpectrum"),std::string("EBP_eSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("E"),10,1);
    if (PU == "40" && bx == "50")    
      DrawSpectra(f,num_events,std::string("EBP"),std::string("etSpectrum"),std::string("EBP_etSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("E_{T}"),10,1); 
    else   
      DrawSpectra(f,num_events,std::string("EBP"),std::string("etSpectrum"),std::string("EBP_etSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("E_{T}"),10,1); 
    DrawSpectra(f,num_events,std::string("EBP"),std::string("calibration"),std::string("EBP_calibSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("Calibration"),0.5,1);

    //EEM
    if ((PU == "20" && bx == "50" && multifit == "no_multifit") || (PU == "40" && bx == "50")) {
      DrawSpectra(f,num_events,std::string("EEM"),std::string("eSpectrum"),std::string("EEM_eSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("E"),30,1);    
      DrawSpectra(f,num_events,std::string("EEM"),std::string("etSpectrum"),std::string("EEM_etSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("E_{T}"),30,1);    
    }
    else  {
      DrawSpectra(f,num_events,std::string("EEM"),std::string("eSpectrum"),std::string("EEM_eSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("E"),30,1);    
      DrawSpectra(f,num_events,std::string("EEM"),std::string("etSpectrum"),std::string("EEM_etSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("E_{T}"),30,1);    
    }
    DrawSpectra(f,num_events,std::string("EEM"),std::string("calibration"),std::string("EEM_calibSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("Calibration"),0.5,1);

    //EEP
    if (PU == "40" && bx == "50" && multifit == "no_multifit") {
      DrawSpectra(f,num_events,std::string("EEP"),std::string("eSpectrum"),std::string("EEP_eSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("E"),30,1);    
      DrawSpectra(f,num_events,std::string("EEP"),std::string("etSpectrum"),std::string("EEP_etSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("E_{T}"),30,1);    
    }
    else  {
      DrawSpectra(f,num_events,std::string("EEP"),std::string("eSpectrum"),std::string("EEP_eSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("E"),30,1);    
      DrawSpectra(f,num_events,std::string("EEP"),std::string("etSpectrum"),std::string("EEP_etSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("E_{T}"),30,1);    
    }
    DrawSpectra(f,num_events,std::string("EEP"),std::string("calibration"),std::string("EEP_calibSpectrum_PU"+PU+"_"+bx+"ns_"+multifit),std::string("Calibration"),0.5,1);

    //Singles crystals energy spectra (EB ring 15, EEM ring 16)
    DrawSpectra(f,num_events,std::string("EB"),std::string("ieta_15_iphi"),std::string("EB_ring15_PU"+PU+"_"+bx+"ns_"+multifit),std::string("Energy"),10,1);
    DrawSpectra(f,num_events,std::string("EB"),std::string("iring_16"),std::string("EEM_ring16_PU"+PU+"_"+bx+"ns_"+multifit),std::string("Energy"),10,1);

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

//function to get the histos from the file and draw them in the same canvas
void DrawSpectra (TFile* file, Double_t n_events, std::string EBEE, std::string h_name, std::string title, std::string xTitle, Double_t xMax, Double_t yMax) {
  
  ostringstream t;

  int n;
  if (h_name == "ieta_15_iphi")
    n = 361;
  else if (h_name == "iring_16")
    n = 221;
  else if (h_name == "eSpectrum" || h_name == "etSpectrum" || h_name == "calibration") {
    if (EBEE == "EBM" || EBEE == "EBP") 
      n = 85;     
    if (EBEE == "EEM" || EBEE == "EEP")
      n = 39;
  }
    
  TCanvas* c = new TCanvas("c","c");

  for (int i=0; i<n; i++) {
    if (h_name == "ieta_15_iphi")  
      t << "makeAlCaPhiSymSpectra/" << EBEE << "_" << h_name << "_" << i;
    else
      t << "makeAlCaPhiSymSpectra/" << EBEE << "_" << h_name << "_" << i+1;
    TH1F* histo = (TH1F*)file->Get(t.str().c_str()); 
    t.str("");

    if (i==0) {
      histo -> SetTitle(title.c_str());
      histo -> GetXaxis() -> SetTitle(xTitle.c_str());
      histo -> GetXaxis() -> SetLabelSize(0.04);
      histo -> GetXaxis() -> SetTitleSize(0.04);
      histo -> GetXaxis() -> SetTitleOffset(1.);
      histo -> Scale(1./n_events);
      Int_t rebin = (histo -> GetNbinsX())/1000;
      histo -> Rebin (rebin);
      histo -> SetAxisRange(0.,xMax, "X");
      histo -> SetAxisRange(0.,yMax, "Y");
      histo -> SetStats(0);
      histo -> SetLineColor(kBlue);
      histo -> SetLineWidth(1.);
      c->cd();
      histo -> Draw();
    }

    else if (i > 0) {
      histo -> Scale(1./n_events);
      Int_t rebin = (histo -> GetNbinsX())/1000;
      histo -> Rebin (rebin); 
      histo -> SetAxisRange(0.,xMax, "X");
      histo -> SetAxisRange(0.,yMax, "Y");
      histo -> SetLineColor(kBlue);
      histo -> SetLineWidth(1.);
      c->cd();
      histo -> Draw("same");
    }
 
  }

  c -> Print((title + ".png").c_str(), "png");
  c -> Print((title + ".pdf").c_str(), "pdf"); 
	
  delete c;

}
