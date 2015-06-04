/*
C++ code to draw occupancy plots.
This program read outputTotal.root, that is the union of the .root files created by bin/MakeAlCaPhiSymOccupancyPlot.cpp

To compile: c++ -o finalPlots `root-config --cflags --glibs` finalPlots.cpp TEndcapRings.cc
To run: ./finalPlots scenariosList.txt

OUTPUT:
- TH2F occupancy plots
- computation of the necessary occupancy to calibrate in 10 hours assuming: rate = 1.5 kHz, statistical precision = 0.1% (---> 1000 hits)
- computation of mean hits per event and of lower occupancy for EB, EEM, EEP
- TH2F time vs ring distributions
- TGraph mean occupancy vs ring
- TGraph lowest occupancy vs ring
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

#include "TEndcapRings.h"

void drawHistos(TH2F* h2, std::string Title, std::string xTitle, std::string yTitle, Double_t zmax, Double_t ymax);

void drawGraphs(TGraph* g1,TGraph* g2, std::string Title, std::string g1_Title, std::string g2_Title, float xmin, float xmax, float ymin, float ymax, Double_t lineValue_1, Double_t lineValue_2, std::string s, std::string time);

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
    

    std::string str = uintToString(ii+1);
    TFile *outputFile = new TFile (("occupancy_rings_" + str + ".root").c_str(),"RECREATE"); 
    std::cout << "New file created: " << "occupancy_rings_" << str << ".root" << std::endl; 

    TFile* f = new TFile(inputFiles.at(ii).c_str());
  
    //if you are reading .root with subdir makeAlCaPhiSymSpectra
    TH1F* h_events = (TH1F*)f->Get("makeAlCaPhiSymSpectra/h_nEvents");
    TH2F* h2_EB = (TH2F*)f->Get("makeAlCaPhiSymSpectra/h2_hitOccupancy_EB_highCut1"); //360 x 170 bins
    TH2F* h2_EEM = (TH2F*)f->Get("makeAlCaPhiSymSpectra/h2_hitOccupancy_EEM_highCut1"); //100 x 100 bins
    TH2F* h2_EEP = (TH2F*)f->Get("makeAlCaPhiSymSpectra/h2_hitOccupancy_EEP_highCut1"); //100 x 100 bins*/

    //if you are reading .root without subdir
    /*TH1F* h_events = (TH1F*)f->Get("h_nEvents");
    TH2F* h2_EB = (TH2F*)f->Get("h2_hitOccupancy_EB"); //360 x 170 bins
    TH2F* h2_EEM = (TH2F*)f->Get("h2_hitOccupancy_EEM"); //100 x 100 bins
    TH2F* h2_EEP = (TH2F*)f->Get("h2_hitOccupancy_EEP"); //100 x 100 bins */
   
    Double_t num_events = h_events->GetBinContent(h_events->FindBin(0));
    cout << "Total number of events = " << num_events << endl;

    h2_EB->Scale(1./num_events);
    h2_EEM->Scale(1./num_events);
    h2_EEP->Scale(1./num_events);
    cout << "Scaling for " << 1./num_events << endl;
        
    //--- draw scaled occupancy graphs
   
    drawHistos(h2_EB, std::string(("EB_PU" + PU + "_" + bx + "ns_" + multifit).c_str()), "iphi", "ieta", 0.001, 0);
    drawHistos(h2_EEM, std::string(("EEM_PU" + PU + "_" + bx + "ns_" + multifit).c_str()), "ix", "iy", 0.001, 0);
    drawHistos(h2_EEP, std::string(("EEP_PU" + PU + "_" + bx + "ns_" + multifit).c_str()), "ix", "iy", 0.001, 0);

    //--- computation of the necessary occupancy to calibrate in 10 hours assuming:
    //  - rate = 1.5 kHz
    //  - statistical precision = 0.1%
    //  - necessary number of hits to reach this precision = 1000

    Double_t hits = 1000;
    Double_t rate = 1.5;
    Double_t time_hours = 10;

    Double_t necessary_occupancy = hits/(time_hours * rate * 1000 * 3600);
    cout << "Necessary occupancy to calibrate in 10 hours = " << necessary_occupancy << endl;

    //---computation of mean hits per event and of lowest occupancy for EB, EEM, EEP
    //---fill TH2 histos with time to calibrate distributions for every ring
 
    Double_t EB_mean = 0.;
    Double_t EEM_mean= 0.;
    Double_t EEP_mean = 0.;
    Double_t EB_lowest = 10.;
    Double_t EEM_lowest = 10.;
    Double_t EEP_lowest = 10.;
    
    int EB_channels = 0;
    for (int i=1; i <= h2_EB->GetNbinsX(); i++) {
      for (int j=1; j <= h2_EB->GetNbinsY(); j++) {
	if (h2_EB->GetBinContent(h2_EB->FindBin(i,j)) > 0) { //find mean occupancy
	  EB_channels = EB_channels + 1;
	  EB_mean = EB_mean + h2_EB->GetBinContent(h2_EB->FindBin(i,j));
	}
      
	if (h2_EB->GetBinContent(h2_EB->FindBin(i,j)) < EB_lowest && h2_EB->GetBinContent(h2_EB->FindBin(i,j)) != 0) {  //Find lowest occupancy
	  EB_lowest = h2_EB->GetBinContent(h2_EB->FindBin(i,j));
	}
      }
    }
    EB_mean = EB_mean/EB_channels;

    int EEM_channels = 0;
    for (int i=1; i <= h2_EEM->GetNbinsX(); i++) {
      for (int j=1; j <= h2_EEM->GetNbinsY(); j++) {
	if (h2_EEM->GetBinContent(h2_EEM->FindBin(i,j)) > 0) { //find lowest occupancy
	  EEM_channels = EEM_channels + 1;
	  EEM_mean = EEM_mean + h2_EEM->GetBinContent(h2_EEM->FindBin(i,j));
	}

	if (h2_EEM->GetBinContent(h2_EEM->FindBin(i,j)) < EEM_lowest && h2_EEM->GetBinContent(h2_EEM->FindBin(i,j)) != 0) //find lowest occupancy
	  EEM_lowest = h2_EEM->GetBinContent(h2_EEM->FindBin(i,j));
      }
    }
    EEM_mean = EEM_mean/EEM_channels;

    int EEP_channels = 0;
    for (int i=1; i <= h2_EEP->GetNbinsX(); i++) {
      for (int j=1; j <= h2_EEP->GetNbinsY(); j++) {
	if (h2_EEP->GetBinContent(h2_EEP->FindBin(i,j)) > 0) { //find mean occupancy
	  EEP_channels = EEP_channels + 1;
	  EEP_mean = EEP_mean + h2_EEP->GetBinContent(h2_EEP->FindBin(i,j));
	}

	if (h2_EEP->GetBinContent(h2_EEP->FindBin(i,j)) < EEP_lowest && h2_EEP->GetBinContent(h2_EEP->FindBin(i,j)) != 0) //find lowest occupancy
	  EEP_lowest = h2_EEP->GetBinContent(h2_EEP->FindBin(i,j));
      }
    }
    EEP_mean = EEP_mean/EEP_channels;

    Double_t EE_mean = (EEM_mean + EEP_mean)/2.;
    Double_t EE_lowest = 0.;
    if (EEM_lowest < EEP_lowest)
      EE_lowest = EEM_lowest;
    else if (EEM_lowest >= EEP_lowest)
      EE_lowest = EEP_lowest;

    cout << "EB mean occupancy = " << EB_mean << endl;
    cout << "EEM mean occupancy = " << EEM_mean << endl;
    cout << "EEP mean occupancy = " << EEP_mean << endl;
    cout << "EE mean occupancy = " << EE_mean << endl;
  
    cout << "EB lowest occupancy = " << EB_lowest  << endl;
    cout << "EEM lowest occupancy = " << EEM_lowest << endl;
    cout << "EEP lowest occupancy = " << EEP_lowest << endl;
    cout << "EE lowest occupancy = " << EE_lowest << endl;

    //--- computation of the necessary time to calibrate with EB(EE)_mean and EB(EE)_lowest occupancy assuming;
    //  - rate = 1.5 kHz
    //  - statistical precision = 0.1%
    //  - necessary number of hits to reach this precision = 1000

    Double_t EB_mean_time = hits/(EB_mean * rate * 1000 * 3600);
    Double_t EE_mean_time = hits/(EE_mean * rate * 1000 * 3600);
    Double_t EB_lowest_time = hits/(EB_lowest * rate * 1000 * 3600);
    Double_t EE_lowest_time = hits/(EE_lowest * rate * 1000 * 3600);

    cout << "Necessary time to calibrate with EB mean occupancy = " << EB_mean_time << " hours" << endl;
    cout << "Necessary time to calibrate with EE mean occupancy = " << EE_mean_time << " hours" << endl;
    cout << "Necessary time to calibrate with EB lowest occupancy = " << EB_lowest_time << " hours" << endl;
    cout << "Necessary time to calibrate with EE lowest occupancy = " << EE_lowest_time << " hours" << endl;
  
    //--- create and fill histos with occupancy for every ring
    static const int EB_rings = 85;
    static const int EE_rings = 39;

    TEndcapRings *eRings = new TEndcapRings(); 

    Int_t nBins = 10000;   
    Double_t occupancyMin = 0.;
    Double_t occupancyMax = 0.02;
    Double_t meanOccupancy = 0.;
    Double_t lowestOccupancy = 0.;

    std::vector<TH1F*> EBM_occupancy_histos;
    std::vector<TH1F*> EBP_occupancy_histos;
    std::vector<TH1F*> EEM_occupancy_histos;
    std::vector<TH1F*> EEP_occupancy_histos;
    
    EBM_occupancy_histos.resize(EB_rings);
    EBP_occupancy_histos.resize(EB_rings);
    EEM_occupancy_histos.resize(EE_rings);
    EEP_occupancy_histos.resize(EE_rings);
  
    ostringstream t;

    for (int i=0; i<EB_rings; i++) { //EB-
      t << "EBM_occupancy_" << i+1;
      EBM_occupancy_histos[i]=new TH1F(t.str().c_str(),";occupancy",nBins,occupancyMin,occupancyMax); 
      t.str("");
    }

    for (int i=0; i<EB_rings; i++) { //EB+
      t << "EBP_occupancy_" << i+1;
      EBP_occupancy_histos[i]=new TH1F(t.str().c_str(),";occupancy",nBins,occupancyMin,occupancyMax); 
      t.str("");
    }

    for (int i=0; i<EE_rings; i++) { //EE-
      t << "EEM_occupancy_" << i+1;
      EEM_occupancy_histos[i]=new TH1F(t.str().c_str(),";occupancy",nBins,occupancyMin,occupancyMax); 
      t.str("");
    }

    for (int i=0; i<EE_rings; i++) { //EE+
      t << "EEP_occupancy_" << i+1;
      EEP_occupancy_histos[i]=new TH1F(t.str().c_str(),";occupancy",nBins,occupancyMin,occupancyMax); 
      t.str("");
    }

    //--- create histos for with time to calibrate distributions for every ring    
    TH2F* h2_EB_time_rings = new TH2F("h2_EB_time_rings","h2_EB_time_rings",170,-85,85,100,0.,2.);
    TH2F* h2_EEM_time_rings = new TH2F("h2_EEM_time_rings","h2_EEM_time_rings",39,1,40,200,0.,30.);
    TH2F* h2_EEP_time_rings = new TH2F("h2_EEP_time_rings","h2_EEP_time_rings",39,1,40,200,0.,30.);

    Double_t EB_time = 0.;
    Double_t EEM_time = 0.;
    Double_t EEP_time = 0.;

    //---create time maps
    TH2F* h2_EB_timeMap = new TH2F("h2_EB_timeMap","h2_EB_timeMap",360,0,360,170,-85,85);
    TH2F* h2_EEM_timeMap = new TH2F("h2_EEM_timeMap","h2_EEM_timeMap",100,1,101,100,1,101);
    TH2F* h2_EEP_timeMap = new TH2F("h2_EEP_timeMap","h2_EEP_timeMap",100,1,101,100,1,101);
  
    //fill rings histos, TH2F time vs ring, TH2F time maps
    for (int iphi=0; iphi<361; iphi++) { //EB-
      for (int ieta=-85; ieta<0; ieta++) {
	if ( h2_EB->GetBinContent(h2_EB->FindBin(iphi,ieta)) != 0 ) {
	  EBM_occupancy_histos[ieta+85]->Fill( h2_EB->GetBinContent(h2_EB->FindBin(iphi,ieta)) ); 
	  EB_time = hits/( h2_EB->GetBinContent(h2_EB->FindBin(iphi,ieta) ) * rate * 1000 * 3600 );
	  h2_EB_time_rings->Fill(ieta,EB_time);
	  h2_EB_timeMap->SetBinContent(h2_EB_timeMap->FindBin(iphi,ieta),EB_time);
	  if (EB_time > 10)
	    std::cout << "---EBM iphi = " << iphi << ", ieta = " << ieta << " ---> time = " << EB_time << std::endl;
	}
      }
    }
  
    for (int iphi=0; iphi<361; iphi++) { //EB+
      for (int ieta=1; ieta<86; ieta++) {
	if ( h2_EB->GetBinContent(h2_EB->FindBin(iphi,ieta)) != 0 ) {
	  EBP_occupancy_histos[ieta-1]->Fill( h2_EB->GetBinContent(h2_EB->FindBin(iphi,ieta)) ); 
	  EB_time = hits/( h2_EB->GetBinContent(h2_EB->FindBin(iphi,ieta) ) * rate * 1000 * 3600 );
	  h2_EB_time_rings->Fill(ieta,EB_time);
	  h2_EB_timeMap->SetBinContent(h2_EB_timeMap->FindBin(iphi,ieta),EB_time);
	  if (EB_time > 10)
	    std::cout << "---EBP iphi = " << iphi << ", ieta = " << ieta << " ---> time = " << EB_time << std::endl;
	}
      }
    }

    for (int ix=1; ix<101; ix++) { //EE-
      for (int iy=1; iy<101; iy++) {
	int iring = eRings->GetEndcapRing(ix,iy,-1);
	if ( h2_EEM->GetBinContent(h2_EEM->FindBin(ix,iy)) != 0 ) {
	  EEM_occupancy_histos[iring]->Fill( h2_EEM->GetBinContent(h2_EEM->FindBin(ix,iy)) );
	  EEM_time = hits/( h2_EEM->GetBinContent(h2_EEM->FindBin(ix,iy)) * rate * 1000 * 3600 );
	  h2_EEM_time_rings->Fill(iring,EEM_time);
	  h2_EEM_timeMap->SetBinContent(h2_EEM_timeMap->FindBin(ix,iy),EEM_time);
	  if (EEM_time > 10)
	    std::cout << "---EEM ix = " << ix << ", iy = " << iy << ", iring = " << iring << " ---> time = " << EEM_time << std::endl;
	} 
      }
    }

    for (int ix=1; ix<101; ix++) { //EE+
      for (int iy=1; iy<101; iy++) {
	int iring = eRings->GetEndcapRing(ix,iy,-1);
	if ( h2_EEP->GetBinContent(h2_EEP->FindBin(ix,iy)) != 0 ) {
	  EEP_occupancy_histos[iring]->Fill( h2_EEP->GetBinContent(h2_EEP->FindBin(ix,iy)) ); 
	  EEP_time = hits/( h2_EEP->GetBinContent(h2_EEP->FindBin(ix,iy)) * rate * 1000 * 3600 );
	  h2_EEP_time_rings->Fill(iring,EEP_time);
	  h2_EEP_timeMap->SetBinContent(h2_EEP_timeMap->FindBin(ix,iy),EEP_time);
	  if (EEP_time > 10)
	    std::cout << "---EEP ix = " << ix << ", iy = " << iy << ", iring = " << iring << " ---> time = " << EEP_time << std::endl;
	  if (EEP_time > 5 && EEP_time < 10)
	    std::cout << "---EEP ix = " << ix << ", iy = " << iy << ", iring = " << iring << " ---> time = " << EEP_time << std::endl;
	  if (EEP_time > 0.5 && EEP_time < 5 && iring == 1)
	    std::cout << "---EEP ix = " << ix << ", iy = " << iy << ", iring = " << iring << " ---> time = " << EEP_time << std::endl;
	}
      }
    }

    //---write ring occupancy histos
    outputFile->cd();

    for(int i=0;i<EB_rings;i++){
      EBM_occupancy_histos[i]->Write();
      EBP_occupancy_histos[i]->Write();
    }

    for(int i=0;i<EE_rings;i++){
      EEM_occupancy_histos[i]->Write();
      EEP_occupancy_histos[i]->Write();
    }

    //---draw time TH2F time distributions

    if (PU == "20") {
      drawHistos(h2_EB_time_rings,std::string(("EB_TimeDistribution_PU"+PU+"_"+bx+"ns_"+multifit).c_str()),"iRing","time to calibrate at 0.1% precision level (h)",140,0.8 );
      drawHistos(h2_EEM_time_rings,std::string(("EEM_TimeDistribution_PU"+PU+"_"+bx+"ns_"+multifit).c_str()),"iRing","time to calibrate at 0.1% precision level (h)",150,10.);
      drawHistos(h2_EEP_time_rings,std::string(("EEP_TimeDistribution_PU"+PU+"_"+bx+"ns_"+multifit).c_str()),"iRing","time to calibrate at 0.1% precision level (h)",150,10.);
    }

    if (PU == "40") {
      drawHistos(h2_EB_time_rings,std::string(("EB_TimeDistribution_PU"+PU+"_"+bx+"ns_"+multifit).c_str()),"iRing","time to calibrate at 0.1% precision level (h)",200,0.4);
      drawHistos(h2_EEM_time_rings,std::string(("EEM_TimeDistribution_PU"+PU+"_"+bx+"ns_"+multifit).c_str()),"iRing","time to calibrate at 0.1% precision level (h)",200,4.);
      drawHistos(h2_EEP_time_rings,std::string(("EEP_TimeDistribution_PU"+PU+"_"+bx+"ns_"+multifit).c_str()),"iRing","time to calibrate at 0.1% precision level (h)",200,4.);
    }

    //draw time maps
    drawHistos(h2_EB_timeMap,std::string(("EB_TimeMap_PU"+PU+"_"+bx+"ns_"+multifit).c_str()),"iphi","ieta",0.6,0);
    drawHistos(h2_EEM_timeMap,std::string(("EEM_TimeMap_PU"+PU+"_"+bx+"ns_"+multifit).c_str()),"ix","iy",4,0);
    drawHistos(h2_EEP_timeMap,std::string(("EEP_TimeMap_PU"+PU+"_"+bx+"ns_"+multifit).c_str()),"ix","iy",4,0);
   
    //graphs for mean and lowest occupancy

    TGraph *EBM_meanOccupancy = new TGraph();
    TGraph *EBP_meanOccupancy = new TGraph(); 
    TGraph *EEM_meanOccupancy = new TGraph(); 
    TGraph *EEP_meanOccupancy = new TGraph(); 

    TGraph *EBM_lowestOccupancy = new TGraph();
    TGraph *EBP_lowestOccupancy = new TGraph(); 
    TGraph *EEM_lowestOccupancy = new TGraph(); 
    TGraph *EEP_lowestOccupancy = new TGraph(); 
        
    for(int i = 0; i < EB_rings; i++) { //EB-
      meanOccupancy = EBM_occupancy_histos[i]->GetMean();      
      EBM_meanOccupancy->SetPoint(i,85-i,meanOccupancy);
      lowestOccupancy = (EBM_occupancy_histos[i]->FindFirstBinAbove(0.)) * (occupancyMax - occupancyMin) / nBins;
      EBM_lowestOccupancy->SetPoint(i,85-i,lowestOccupancy);
    }

    for(int i = 0; i < EB_rings; i++) { //EB+
      meanOccupancy = EBP_occupancy_histos[i]->GetMean();      
      EBP_meanOccupancy->SetPoint(i,i+1,meanOccupancy);
      lowestOccupancy = (EBP_occupancy_histos[i]->FindFirstBinAbove(0.)) * (occupancyMax - occupancyMin) / nBins;
      EBP_lowestOccupancy->SetPoint(i,i+1,lowestOccupancy);
    }

    for(int i = 0; i < EE_rings; i++) { //EE-
      meanOccupancy = EEM_occupancy_histos[i]->GetMean();      
      EEM_meanOccupancy->SetPoint(i,i+1,meanOccupancy);
      lowestOccupancy = (EEM_occupancy_histos[i]->FindFirstBinAbove(0.)) * (occupancyMax - occupancyMin) / nBins;
      EEM_lowestOccupancy->SetPoint(i,i+1,lowestOccupancy);
    }

    for(int i = 0; i < EE_rings; i++) { //EE+
      meanOccupancy = EEP_occupancy_histos[i]->GetMean();      
      EEP_meanOccupancy->SetPoint(i,i+1,meanOccupancy);
      lowestOccupancy = (EEP_occupancy_histos[i]->FindFirstBinAbove(0.)) * (occupancyMax - occupancyMin) / nBins;
      EEP_lowestOccupancy->SetPoint(i,i+1,lowestOccupancy);
    }

    outputFile->Close();

    //draw graphs
   
    if (PU == "20") {
      drawGraphs(EBM_meanOccupancy,EBP_meanOccupancy,std::string("EB_PU" + PU + "_" + bx + "ns_" + multifit),std::string("EBM"),std::string("EBP"),0,87,0,0.0014, necessary_occupancy, EB_mean, std::string("mean"), std::string(Double_tToString(EB_mean_time))); 

      drawGraphs(EEM_meanOccupancy,EEP_meanOccupancy,std::string("EE_PU" + PU + "_" + bx + "ns_" + multifit),std::string("EEM"),std::string("EEP"),0,40,0,0.0014, necessary_occupancy, EE_mean, std::string("mean"), std::string(Double_tToString(EE_mean_time)));

      drawGraphs(EBM_lowestOccupancy,EBP_lowestOccupancy,std::string("EB_PU" + PU + "_" + bx + "ns_" + multifit),std::string("EBM"),std::string("EBP"),0,87,0,0.0014, necessary_occupancy, EB_lowest, std::string("lowest"), std::string(Double_tToString(EB_lowest_time)));

      drawGraphs(EEM_lowestOccupancy,EEP_lowestOccupancy,std::string("EE_PU" + PU + "_" + bx + "ns_" + multifit),std::string("EEM"),std::string("EEP"),0,40,0,0.0014, necessary_occupancy,EE_lowest, std::string("lowest"), std::string(Double_tToString(EE_lowest_time)));
    }

    if (PU == "40") {
      drawGraphs(EBM_meanOccupancy,EBP_meanOccupancy,std::string("EB_PU" + PU + "_" + bx + "ns_" + multifit),std::string("EBM"),std::string("EBP"),0,87,0,0.004, necessary_occupancy, EB_mean, std::string("mean"), std::string(Double_tToString(EB_mean_time))); 

      drawGraphs(EEM_meanOccupancy,EEP_meanOccupancy,std::string("EE_PU" + PU + "_" + bx + "ns_" + multifit),std::string("EEM"),std::string("EEP"),0,40,0,0.004, necessary_occupancy, EE_mean, std::string("mean"), std::string(Double_tToString(EE_mean_time)));

      drawGraphs(EBM_lowestOccupancy,EBP_lowestOccupancy,std::string("EB_PU" + PU + "_" + bx + "ns_" + multifit),std::string("EBM"),std::string("EBP"),0,87,0,0.004, necessary_occupancy, EB_lowest, std::string("lowest"), std::string(Double_tToString(EB_lowest_time)));

      drawGraphs(EEM_lowestOccupancy,EEP_lowestOccupancy,std::string("EE_PU" + PU + "_" + bx + "ns_" + multifit),std::string("EEM"),std::string("EEP"),0,40,0,0.004, necessary_occupancy,EE_lowest, std::string("lowest"), std::string(Double_tToString(EE_lowest_time)));
    }

  } //files list
  
} //main

void drawHistos(TH2F* h2, std::string Title, std::string xTitle, std::string yTitle, Double_t zmax, Double_t ymax){

  h2 -> SetTitle(Title.c_str());
  h2 -> GetXaxis() -> SetTitle(xTitle.c_str());
  h2 -> GetYaxis() -> SetTitle(yTitle.c_str());
  //h2 -> GetZaxis() -> SetTitle("");
    
  h2 -> GetXaxis() -> SetLabelSize(0.04);
  h2 -> GetXaxis() -> SetTitleSize(0.04);
  h2 -> GetXaxis() -> SetTitleOffset(0.9);
  h2 -> GetYaxis() -> SetLabelSize(0.04);
  h2 -> GetYaxis() -> SetTitleSize(0.04);
  h2 -> GetYaxis() -> SetTitleOffset(1.2);
  h2 -> GetZaxis() -> SetLabelSize(0.04);
  h2 -> GetZaxis() -> SetTitleSize(0.04);
  h2 -> GetZaxis() -> SetTitleOffset(1.2);
  h2-> SetStats(0);

  if (zmax != 0) {
    //h2 -> GetZaxis() -> SetRangeUser(0,zmax);
    h2 -> SetAxisRange(0.,zmax, "Z");
  }
 
  if (ymax != 0) {
    //h2 -> GetYaxis() -> SetRangeUser(0,ymax);
    h2 -> SetAxisRange(0.,ymax, "Y");
  }

  TCanvas* c = new TCanvas("c","c");
  c->SetRightMargin(0.2);
  c -> cd();
  h2 -> Draw("colz");
  c -> Print((Title + ".png").c_str(), "png");
  c -> Print((Title + ".pdf").c_str(), "pdf"); 
  delete c;
}

void drawGraphs(TGraph* g1,TGraph* g2, std::string Title, std::string g1_Title, std::string g2_Title, float xmin, float xmax, float ymin, float ymax, Double_t lineValue_1, Double_t lineValue_2, std::string s, std::string time) {
  
  gStyle -> SetOptFit (00111);
  gStyle -> SetOptStat ("");
  gStyle -> SetStatX (.90);
  gStyle -> SetStatY (.90);
  gStyle -> SetStatW (.15);

  g1 -> SetTitle(Title.c_str());
  g1 -> GetXaxis() -> SetLabelSize(0.04);
  g1 -> GetYaxis() -> SetLabelSize(0.04);
  g1 -> GetXaxis() -> SetTitleSize(0.05);
  g1 -> GetYaxis() -> SetTitleSize(0.05);
  g1 -> GetXaxis() -> SetTitleOffset(0.8);
  g1 -> GetYaxis() -> SetTitleOffset(1.2);
  g1 -> GetYaxis() -> SetRangeUser(ymin,ymax);
  g1 -> GetXaxis() -> SetRangeUser(xmin,xmax);
 
  g1 -> GetXaxis() -> SetTitle("iRing");
  g1 -> GetYaxis() -> SetTitle((s + " occupancy (hits/events)").c_str());
   
  g1 -> SetMarkerStyle(20);
  g1 -> SetMarkerSize(0.6);
  g1 -> SetMarkerColor(kBlack);
  g1 -> SetLineColor(kBlack);
  g1 -> SetLineWidth(1.8);

  g2 -> SetMarkerStyle(20);
  g2 -> SetMarkerSize(0.6);
  g2 -> SetMarkerColor(kBlue);
  g2 -> SetLineColor(kBlue);
  g2 -> SetLineWidth(1.8);

  TLegend* legend = new TLegend(.75, 0.84, 1., 0.96);
  legend -> SetFillColor(kWhite);
  legend -> SetFillStyle(1000);
  legend -> SetLineWidth(0); 
  legend -> SetLineColor(kWhite);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.04);
  legend -> AddEntry(g1,g1_Title.c_str(),"L");
  legend -> AddEntry(g2,g2_Title.c_str(),"L");

  TCanvas* c1 = new TCanvas("c1","c1");
  c1 -> cd();
  c1->SetLeftMargin(0.15);

  TGaxis::SetMaxDigits(2);

  g1 -> Draw("APL");
  g2 -> Draw("PL");
  legend -> Draw("same");

  TLine *line_1 = new TLine(xmin, lineValue_1, xmax, lineValue_1);
  line_1->DrawLine(xmin, lineValue_1, xmax, lineValue_1);
  line_1->SetLineColor(kGreen+2);
  line_1->SetLineWidth(2);
  line_1->Draw("same");

  TLine *line_2 = new TLine(xmin, lineValue_2, xmax, lineValue_2);
  line_2->DrawLine(xmin, lineValue_2, xmax, lineValue_2);
  line_2->SetLineColor(kRed);
  line_2->SetLineWidth(2);
  line_2->Draw("same");

  TLatex *text_1 = new TLatex();
  text_1->SetNDC();
  text_1->SetTextColor(kGreen+2);
  text_1->SetTextSize(0.034);
  text_1->DrawLatex(0.22,0.81,"#splitline{green line: occupancy value to calibrate in 10 hours}{(statistical precision = 0.1%, rate = 1.5 kHz)}");
  text_1->Draw("same");
  
  TLatex *text_2 = new TLatex();
  text_2->SetNDC();
  text_2->SetTextColor(kRed);
  text_2->SetTextSize(0.034);
  if (s=="mean")
    text_2->DrawLatex(0.22,0.71,("#splitline{red line: total mean occupancy}{#Rightarrow time to calibrate = " + time + " h}").c_str());
  if (s=="lowest")
    text_2->DrawLatex(0.22,0.71,("#splitline{red line: lowest occupancy}{#Rightarrow time to calibrate = " + time + " h}").c_str());
  text_2->Draw("same");
  
  c1 -> Print((s + "OccupancyVSring_" + Title + ".png").c_str(),"png");
  c1 -> Print((s + "OccupancyVSring_" + Title + ".pdf").c_str(),"pdf");

  delete c1;

}

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
