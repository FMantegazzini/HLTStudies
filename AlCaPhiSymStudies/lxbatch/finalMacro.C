/*
Macro to draw occupancy plots.
This macro read outputTotal.root, that is the union of the .root files created by the analyzer.
*/

#include <iostream>
#include <stdio>
#include <fstream>
#include <iostream>
#include <string.h>
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"

void finalMacro () {

  TFile* f = new TFile("file:outputOccupancy_total.root"); 
  
  TH1F* h_events = (TH1F*)f->Get("h_nEvents");
  TH2F* h_EB = (TH2F*)f->Get("h2_hitOccupancy_EB");
  TH2F* h_EEM = (TH2F*)f->Get("h2_hitOccupancy_EEM");
  TH2F* h_EEP = (TH2F*)f->Get("h2_hitOccupancy_EEP");

  Double_t num_events = h_events->GetBinContent(h_events->FindBin(0));
  cout << "Total number of events = " << num_events << endl;

  h_EB->Scale(1./num_events);
  h_EEM->Scale(1./num_events);
  h_EEP->Scale(1./num_events);
  cout << "Scaling for " << 1./num_events << endl;
  
  std::string PU;
  std::string bx;
  std::string yesno;
  std::string multifit;
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
    
  h_EB -> SetTitle(("EB_PU" + PU + "_" + bx + "ns_" + multifit).c_str());
  h_EB -> GetXaxis() -> SetTitle("iphi");
  h_EB -> GetYaxis() -> SetTitle("ieta");
  h_EB -> GetZaxis() -> SetTitle("hits/event");

  h_EB -> GetXaxis() -> SetLabelSize(0.04);
  h_EB -> GetXaxis() -> SetTitleSize(0.05);
  h_EB -> GetXaxis() -> SetTitleOffset(0.9);
  h_EB -> GetYaxis() -> SetLabelSize(0.04);
  h_EB -> GetYaxis() -> SetTitleSize(0.05);
  h_EB -> GetYaxis() -> SetTitleOffset(0.9);
  h_EB -> GetZaxis() -> SetLabelSize(0.04);
  h_EB -> GetZaxis() -> SetTitleSize(0.05);
  h_EB -> GetZaxis() -> SetTitleOffset(1.2);
  h_EB -> GetZaxis() -> SetRangeUser(0,0.001);
  h_EB -> SetStats(0);

  h_EEM -> SetTitle(("EEM_PU" + PU + "_" + bx + "ns_" + multifit).c_str());
  h_EEM -> GetXaxis() -> SetTitle("ix");
  h_EEM -> GetYaxis() -> SetTitle("iy");
  h_EEM -> GetZaxis() -> SetTitle("hits/event");

  h_EEM -> GetXaxis() -> SetLabelSize(0.04);
  h_EEM -> GetXaxis() -> SetTitleSize(0.05);
  h_EEM -> GetXaxis() -> SetTitleOffset(0.9);
  h_EEM -> GetYaxis() -> SetLabelSize(0.04);
  h_EEM -> GetYaxis() -> SetTitleSize(0.05);
  h_EEM -> GetYaxis() -> SetTitleOffset(0.9);
  h_EEM -> GetZaxis() -> SetLabelSize(0.04);
  h_EEM -> GetZaxis() -> SetTitleSize(0.05);
  h_EEM -> GetZaxis() -> SetTitleOffset(1.);
  h_EEM -> GetZaxis() -> SetRangeUser(0,0.0002);
  h_EEM -> SetStats(0);

  h_EEP -> SetTitle(("EEP_PU" + PU + "_" + bx + "ns_" + multifit).c_str());
  h_EEP -> GetXaxis() -> SetTitle("ix");
  h_EEP -> GetYaxis() -> SetTitle("iy");
  h_EEP -> GetZaxis() -> SetTitle("hits/event");
  h_EEP -> GetXaxis() -> SetLabelSize(0.04);
  h_EEP -> GetXaxis() -> SetTitleSize(0.05);
  h_EEP -> GetXaxis() -> SetTitleOffset(0.9);
  h_EEP -> GetYaxis() -> SetLabelSize(0.04);
  h_EEP -> GetYaxis() -> SetTitleSize(0.05);
  h_EEP -> GetYaxis() -> SetTitleOffset(0.9);
  h_EEP -> GetZaxis() -> SetLabelSize(0.04);
  h_EEP -> GetZaxis() -> SetTitleSize(0.05);
  h_EEP -> GetZaxis() -> SetTitleOffset(1.);
  h_EEP -> GetZaxis() -> SetRangeUser(0,0.0002);
  h_EEP -> SetStats(0);

  TCanvas* c1 = new TCanvas("c1","c1");
  c1->SetRightMargin(0.2);
  c1 -> cd();
  h_EB -> Draw("colz");
  c1 -> Print(("EB_PU" + PU + "_" + bx + "ns_" + multifit + ".png").c_str(), "png");
  c1 -> Print(("EB_PU" + PU + "_" + bx + "ns_" + multifit + ".pdf").c_str(), "pdf"); 
  
  TCanvas* c2 = new TCanvas("c2","c2");
  c2->SetRightMargin(0.2);
  c2 -> cd();
  h_EEM -> Draw("colz");
  c2 -> Print(("EEM_PU" + PU + "_" + bx + "ns_" + multifit + ".png").c_str(), "png");
  c2 -> Print(("EEM_PU" + PU + "_" + bx + "ns_" + multifit + ".pdf").c_str(), "pdf"); 
  
  TCanvas* c3 = new TCanvas("c3","c3");
  c3->SetRightMargin(0.2);
  c3 -> cd();
  h_EEP -> Draw("colz");
  c3 -> Print(("EEP_PU" + PU + "_" + bx + "ns_" + multifit + ".png").c_str(), "png");
  c3 -> Print(("EEP_PU" + PU + "_" + bx + "ns_" + multifit + ".pdf").c_str(), "pdf"); 
  
  //computation of mean hits per event and of lower occupancy

  Double_t EB_mean = 0.;
  Double_t EEM_mean= 0.;
  Double_t EEP_mean = 0.;
  Double_t EB_lowerOccupancy = 10.;
  Double_t EEM_lowerOccupancy = 10.;
  Double_t EEP_lowerOccupancy = 10.;

  int EB_channels = 0;
  for (int i=1; i <= h_EB->GetNbinsX(); i++) {
    for (int j=1; j <= h_EB->GetNbinsY(); j++) {
      if (h_EB->GetBinContent(h_EB->FindBin(i,j)) > 0) { //find mean occupancy
	EB_channels = EB_channels + 1;
	EB_mean = EB_mean + h_EB->GetBinContent(h_EB->FindBin(i,j));
      }
      
      if (h_EB->GetBinContent(h_EB->FindBin(i,j)) < EB_lowerOccupancy && h_EB->GetBinContent(h_EB->FindBin(i,j)) != 0) {  //Find lower occupancy
	EB_lowerOccupancy = h_EB->GetBinContent(h_EB->FindBin(i,j));
      }
    }
  }
  EB_mean = EB_mean/EB_channels;

 int EEM_channels = 0;
  for (int i=1; i <= h_EEM->GetNbinsX(); i++) {
    for (int j=1; j <= h_EEM->GetNbinsY(); j++) {
      if (h_EEM->GetBinContent(h_EEM->FindBin(i,j)) > 0) { //find lower occupancy
	EEM_channels = EEM_channels + 1;
	EEM_mean = EEM_mean + h_EEM->GetBinContent(h_EEM->FindBin(i,j));
      }
      if (h_EEM->GetBinContent(h_EEM->FindBin(i,j)) < EEM_lowerOccupancy && h_EEM->GetBinContent(h_EEM->FindBin(i,j)) != 0) //find lower occupancy
	EEM_lowerOccupancy = h_EEM->GetBinContent(h_EEM->FindBin(i,j));
    }
  }
  EEM_mean = EEM_mean/EEM_channels;

 int EEP_channels = 0;
  for (int i=1; i <= h_EEP->GetNbinsX(); i++) {
    for (int j=1; j <= h_EEP->GetNbinsY(); j++) {
      if (h_EEP->GetBinContent(h_EEP->FindBin(i,j)) > 0) { //find mean occupancy
	EEP_channels = EEP_channels + 1;
	EEP_mean = EEP_mean + h_EEP->GetBinContent(h_EEP->FindBin(i,j));
      }
      if (h_EEP->GetBinContent(h_EEP->FindBin(i,j)) < EEP_lowerOccupancy && h_EEP->GetBinContent(h_EEP->FindBin(i,j)) != 0) //find lower occupancy
	EEP_lowerOccupancy = h_EEP->GetBinContent(h_EEP->FindBin(i,j));
    }
  }
  EEP_mean = EEP_mean/EEP_channels;

  cout << "EB mean occupancy = " << EB_mean << endl;
  cout << "EEM mean occupancy = " << EEM_mean << endl;
  cout << "EEP mean occupancy = " << EEP_mean << endl;
  
  cout << "EB lower occupancy = " << EB_lowerOccupancy << endl;
  cout << "EEM lower occupancy = " << EEM_lowerOccupancy << endl;
  cout << "EEP lower occupancy = " << EEP_lowerOccupancy << endl;
}
