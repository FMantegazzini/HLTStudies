/*
Macro to draw and fit event_size vs PU
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

void EvSizePlots () {

  TGraph *g1 = new TGraph();

  g1->SetPoint(0,20,9.7);
  g1->SetPoint(1,40,16.1);
  
  g1 -> SetTitle("bx50ns_EB_3_EE_5");
  g1 -> GetXaxis() -> SetLabelSize(0.04);
  g1 -> GetYaxis() -> SetLabelSize(0.04);
  g1 -> GetXaxis() -> SetTitleSize(0.05);
  g1 -> GetYaxis() -> SetTitleSize(0.05);
  g1 -> GetYaxis() -> SetTitleOffset(1.);
  g1 -> GetYaxis() -> SetRangeUser(4.,30.);
  g1 -> GetXaxis() -> SetRangeUser(-1.,50.);
 
  g1 -> GetXaxis() -> SetTitle("PU");
  g1 -> GetYaxis() -> SetTitle("Event Size(kB)");
   
  g1 -> SetMarkerStyle(20);
  g1 -> SetMarkerSize(1.5);
  g1 -> SetMarkerColor(kBlack);
  g1 -> SetLineColor(kBlack);
  g1 -> SetLineWidth(1.8);

  gStyle->SetOptFit(1111);
  
  TF1 *f = new TF1("linear", "[0]+[1]*x", -1, 50);
  g1 -> Fit("linear","","", -1, 50);

  TCanvas* c1 = new TCanvas("c1","c1");
  c1 -> cd();
  g1 -> Draw("APL");
  c1 -> Print("evsize_PU.png","png");
  c1 -> Print("evsize_PU.pdf","pdf");

}
