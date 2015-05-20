#ifndef MakeAlCaPhiSymSpectra_h
#define MakeAlCaPhiSymSpectra_h

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TProfile.h"

// Framework
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ProducerBase.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HLTStudies/AlCaPhiSymStudies/interface/TEndcapRings.h"

class MakeAlCaPhiSymSpectra :  public edm::EDAnalyzer
{

 public:

  /// Constructor
  MakeAlCaPhiSymSpectra( const edm::ParameterSet& iConfig );
  
  /// Destructor
  ~MakeAlCaPhiSymSpectra();

 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------additional functions------------------
  // nothing

  // ----------member data ---------------------------
  edm::InputTag recHitCollection_EB_;
  edm::InputTag recHitCollection_EE_;

  const int EB_rings = 85;
  const int EE_rings = 39;

  float ebCut_ADC = 13.;
  float eeCut_ADC = 20.;

  float enMin;
  float enMax_EB;
  float enMax_EE;
  float calMin;
  float calMax;

  int nBins;
  int naiveId_;

  TEndcapRings* eRings;

  std::map<int,std::map<int,std::map<int,TH1F*> > > EBmap;
  std::map<int,std::map<int,std::map<int,TH1F*> > > EEmap;
 
  // ------------- HISTOGRAMS ------------------------------------
  TH1F* h_nEvents;

  TH2F* h2_hitOccupancy_EB = new TH2F("h2_hitOccupancy_EB_highCut1","h2_hitOccupancy_EB_highCut1",360,0,360,170,-85,85);
  TH2F* h2_hitOccupancy_EB = new TH2F("h2_hitOccupancy_EB_highCut2","h2_hitOccupancy_EB_highCut2",360,0,360,170,-85,85);
  TH2F* h2_hitOccupancy_EEP = new TH2F("h2_hitOccupancy_EEP_highCut1","h2_hitOccupancy_EEP_highCut1",100,1,101,100,1,101);
  TH2F* h2_hitOccupancy_EEP = new TH2F("h2_hitOccupancy_EEP_highCut2","h2_hitOccupancy_EEP_highCut2",100,1,101,100,1,101);
  TH2F* h2_hitOccupancy_EEM = new TH2F("h2_hitOccupancy_EEM_highCut1","h2_hitOccupancy_EEM_highCut1",100,1,101,100,1,101);
  TH2F* h2_hitOccupancy_EEM = new TH2F("h2_hitOccupancy_EEM_highCut2","h2_hitOccupancy_EEM_highCut2",100,1,101,100,1,101);
  
  std::vector<TH1F*> EBM_eSpectrum_histos;
  std::vector<TH1F*> EBP_eSpectrum_histos;
  std::vector<TH1F*> EEM_eSpectrum_histos;
  std::vector<TH1F*> EEP_eSpectrum_histos;

  std::vector<TH1F*> EBM_etSpectrum_histos;
  std::vector<TH1F*> EBP_etSpectrum_histos;
  std::vector<TH1F*> EEM_etSpectrum_histos;
  std::vector<TH1F*> EEP_etSpectrum_histos;

  std::vector<TH1F*> EBM_calibration_histos;
  std::vector<TH1F*> EBP_calibration_histos;
  std::vector<TH1F*> EEM_calibration_histos;
  std::vector<TH1F*> EEP_calibration_histos;

};

#endif
