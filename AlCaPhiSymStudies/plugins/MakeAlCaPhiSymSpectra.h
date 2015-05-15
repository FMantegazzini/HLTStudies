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

  const int nBins = 1000;
  const int EB_rings = 85;
  const int EE_rings = 38;

  float enMin;
  float enMax;
  float calMin;
  float calMax;

  int naiveId_;

  //ostringstream t;

  // ------------- HISTOGRAMS ------------------------------------
  TH1F* h_nEvents;
  
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
