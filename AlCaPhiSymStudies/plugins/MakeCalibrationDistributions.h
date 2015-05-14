#ifndef MakeCalibrationDistributions_h
#define MakeCalibrationDistributions_h

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

class MakeCalibrationDistributions :  public edm::EDAnalyzer
{

 public:

  /// Constructor
  MakeCalibrationDistributions( const edm::ParameterSet& iConfig );
  
  /// Destructor
  ~MakeCalibrationDistributions();

 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------additional functions------------------
  // nothing

  // ----------member data ---------------------------
  edm::InputTag recHitCollection_EB_;
  edm::InputTag recHitCollection_EE_;

  const float ADCtoGeV_EB = 0.038940; //??
  const float ADCtoGeV_EE = 0.062850; //??
  int naiveId_;      

  // ------------- HISTOGRAMS ------------------------------------
  TH1F* h_nEvents;
  //vettori di histo--> da aggiungere
};

#endif
