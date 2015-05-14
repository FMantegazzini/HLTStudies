#include "MakeCalibrationDistributions.h"
#include "HLTStudies/AlCaPhiSymStudies/interface/TEndcapRings.h" //??

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/EventBase.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFramePlaybackInfoExtended.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/PCrossingFrame.h"


#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include <fstream>
#include <iostream>
#include <memory>

MakeCalibrationDistributions::MakeCalibrationDistributions(const edm::ParameterSet& ps){

  recHitCollection_EB_                   = ps.getParameter<edm::InputTag>("recHitCollection_EB");
  recHitCollection_EE_                   = ps.getParameter<edm::InputTag>("recHitCollection_EE");

  naiveId_ = 0;

  edm::Service<TFileService> fs;
  h_nEvents                    = fs->make<TH1F>("h_nEvents","h_nEvents",3,-1,2);

  h_calibration_EBM = fs->make<TH1F>("h_calibration_EBM","h_calibration_EBM",1000,0.,30.); // ma io voglio un vettore!!!
  h_calibration_EBP = fs->make<TH1F>("h_calibration_EBP","h_calibration_EBP",1000,0.,30.); // ma io voglio un vettore!!!
  h_calibration_EEM = fs->make<TH1F>("h_calibration_EEM","h_calibration_EEM",1000,0.,30.); // ma io voglio un vettore!!!
  h_calibration_EEP = fs->make<TH1F>("h_calibration_EEP","h_calibration_EEP",1000,0.,30.); // ma io voglio un vettore!!!

  //resize vettori di histos
  //rinomina histos
}
MakeCalibrationDistributions::~MakeCalibrationDistributions()
{
}

// ------------ method called once each job just before starting event loop  ------------

void MakeCalibrationDistributions::beginJob()
{
}

// ------------ method called to for each event  ------------
void MakeCalibrationDistributions::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
  naiveId_++;
  std::cout << "Event = " << naiveId_ << std::endl;
  
  using namespace edm;
  
  //LaserCorrections
  edm::ESHandle<EcalLaserDbService> theLaser;
  iSetup.get<EcalLaserDbRecord>().get(theLaser);
  
  //InterCalibration constants
  edm::ESHandle<EcalIntercalibConstants> pIcal;
  iSetup.get<EcalIntercalibConstantsRcd>().get(pIcal);
  const EcalIntercalibConstants* Mcal = pIcal.product();

  const edm::Timestamp& evtTimeStamp = edm::Timestamp(0);

  //rechitsEB
  edm::Handle<EcalRecHitCollection> recHitsEB;
  ev.getByLabel( recHitCollection_EB_, recHitsEB );
  const EcalRecHitCollection* theBarrelEcalRecHits = recHitsEB.product () ;
  if ( ! recHitsEB.isValid() ) {
    std::cerr << "EcalValidation::analyze --> recHitsEB not found" << std::endl; 
  }
  
  EBRecHitCollection::const_iterator itb;
  for (itb = ebRecHitsHandle->begin(); itb != ebRecHitsHandle->end(); ++itb)
    {
      EBDetId id_crystal(itb->id());          
      
      float LaserCorrection = theLaser->getLaserCorrection(id_crystal, evtTimeStamp);
      float InterCalibConst = 0.; // ??
      //float ADCtoGeV_EB = 0.;
      float Calibration = LaserCorrection * InterCalibConst * ADCToGeV;
      float ieta = id_crystal.ieta();
      if (ieta < 0) //EBM
	EBM_calibration_histos[ieta+85]->Fill(Calibration);
      else if (ieta > 0) //EBP
	EBP_calibration_histos[ieta-1]->Fill(Calibration);
    }

  //rechitsEE
  edm::Handle<EcalRecHitCollection> recHitsEE;
  ev.getByLabel( recHitCollection_EE_, recHitsEE );
  const EcalRecHitCollection* theEndcapEcalRecHits = recHitsEE.product () ;
  if ( ! recHitsEE.isValid() ) {
    std::cerr << "EcalValidation::analyze --> recHitsEE not found" << std::endl; 
  }
 
  EERecHitCollection::const_iterator ite;
  for (ite = eeRecHitsHandle->begin(); ite != eeRecHitsHandle->end(); ++ite)
    {
      EEDetId id_crystal(ite->id());

      float LaserCorrection = theLaser->getLaserCorrection(id_crystal, evtTimeStamp);
      float InterCalibConst = 0.; // ??
       //float ADCtoGeV_EE = 0.;
      float Calibration = LaserCorrection * InterCalibConst * ADCToGeV;
       
      float iring = eRings->GetEndcapRing( id_crystal.ix(), id_crystal.iy(), id_crystal.zside() );
            
      if (id_crystal.zside() < 0) //EEM
	EEM_calibration_histos[iring]->Fill(Calibration);
      else if (id_crystal.zside > 0) //EEP
	EEP_calibration_histos[iring]->Fill(Calibration);
    }
}

void MakeCalibrationDistributions::endJob()
{
  h_nEvents->SetBinContent(h_nEvents->FindBin(0),naiveId_); 
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeCalibrationDistributions);
