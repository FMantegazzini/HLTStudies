#include "MakeAlCaPhiSymSpectra.h"
#include "HLTStudies/AlCaPhiSymStudies/interface/TEndcapRings.h"

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

#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsMCRcd.h"

#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"

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

MakeAlCaPhiSymSpectra::MakeAlCaPhiSymSpectra(const edm::ParameterSet& ps){

  recHitCollection_EB_                   = ps.getParameter<edm::InputTag>("recHitCollection_EB");
  recHitCollection_EE_                   = ps.getParameter<edm::InputTag>("recHitCollection_EE");

  enMin = 0.;
  enMax = 8.;
  calMin = 0.;
  calMax = 20;
  
  naiveId_ = 0;

  edm::Service<TFileService> fs;

  h_nEvents = fs->make<TH1F>("h_nEvents","h_nEvents",3,-1,2);
  std::cout << "h_nEvents created" << std::endl;

  EBM_eSpectrum_histos.resize(EB_rings);
  EBP_eSpectrum_histos.resize(EB_rings);
  EEM_eSpectrum_histos.resize(EE_rings);
  EEP_eSpectrum_histos.resize(EE_rings);

  EBM_etSpectrum_histos.resize(EB_rings);
  EBP_etSpectrum_histos.resize(EB_rings);
  EEM_etSpectrum_histos.resize(EE_rings);
  EEP_etSpectrum_histos.resize(EE_rings);

  EBM_calibration_histos.resize(EB_rings);
  EBP_calibration_histos.resize(EB_rings);
  EEM_calibration_histos.resize(EE_rings);
  EEP_calibration_histos.resize(EE_rings);

  std::cout << "histos vectors resized" << std::endl;

  std::ostringstream t;
  for (int i=0; i<EB_rings; i++) { //EB

    t << "EBM_eSpectrum_" << i+1;
    EBM_eSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax); 
    t.str("");
    t << "EBM_etSpectrum_" << i+1;
    EBM_etSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax); 
    t.str("");
    t << "EBM_calibration_" << i+1;
    EBM_calibration_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,calMin,calMax); 
    t.str("");

    t << "EBP_eSpectrum_" << i+1;
    EBP_eSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax); 
    t.str("");
    t << "EBP_etSpectrum_" << i+1;
    EBP_etSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax); 
    t.str("");
    t << "EBP_calibration_" << i+1;
    EBM_calibration_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax); 
    t.str("");

  }

  for (int i=0; i<EE_rings; i++) { //EE
    t << "EEM_eSpectrum_" << i+1;
    EEM_eSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax); 
    t.str("");
    t << "EEM_etSpectrum_" << i+1;
    EEM_etSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax); 
    t.str("");
    t << "EEM_calibration_" << i+1;
    EEM_calibration_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,calMin,calMax); 
    t.str("");

    t << "EEP_eSpectrum_" << i+1;
    EEP_eSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax); 
    t.str("");
    t << "EEP_etSpectrum_" << i+1;
    EEP_etSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax); 
    t.str("");
    t << "EEP_calibration_" << i+1;
    EEP_calibration_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,calMin,calMax); 
    t.str("");

  }
 
  std::cout << "histos created" << std::endl; 
  
}

MakeAlCaPhiSymSpectra::~MakeAlCaPhiSymSpectra()
{
}

// ------------ method called once each job just before starting event loop  ------------

void MakeAlCaPhiSymSpectra::beginJob()
{
}

// ------------ method called to for each event  ------------
void MakeAlCaPhiSymSpectra::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
  naiveId_++;
  std::cout << "Event = " << naiveId_ << std::endl;
  
  using namespace edm;
  /*
  //---LaserCorrections
  edm::ESHandle<EcalLaserDbService> theLaser;
  iSetup.get<EcalLaserDbRecord>().get(theLaser);
  std::cout << ">>> Handle LaserCorrections" << std::endl;
 
  //---InterCalibration constant s
  edm::ESHandle<EcalIntercalibConstantsMC> theIC;
  iSetup.get<EcalIntercalibConstantsMCRcd>().get(theIC) ;
  const EcalIntercalibConstantsMC* ical = theIC.product();
  const EcalIntercalibConstantMCMap &icalMap = ical->getMap();
  std::cout << ">>> Handle IC" << std::endl;

  //---ADCToGeV constants
  edm::ESHandle<EcalADCToGeVConstant> theADCToGeV;
  iSetup.get<EcalADCToGeVConstantRcd>().get(theADCToGeV);
  const EcalADCToGeVConstant* agc = theADCToGeV.product();
  std::cout << ">>> Handle ADCToGeV" << std::endl;

  const edm::Timestamp& evtTimeStamp = edm::Timestamp(0);
  */
  TEndcapRings *eRings = new TEndcapRings(); 
  
  //---rechitsEB
  edm::Handle<EcalRecHitCollection> recHitsEB;
  ev.getByLabel( recHitCollection_EB_, recHitsEB );
  const EcalRecHitCollection* theBarrelEcalRecHits = recHitsEB.product () ;
  if ( ! recHitsEB.isValid() ) {
    std::cerr << "EcalValidation::analyze --> recHitsEB not found" << std::endl; 
  }
  std::cout << ">>> Handle EBRechits" << std::endl;
  
  //---rechitsEE
  edm::Handle<EcalRecHitCollection> recHitsEE;
  ev.getByLabel( recHitCollection_EE_, recHitsEE );
  const EcalRecHitCollection* theEndcapEcalRecHits = recHitsEE.product () ;
  if ( ! recHitsEE.isValid() ) {
    std::cerr << "EcalValidation::analyze --> recHitsEE not found" << std::endl; 
  }
  std::cout << ">>> Handle EERechit" << std::endl;
    
  //---EBRechits iteration
  EBRecHitCollection::const_iterator itb;
  std::cout << ">>> start EBRechits iteration " << std::endl;
  for (itb = theBarrelEcalRecHits->begin(); itb != theBarrelEcalRecHits->end(); ++itb)
    {
      EBDetId id_crystal(itb->id());

      float ieta = id_crystal.ieta();
      float eta = eRings->GetEtaFromIRing (ieta);
      std::cout << "ieta = " << ieta << ", eta = " << eta << std::endl;
            
      float e  = itb->energy();
      float et = itb->energy()/cosh(eta);
      std::cout << "e = " << e << ", et = " << et << std::endl;
      /*    
      float LaserCorrection = theLaser->getLaserCorrection(id_crystal, evtTimeStamp);
      std::cout << "LaserCorrection = " << LaserCorrection << std::endl;
      float InterCalibConst = 1.;
      EcalIntercalibConstantMCMap::const_iterator icalit = icalMap.find(id_crystal);
      if( icalit!=icalMap.end() )
	{
	  InterCalibConst = (*icalit);
	}
      std::cout << "IC = " << InterCalibConst << std::endl;

      float ADCToGeV_EB = agc->getEBValue();
      std::cout << "ADCToGeV = " << ADCToGeV_EB << std::endl;

      float Calibration = LaserCorrection * InterCalibConst * ADCToGeV_EB;
      std::cout << "Calibration = " << Calibration << std::endl;
      
      if (ieta < 0) { //EBM
	EBM_eSpectrum_histos[ieta+85]->Fill(e);
	EBM_etSpectrum_histos[ieta+85]->Fill(et);
	//EBM_calibration_histos[ieta+85]->Fill(Calibration);
      }
      else if (ieta > 0) { //EBP
	EBP_eSpectrum_histos[ieta-1]->Fill(e);
	EBP_etSpectrum_histos[ieta-1]->Fill(et);
	//EBP_calibration_histos[ieta-1]->Fill(Calibration);
      }
      std::cout << ">>> histos filled" << std::endl;*/
    }
  
  std::cout << "EBRechits iteration finished" << std::endl;

  //---EERechits iteration
  EERecHitCollection::const_iterator ite;
  std::cout << ">>> start EERechits iteration " << std::endl;
  for (ite = theEndcapEcalRecHits->begin(); ite != theEndcapEcalRecHits->end(); ++ite)
    {
      EEDetId id_crystal(ite->id());

      float iring = eRings->GetEndcapRing( id_crystal.ix(), id_crystal.iy(), id_crystal.zside() );
      float eta = eRings->GetEtaFromIRing (iring);
      std::cout << "iring = " << iring << ", eta = " << eta << std::endl;
      
      float e  = ite->energy();
      float et = ite->energy()/cosh(eta);
      std::cout << "e = " << e << ", et = " << et << std::endl;

      /*
      float LaserCorrection = theLaser->getLaserCorrection(id_crystal, evtTimeStamp);
      
      float InterCalibConst = 1.;
      EcalIntercalibConstantMCMap::const_iterator icalit = icalMap.find(id_crystal);
      if( icalit!=icalMap.end() )
	{
	  InterCalibConst = (*icalit);
	}
  
      float ADCToGeV_EE = agc->getEEValue();
      float Calibration = LaserCorrection * InterCalibConst * ADCToGeV_EE;
      */
      /*
      if (id_crystal.zside() < 0) { //EEM
	EEM_eSpectrum_histos[iring]->Fill(e);
	EEM_etSpectrum_histos[iring]->Fill(et);
	//EEM_calibration_histos[iring]->Fill(Calibration);
      }
      else if (id_crystal.zside() > 0) { //EEP
	EEP_eSpectrum_histos[iring]->Fill(e);
	EEP_etSpectrum_histos[iring]->Fill(et);
	//EEP_calibration_histos[iring]->Fill(Calibration);
      }
      */
    }  
  
  std::cout << "EERechits iteration finished" << std::endl;

}
  

// ------------ method called once each job just after ending the event loop  ------------
void MakeAlCaPhiSymSpectra::endJob()
{
  h_nEvents->SetBinContent(h_nEvents->FindBin(0),naiveId_); 
}


//define this as a plug-in
DEFINE_FWK_MODULE(MakeAlCaPhiSymSpectra);
