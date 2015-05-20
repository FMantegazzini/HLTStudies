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

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatusCode.h"

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
  enMax = 100.;
  calMin = 0.;
  calMax = 50.;

  nBins = 100000;  
  naiveId_ = 0;

  edm::Service<TFileService> fs;

  h_nEvents = fs->make<TH1F>("h_nEvents","h_nEvents",3,-1,2);

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

  std::ostringstream t;
  for (int i=0; i<EB_rings; i++) { //EB

    t << "EBM_eSpectrum_" << i+1;
    EBM_eSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax_EB); 
    t.str("");
    t << "EBM_etSpectrum_" << i+1;
    EBM_etSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax_EB); 
    t.str("");
    t << "EBM_calibration_" << i+1;
    EBM_calibration_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,calMin,calMax); 
    t.str("");

    t << "EBP_eSpectrum_" << i+1;
    EBP_eSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax_EB); 
    t.str("");
    t << "EBP_etSpectrum_" << i+1;
    EBP_etSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax_EB); 
    t.str("");
    t << "EBP_calibration_" << i+1;
    EBP_calibration_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,calMin,calMax); 
    t.str("");

  }

  for (int i=0; i<EE_rings; i++) { //EE
    t << "EEM_eSpectrum_" << i+1;
    EEM_eSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax_EE); 
    t.str("");
    t << "EEM_etSpectrum_" << i+1;
    EEM_etSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax_EE); 
    t.str("");
    t << "EEM_calibration_" << i+1;
    EEM_calibration_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,calMin,calMax); 
    t.str("");

    t << "EEP_eSpectrum_" << i+1;
    EEP_eSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax_EE); 
    t.str("");
    t << "EEP_etSpectrum_" << i+1;
    EEP_etSpectrum_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,enMin,enMax_EE); 
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
  eRings = new TEndcapRings();
}

// ------------ method called to for each event  ------------
void MakeAlCaPhiSymSpectra::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
  naiveId_++;
  std::cout << "Event = " << naiveId_ << std::endl;
  
  using namespace edm;
   
  //---LaserCorrections
  edm::ESHandle<EcalLaserDbService> theLaser;
  iSetup.get<EcalLaserDbRecord>().get(theLaser);
  std::cout << ">>> Handle LaserCorrections" << std::endl;
 
  //---InterCalibration constant
  edm::ESHandle<EcalIntercalibConstantsMC> theIC;
  iSetup.get<EcalIntercalibConstantsMCRcd>().get(theIC) ;
  const EcalIntercalibConstantsMC* ical = theIC.product();
  const EcalIntercalibConstantMCMap &icalMap = ical->getMap();
  std::cout << ">>> Handle IC" << std::endl;

  /*---ADCToGeV constants
  edm::ESHandle<EcalADCToGeVConstant> theADCToGeV;
  iSetup.get<EcalADCToGeVConstantRcd>().get(theADCToGeV);
  const EcalADCToGeVConstant* agc = theADCToGeV.product();
  std::cout << ">>> Handle ADCToGeV" << std::endl;
  */

  const edm::Timestamp& evtTimeStamp = edm::Timestamp(0);

  //---Channel Status
  edm::ESHandle<EcalChannelStatus> csHandle;
  iSetup.get<EcalChannelStatusRcd>().get(csHandle);
  // const EcalChannelStatus* theEcalChStatus = csHandle.product();
  const EcalChannelStatus& channelStatus = *csHandle; 
  
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

      uint16_t statusCode = 0;
      statusCode = channelStatus[itb->id().rawId()].getStatusCode();
      std::cout << "ChStatus = " << statusCode << std::endl;

      if (statusCode != 0) continue;

      int ieta = id_crystal.ieta();
      int iphi = id_crystal.iphi();
      int iz = id_crystal.iz(); //??
      float eta = eRings->GetEtaFromIRing (ieta);
            
      float e  = itb->energy();
      float et = itb->energy()/cosh(eta);

      float ebCut_GeV = ebCut_ADC*0.04;
      float et_thr_1 = ebCut_GeV/cosh(eta) + 1.;
      float et_thr_2 = ebCut_GeV/cosh(eta) + 2.;

      //occupancy plots
      if (e > ebCut_GeV && et < et_thr_1)
	h2_hitOccupancy_EB_highCut1->Fill(id_crystal.iphi(),id_crystal.ieta());
      if (e > ebCut_GeV && et < et_thr_2)
	h2_hitOccupancy_EB_highCut2->Fill(id_crystal.iphi(),id_crystal.ieta()); 
      
      //calibration
      float LaserCorrection = theLaser->getLaserCorrection(id_crystal, evtTimeStamp);
      
      float InterCalibConst = 1.;
      EcalIntercalibConstantMCMap::const_iterator icalit = icalMap.find(id_crystal);
      if( icalit!=icalMap.end() )
	{
	  InterCalibConst = (*icalit);
	}
      
      //float ADCToGeV_EB = agc->getEBValue();      
      float Calibration = LaserCorrection * InterCalibConst;
      //std::cout << "Calibration = " << Calibration << std::endl;
      
      //spectra
      if (ieta < 0) { //EBM
	EBM_eSpectrum_histos[ieta+85]->Fill(e);
	EBM_etSpectrum_histos[ieta+85]->Fill(et);
	EBM_calibration_histos[ieta+85]->Fill(Calibration);
      }
      else if (ieta > 0) { //EBP
	EBP_eSpectrum_histos[ieta-1]->Fill(e);
	EBP_etSpectrum_histos[ieta-1]->Fill(et);
	EBP_calibration_histos[ieta-1]->Fill(Calibration);
	}

      //crystals energy distributions for ring ieta = 15
      if (ieta == 15) {
	TH1F *h = new TH1F("h","h",100000,0,100); 
	EBmap[ieta][iphi][iz] = h;
	h->Fill(e);
      }
    }
  
  std::cout << "EBRechits iteration finished" << std::endl;

  //---EERechits iteration
  EERecHitCollection::const_iterator ite;
  std::cout << ">>> start EERechits iteration " << std::endl;
  for (ite = theEndcapEcalRecHits->begin(); ite != theEndcapEcalRecHits->end(); ++ite)
    {
      EEDetId id_crystal(ite->id());

      uint16_t statusCode = 0;
      statusCode = channelStatus[ite->id().rawId()].getStatusCode();
      std::cout << "ChStatus = " << statusCode << std::endl;

      if (statusCode != 0) continue;

      int ix = id_crystal.ix();
      int iy = id_crystal.iy();
      int iz =  id_crystal.zside());
      float iring = eRings->GetEndcapRing( id_crystal.ix(), id_crystal.iy(), id_crystal.zside() );
      float eta = eRings->GetEtaFromIRing (iring);
      
      float e  = ite->energy();
      float et = ite->energy()/cosh(eta);

      //occupancy plots
      if (id_crystal.zside() > 0) { //EEP
	float eepCut_GeV = eeCut_ADC*( 72.92+(3.549)*iring + (0.2262)*iring*iring )/1000.;
	float et_thr_1 = eepCut_GeV/cosh(eta) + 1.;
	float et_thr_2 = eepCut_GeV/cosh(eta) + 2.;
	if (e > eepCut_GeV && et < et_thr_1)
	  h2_hitOccupancy_EEP_highCut1->Fill(id_crystal.ix(),id_crystal.iy());
	if (e > eepCut_GeV && et < et_thr_2)
	  h2_hitOccupancy_EEP_highCut2->Fill(id_crystal.ix(),id_crystal.iy());
      }
            
      if (id_crystal.zside() < 0) { //EEM
	float eemCut_GeV = eeCut_ADC*( 79.29+(4.148)*iring + (0.2442)*iring*iring )/1000.;
	float et_thr_1 = eemCut_GeV/cosh(eta) + 1.;
	float et_thr_2 = eemCut_GeV/cosh(eta) + 2.;
	if (e > eemCut_GeV && et < et_thr_1)
	  h2_hitOccupancy_EEM_highCut1->Fill(id_crystal.ix(),id_crystal.iy());
	if (e > eemCut_GeV && et < et_thr_2)
	  h2_hitOccupancy_EEM_highCut2->Fill(id_crystal.ix(),id_crystal.iy());
      }

      //calibration
      float LaserCorrection = theLaser->getLaserCorrection(id_crystal, evtTimeStamp);
      
      float InterCalibConst = 1.;
      EcalIntercalibConstantMCMap::const_iterator icalit = icalMap.find(id_crystal);
      if( icalit!=icalMap.end() )
	{
	  InterCalibConst = (*icalit);
	}
  
      //float ADCToGeV_EE = agc->getEEValue();
      float Calibration = LaserCorrection * InterCalibConst;
      //std::cout << "Calibration = " << Calibration << std::endl; 
           
      //spectra
      if (id_crystal.zside() < 0) { //EEM
	EEM_eSpectrum_histos[iring]->Fill(e);
	EEM_etSpectrum_histos[iring]->Fill(et);
	EEM_calibration_histos[iring]->Fill(Calibration);
      }
      else if (id_crystal.zside() > 0) { //EEP
	EEP_eSpectrum_histos[iring]->Fill(e);
	EEP_etSpectrum_histos[iring]->Fill(et);
	EEP_calibration_histos[iring]->Fill(Calibration);
      }

      //crystals energy distributions for iring = 16
      if (iring == 16) {
	TH1F *h = new TH1F("h","h",100000,0,100); 
	EBmap[ix][iy][iz] = h;
	h->Fill(e);
      }
      
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
