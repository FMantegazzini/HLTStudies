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
#include <map>

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

  h2_hitOccupancy_EB_highCut1 = fs->make<TH2F>("h2_hitOccupancy_EB_highCut1","h2_hitOccupancy_EB_highCut1",360,0,360,170,-85,85);
  h2_hitOccupancy_EB_highCut2 = fs->make<TH2F>("h2_hitOccupancy_EB_highCut2","h2_hitOccupancy_EB_highCut2",360,0,360,170,-85,85);
  h2_hitOccupancy_EEM_highCut1 = fs->make<TH2F>("h2_hitOccupancy_EEM_highCut1","h2_hitOccupancy_EEM_highCut1",100,1,101,100,1,101);
  h2_hitOccupancy_EEM_highCut2 = fs->make<TH2F>("h2_hitOccupancy_EEM_highCut2","h2_hitOccupancy_EEM_highCut2",100,1,101,100,1,101);
  h2_hitOccupancy_EEP_highCut1 = fs->make<TH2F>("h2_hitOccupancy_EEP_highCut1","h2_hitOccupancy_EEP_highCut1",100,1,101,100,1,101);
  h2_hitOccupancy_EEP_highCut2 = fs->make<TH2F>("h2_hitOccupancy_EEP_highCut2","h2_hitOccupancy_EEP_highCut2",100,1,101,100,1,101);

  h2_calib_EB = fs->make<TH2F>("h2_calib_EB","calibration_EB",360,0,360,170,-85,85);
  h2_calib_EEM = fs->make<TH2F>("h2_calib_EEM","calibration_EEM",100,1,101,100,1,101);
  h2_calib_EEP = fs->make<TH2F>("h2_calib_EEP","calibration_EEP",100,1,101,100,1,101);

  h2_LC_EB = fs->make<TH2F>("h2_LC_EB","LC_EB",360,0,360,170,-85,85);
  h2_LC_EEM = fs->make<TH2F>("h2_LC_EEM","LC_EEM",100,1,101,100,1,101);
  h2_LC_EEP = fs->make<TH2F>("h2_LC_EEP","LC_EEP",100,1,101,100,1,101);

  h2_IC_EB = fs->make<TH2F>("h2_IC_EB","IC_EB",360,0,360,170,-85,85);
  h2_IC_EEM = fs->make<TH2F>("h2_IC_EEM","IC_EEM",100,1,101,100,1,101);
  h2_IC_EEP = fs->make<TH2F>("h2_IC_EEP","IC_EEP",100,1,101,100,1,101);

  EEM_ix66_iy26_eSpectrum = fs->make<TH1F>("EEM_ix66_iy26_eSpectrum","EEM_ix66_iy26",1000,0.,15.);
  EEM_ix54_iy25_eSpectrum = fs->make<TH1F>("EEM_ix54_iy25_eSpectrum","EEM_ix54_iy25",1000,0.,15.);
  EEM_ix100_iy57_eSpectrum = fs->make<TH1F>("EEM_ix100_iy57_eSpectrum","EEM_ix100_iy57",1000,0.,15.);

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
  std::ostringstream title;
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
    EBP_calibration_histos[i] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,calMin,calMax); 
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
   
  //crystals energy spectra for EBring = 15
  for (int iphi=0; iphi<361; iphi++) { 
    int ieta = 15;
    int iz = 0;
    t << "EB_ieta_" << ieta << "_iphi_" << iphi;
    EBmap[ieta][iphi][iz] = fs->make<TH1F>(t.str().c_str(),t.str().c_str(),nBins,0.,10.);
    t.str("");
  }  
  
  //crystals energy spectra and calibration spectra for EEM ring 23
  eRings = new TEndcapRings();
  std::vector<int> ix_vector_r23;
  std::vector<int> iy_vector_r23; 

  for (int ix=1; ix<101; ix++) { 
    for (int iy=1; iy<101; iy++) {
      int iring = eRings->GetEndcapRing(ix,iy,-1);
      if (iring==23) {
	ix_vector_r23.push_back(ix);
	iy_vector_r23.push_back(iy);
      }
    }
  }
  
  int i = 0;
  for (unsigned int ii=0; ii<ix_vector_r23.size(); ii++) {
    int iz = -1;
 
    t << "EEM_energy_ring_23_" << i+1; //calibration spectra
    title << "EEM_ring23_ix_" << ix_vector_r23.at(ii) << "_iy_" << iy_vector_r23.at(ii);
    EEmap_energy[ix_vector_r23.at(ii)][iy_vector_r23.at(ii)][iz] = fs->make<TH1F>(t.str().c_str(),title.str().c_str(),nBins,0.,80.);
    t.str("");
    title.str("");

    t << "EEM_calib_ring_23_" << i+1; //calibration spectra
    title << "EEM_ring23_ix_" << ix_vector_r23.at(ii) << "_iy_" << iy_vector_r23.at(ii);
    EEmap_calib[ix_vector_r23.at(ii)][iy_vector_r23.at(ii)][iz] = fs->make<TH1F>(t.str().c_str(),title.str().c_str(),nBins,0.,80.);
    t.str("");
    title.str("");

    i++;
  }
   
  //std::cout << "histos created" << std::endl; 
  
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
 
  //---InterCalibration constant
  edm::ESHandle<EcalIntercalibConstantsMC> theIC;
  iSetup.get<EcalIntercalibConstantsMCRcd>().get(theIC) ;
  const EcalIntercalibConstantsMC* ical = theIC.product();
  const EcalIntercalibConstantMCMap &icalMap = ical->getMap();
 
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
  const EcalChannelStatus& channelStatus = *csHandle; 
  
  //---rechitsEB
  edm::Handle<EcalRecHitCollection> recHitsEB;
  ev.getByLabel( recHitCollection_EB_, recHitsEB );
  const EcalRecHitCollection* theBarrelEcalRecHits = recHitsEB.product () ;
  if ( ! recHitsEB.isValid() ) {
    std::cerr << "EcalValidation::analyze --> recHitsEB not found" << std::endl; 
  }
  
  //---rechitsEE
  edm::Handle<EcalRecHitCollection> recHitsEE;
  ev.getByLabel( recHitCollection_EE_, recHitsEE );
  const EcalRecHitCollection* theEndcapEcalRecHits = recHitsEE.product () ;
  if ( ! recHitsEE.isValid() ) {
    std::cerr << "EcalValidation::analyze --> recHitsEE not found" << std::endl; 
  }
    
  //---EBRechits iteration
  EBRecHitCollection::const_iterator itb;
  std::cout << ">>> start EBRechits iteration " << std::endl;
  for (itb = theBarrelEcalRecHits->begin(); itb != theBarrelEcalRecHits->end(); ++itb)
    {
      EBDetId id_crystal(itb->id());

      int ieta = id_crystal.ieta();
      int iphi = id_crystal.iphi();
      int iz = 0;
      float eta = eRings->GetEtaFromIRing (ieta);

      uint16_t statusCode = 0;
      statusCode = channelStatus[itb->id().rawId()].getStatusCode();
      if (statusCode != 0)
	std::cout << "BAD CHANNEL STATUS: EB ieta = " << ieta << ", iphi = " << iphi << " --> ChStatus = " << statusCode << std::endl;
      if (statusCode != 0) continue;
            
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

      //calibration TH2F
      if ( h2_LC_EB->GetBinContent(h2_LC_EB->FindBin(id_crystal.iphi(),id_crystal.ieta())) == 0 )
	h2_LC_EB->SetBinContent(h2_LC_EB->FindBin(id_crystal.iphi(),id_crystal.ieta()),LaserCorrection);

      if ( h2_IC_EB->GetBinContent(h2_IC_EB->FindBin(id_crystal.iphi(),id_crystal.ieta())) == 0 )
     	h2_IC_EB->SetBinContent(h2_IC_EB->FindBin(id_crystal.iphi(),id_crystal.ieta()),InterCalibConst);

      if ( h2_calib_EB->GetBinContent(h2_calib_EB->FindBin(id_crystal.iphi(),id_crystal.ieta())) == 0 )
     	h2_calib_EB->SetBinContent(h2_calib_EB->FindBin(id_crystal.iphi(),id_crystal.ieta()),Calibration);
      
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
	EBmap[ieta][iphi][iz] -> Fill(e);
      }
    }
  
  std::cout << "EBRechits iteration finished" << std::endl;

  //---EERechits iteration
  EERecHitCollection::const_iterator ite;
  std::cout << ">>> start EERechits iteration " << std::endl;
  for (ite = theEndcapEcalRecHits->begin(); ite != theEndcapEcalRecHits->end(); ++ite)
    {
      EEDetId id_crystal(ite->id());
      
      int ix = id_crystal.ix();
      int iy = id_crystal.iy();
      int iz =  id_crystal.zside();
      int iring = eRings->GetEndcapRing(ix,iy,iz);
      float eta = eRings->GetEtaFromIRing (iring);

      uint16_t statusCode = 0;
      statusCode = channelStatus[ite->id().rawId()].getStatusCode();
      if (statusCode != 0)
	std::cout << "BAD CHANNEL STATUS: EE ix = " << ix << ", iy = " << iy << ", iz = " << iz << " --> ChStatus = " << statusCode << std::endl;
      if (statusCode != 0) continue;
     
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

      //calibration TH2F
     
      if (id_crystal.zside() < 0) { //EEM
	  
	  if ( h2_LC_EEM->GetBinContent(h2_LC_EEM->FindBin(id_crystal.ix(),id_crystal.iy())) == 0 )
	    h2_LC_EEM->SetBinContent(h2_LC_EEM->FindBin(id_crystal.ix(),id_crystal.iy()),LaserCorrection);

	  if ( h2_IC_EEM->GetBinContent(h2_IC_EEM->FindBin(id_crystal.ix(),id_crystal.iy())) == 0 )
	    h2_IC_EEM->SetBinContent(h2_IC_EEM->FindBin(id_crystal.ix(),id_crystal.iy()),InterCalibConst);

	  if ( h2_calib_EEM->GetBinContent(h2_calib_EEM->FindBin(id_crystal.ix(),id_crystal.iy())) == 0 )
	    h2_calib_EEM->SetBinContent(h2_calib_EEM->FindBin(id_crystal.ix(),id_crystal.iy()),Calibration);
	}

      if (id_crystal.zside() > 0) { //EEP

	  if ( h2_LC_EEP->GetBinContent(h2_LC_EEP->FindBin(id_crystal.ix(),id_crystal.iy())) == 0 )
	    h2_LC_EEP->SetBinContent(h2_LC_EEP->FindBin(id_crystal.ix(),id_crystal.iy()),LaserCorrection);
	 
	  if ( h2_IC_EEP->GetBinContent(h2_IC_EEP->FindBin(id_crystal.ix(),id_crystal.iy())) == 0 )
	    h2_IC_EEP->SetBinContent(h2_IC_EEP->FindBin(id_crystal.ix(),id_crystal.iy()),InterCalibConst);

	  if ( h2_calib_EEP->GetBinContent(h2_calib_EEP->FindBin(id_crystal.ix(),id_crystal.iy())) == 0 )
	    h2_calib_EEP->SetBinContent(h2_calib_EEP->FindBin(id_crystal.ix(),id_crystal.iy()),Calibration);
	}
	            
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

      //crystals energy distributions for ring = 23
      if (iring == 23 && iz == -1) {
	EEmap_energy[ix][iy][iz] -> Fill(e);
	}

      //crystals calibration distributions for ring = 23
      if (iring == 23 && iz == -1) {
	EEmap_calib[ix][iy][iz] -> Fill(Calibration);
      }

      //crystals energy distributions for some EEM crystals
      if (id_crystal.ix() == 66 && id_crystal.iy()== 26 && id_crystal.zside() < 0) //time = 62.529
	EEM_ix66_iy26_eSpectrum -> Fill (e);
      if (id_crystal.ix() == 54 && id_crystal.iy()== 25 && id_crystal.zside() < 0) //time = 31.2645
	EEM_ix54_iy25_eSpectrum -> Fill (e);
      if (id_crystal.ix() == 100 && id_crystal.iy()== 57 && id_crystal.zside() < 0) //time = 1.35933
	EEM_ix100_iy57_eSpectrum -> Fill (e);

      //print calibration values for some crystals (time to calibrated > 60 hours)
      if (id_crystal.ix() == 66 && id_crystal.iy()== 26 && id_crystal.zside() < 0)
	std::cout << "*** EEM ix = 66, iy = 26 --> IC = " << InterCalibConst << ", LC = " << LaserCorrection << ", calibration = " << Calibration << std::endl;
      if (id_crystal.ix() == 60 && id_crystal.iy()== 24 && id_crystal.zside() < 0)
	std::cout << "*** EEM ix = 60, iy = 24 --> IC = " << InterCalibConst << ", LC = " << LaserCorrection << ", calibration = " << Calibration << std::endl;
      if (id_crystal.ix() == 59 && id_crystal.iy()== 25 && id_crystal.zside() < 0)
	std::cout << "*** EEM ix = 59, iy = 25 --> IC = " << InterCalibConst << ", LC = " << LaserCorrection << ", calibration = " << Calibration << std::endl;
      
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
