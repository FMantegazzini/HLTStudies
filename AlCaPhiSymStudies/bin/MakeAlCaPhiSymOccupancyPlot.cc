#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "PhysicsTools/Utilities/macros/setTDRStyle.C"

using namespace std;

//function to get eta from iRing (i.e. ieta)
float TEndcapRings::GetEtaFromIRing(const int& iRing);

//**********MAIN**************************************************************************
int main( int argc, char *argv[] )
{
    setTDRStyle();
    gStyle->SetPadTickX( 1 );
    gStyle->SetPadTickY( 1 );
    gStyle->SetOptTitle( 0 );
    gStyle->SetOptFit( 1 );

    gSystem->Load( "libFWCoreFWLite" );
    AutoLibraryLoader::enable();

    if( argc < 2 ) {
        cout << "Usage : " << argv[0] << " [parameters.py]" << endl;
        return 0;
    }
    if( !edm::readPSetsFrom( argv[1] )->existsAs<edm::ParameterSet>( "process" ) ) {
        cout << " ERROR: ParametersSet 'process' is missing in your configuration file"
             << endl;
        return 0;
    }

    //-----get the python configuration-----
    const edm::ParameterSet &process = edm::readPSetsFrom( argv[1] )->getParameter<edm::ParameterSet>( "process" );
    const edm::ParameterSet &filesOpt = process.getParameter<edm::ParameterSet>( "ioFiles" );

    //---io files option---
    vector<string> filesList = filesOpt.getParameter<vector<string> >("inputFiles"); 
    std::cout << filesOpt.getParameter<string>("outputFile").c_str() << std::endl;
    TFile* outFile = new TFile(filesOpt.getParameter<string>("outputFile").c_str(), "recreate");
    
    TH2F* h2_hitOccupancy_EB = new TH2F("h2_hitOccupancy_EB","h2_hitOccupancy_EB",360,0,360,170,-85,85);
    TH2F* h2_hitOccupancy_EEP = new TH2F("h2_hitOccupancy_EEP","h2_hitOccupancy_EEP",100,1,101,100,1,101);
    TH2F* h2_hitOccupancy_EEM = new TH2F("h2_hitOccupancy_EEM","h2_hitOccupancy_EEM",100,1,101,100,1,101);
    TH1F* h_nEvents = new TH1F("h_nEvents","h_nEvents",3,-1,2); 

    float ebCut_ADC = 13.;
    float eeCut_ADC = 20.;

    int iEvent=0;
    for(unsigned int iFile=0; iFile<filesList.size(); iFile++)
    {
        TFile* inFile = TFile::Open(filesList.at(iFile).c_str());
	std::cout << " >>> " << filesList.at(iFile) << std::endl;

        //---FWLite interfaces---
        fwlite::Event event(inFile);
	fwlite::Handle<EcalRecHitCollection> ebRechitHandle;
	fwlite::Handle<EcalRecHitCollection> eeRechitHandle;

        //---events loop---
        for(event.toBegin(); !event.atEnd(); ++event) {
	  iEvent++;
        
	  // get EBRecHits
	  ebRecHitsHandle.getByLabel(event,"ecalRecHit","EcalRecHitsEB","PHISYM" );
	  ebRecHitsHandle.product();
	  
	  EBRecHitCollection::const_iterator itb;
	  for (itb = ebRecHitsHandle->begin(); itb != ebRecHitsHandle->end(); ++itb)
	    {
	      EBDetId id_crystal(itb->id());
	      
	      float ieta = id_crystal.ieta();
	      float eta = GetEtaFromIRing (ieta);
	      
	      float e  = itb->energy();
	      float et = itb->energy()/cosh(eta);
	      
	      float ebCut_GeV = ebCut_ADC*0.4*1.15;
	      float et_thr = ebCut_GeV/cosh(eta) + 1.;
	      
	      if (e > ebCut_GeV && et < et_thr) {
		h2_hitOccupancy_EB->Fill(id_crystal.iphi(),id_crystal.ieta()); 
	      }
	    }
	  
	  // get EERecHits
	  eeRecHitsHandle.getByLabel(event,"ecalRecHit","EcalRecHitsEE","PHISYM");
	  eeRecHitsHandle.product();
	  
	  EERecHitCollection::const_iterator ite;
	  for (ite = eeRecHitsHandle->begin(); ite != eeRecHitsHandle->end(); ++ite)
	    {
	      EEDetId id_crystal(ite->id());
	      
	      float ieta = id_crystal.ieta();
	      float eta = GetEtaFromIRing (ieta);
	      
	      float e  = ite->energy();
	      float et = ite->energy()/cosh(eta);
	      
	      if (id_crystal.zside() > 0) { //EEP
		float eepCut_GeV = eeCut_ADC*(72.92+(3.549)*ieta + (0.2262)*ieta*ieta)*1.15/1000.;
		float et_thr = eepCut_GeV/cosh(eta) + 1.;
		if (e > eepCut_GeV && et < et_thr)
		  h2_hitOccupancy_EEP->Fill(id_crystal.ix(),id_crystal.iy());
	      }
	      
	      if (id_crystal.zside() < 0) { //EEM
		float eemCut_GeV = eeCut_GeV*(79.29+(4.148)*ieta + (0.2442)*ieta*ieta)*1.15/1000.;
		float et_thr = eemCut_GeV/cosh(eta) + 1.;
		if (e > eemCut_GeV && et < et_thr)
		  h2_hitOccupancy_EEM->Fill(id_crystal.ix(),id_crystal.iy());
	      }
	      
	    }

	} //event loop
	        
    } //files list
    
    std::cout << "nEvents = " << iEvent << std::endl;
    
    h_nEvents->SetBinContent(h_nEvents->FindBin(0),iEvent);
    
    outFile->cd();

    h2_hitOccupancy_EB->Write();
    h2_hitOccupancy_EEP->Write();
    h2_hitOccupancy_EEM->Write();
    h_nEvents->Write();

    outFile->Close();

} //main

float GetEtaFromIRing(const int& iRing)
{
  if( (iRing >= 1) && (iRing <= 85) )
    return ( 0.0174*(iRing-0.5) );
  if( (iRing >= -85) && (iRing <= -1) )
    return ( 0.0174*(iRing+0.5) );
  if( iRing >= 86 )
    return ( 1.47135 + 0.0281398*(iRing-86) - 0.00059926*pow(iRing-86,2) + 2.23632e-05*pow(iRing-86,3) );
  if( iRing <= -86 )
    return -1 * ( 1.47135 + 0.0281398*(fabs(iRing)-86) - 0.00059926*pow(fabs(iRing)-86,2) + 2.23632e-05*pow(fabs(iRing)-86,3) );
  return -1.;
}
