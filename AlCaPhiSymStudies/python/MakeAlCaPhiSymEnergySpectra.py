import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.ioFiles = cms.PSet(

    outputFile = cms.string('/afs/cern.ch/work/f/fedmante/PhiSym_HLTStudies/CMSSW_7_4_0_pre9/src/HLTStudies/AlCaPhiSymStudies/output.root'),

    inputFiles = cms.vstring(
        "root://eoscms.cern.ch//store/group/dpg_ecal/alca_ecalcalib/bmarzocc/AlCaPhiSym_RECO/PU20_25ns_eb6_ee9/outputALCAPHISYM_RECO_27.root",
    )
)
