import FWCore.ParameterSet.Config as cms

makeAlCaPhiSymSpectra = cms.EDAnalyzer("MakeAlCaPhiSymSpectra",
    recHitCollection_EB                   = cms.InputTag("ecalRechitMultifit","EcalMultifitRecHitsEB"),
    recHitCollection_EE                   = cms.InputTag("ecalRechitMultifit","EcalMultifitRecHitsEE"),
)
