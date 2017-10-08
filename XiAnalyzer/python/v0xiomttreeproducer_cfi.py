import FWCore.ParameterSet.Config as cms

TreeProducer = cms.EDAnalyzer('V0XiOmTTreeProducer',
     vertexCollName = cms.InputTag('offlinePrimaryVertices'),
     v0CollName   = cms.string('generalV0CandidatesNew'),
     trkSrc         = cms.InputTag('generalTracks'),
     towerSrc       = cms.InputTag('towerMaker'),
     multLow        = cms.double(0),
     multHigh        = cms.double(220),
     doRap          = cms.bool(False),
     useCentrality  = cms.bool(False),
     rapMin         = cms.double(-1.0),
     rapMax         = cms.double(1.0),
     misIDMassCut = cms.double(0.020),
     misIDMassCutEE = cms.double(0.015),
     cent_bin_low   = cms.int32(60),
     cent_bin_high  = cms.int32(100),
     nHitCut1       = cms.int32(3),
     nHitCut2       = cms.int32(3),
     ptCut1         = cms.double(0),
     ptCut2         = cms.double(0)
)
