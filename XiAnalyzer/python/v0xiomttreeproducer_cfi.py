import FWCore.ParameterSet.Config as cms

TreeProducer = cms.EDAnalyzer('V0XiOmTTreeProducer',
     vertexCollName  = cms.InputTag('offlinePrimaryVertices'),
     trkSrc = cms.InputTag('generalTracks'),
     multmin = cms.double(185),
     multmax = cms.double(220),
     doRap = cms.bool(False),
     doXi = cms.bool(False),
     doKs = cms.bool(False),
     doLa = cms.bool(False),
     doOm = cms.bool(False),
     rapMin = cms.double(-1.0),
     rapMax = cms.double(1.0)
)
