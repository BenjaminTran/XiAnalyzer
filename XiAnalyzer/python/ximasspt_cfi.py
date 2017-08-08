import FWCore.ParameterSet.Config as cms

xiMassPt = cms.EDAnalyzer('XiMassPt',
        ksCollection   = cms.InputTag('selectV0CandidatesNewkshort:Kshort'),
        laCollection   = cms.InputTag('selectV0CandidatesNewlambdatight:Lambda'),
        trkSrc         = cms.InputTag('generalTracks'),
        vertexCollName = cms.InputTag('offlinePrimaryVertices'),
        xiCollection   = cms.InputTag('selectV0CandidatesLowXi:Xi'),
        multHigh       = cms.double(250),
        multLow        = cms.double(185),
        zVtxHigh       = cms.double(15.0),
        zVtxLow        = cms.double(-15.0),
        ks             = cms.untracked.bool(True),
        la             = cms.untracked.bool(True),
        xi             = cms.untracked.bool(True)
)
