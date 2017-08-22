import FWCore.ParameterSet.Config as cms

MassPt = cms.EDAnalyzer('MassPtProducer',
        ksCollection   = cms.InputTag('selectV0CandidatesNewkshort:Kshort'),
        laCollection   = cms.InputTag('selectV0CandidatesNewlambda:Lambda'),
        trkSrc         = cms.InputTag('generalTracks'),
        vertexCollName = cms.InputTag('offlinePrimaryVertices'),
        xiCollection   = cms.InputTag('selectV0CandidatesLowXi:Xi'),
        omCollection   = cms.InputTag('selectOmCandidatesNew:Omega'),
        multHigh       = cms.double(220),
        multLow        = cms.double(185),
        zVtxHigh       = cms.double(15.0),
        zVtxLow        = cms.double(-15.0),
        rapMin         = cms.double(-1.0),
        rapMax         = cms.double(1.0),
        MC             = cms.untracked.bool(False),
        ks             = cms.untracked.bool(True),
        la             = cms.untracked.bool(True),
        xi             = cms.untracked.bool(True),
        om             = cms.untracked.bool(True)
)
