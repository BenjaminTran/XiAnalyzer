import FWCore.ParameterSet.Config as cms

MatchCandidates = cms.EDProducer('MatchSelector',
        #gnCollection = cms.InputTag('genParticles'),
        vertexCollName      = cms.InputTag('offlinePrimaryVertices'),
        v0CollName          = cms.string("generalV0CandidatesNew"),
        v0IDName            = cms.string("Kshort"),
        gnV0Collection = cms.InputTag('selectGenCandidates:Kshort'),
        etaCutMin           = cms.double(-2.4),
        etaCutMax           = cms.double(2.4),
        ptCut1              = cms.double(0.0),
        ptCut2              = cms.double(0.0),
        dorap               = cms.bool(True),
        rapMax              = cms.double(1.0),
        rapMin              = cms.double(-1.0)
)
