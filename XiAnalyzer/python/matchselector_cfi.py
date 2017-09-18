import FWCore.ParameterSet.Config as cms

selectV0CandidatesNew = cms.EDProducer('V0Selector',
        gnCollection = cms.InputTag('genParticles'),
        vertexCollName      = cms.InputTag('offlinePrimaryVertices'),
        v0CollName          = cms.string("generalV0CandidatesNew"),
        v0IDName            = cms.string("Kshort"),
        etaCutMin           = cms.double(-2.4),
        etaCutMax           = cms.double(2.4),
        ptCut1              = cms.double(0.0),
        ptCut2              = cms.double(0.0),
        dorap               = cms.bool(False),
        rapMax              = cms.double(1.0),
        rapMin              = cms.double(-1.0)
)
