import FWCore.ParameterSet.Config as cms

selectOmegaCandidatesNew = cms.EDAnalyzer('OmegaSelector',
        vertexCollName   = cms.InputTag('offlinePrimaryVertices'),
        #v0CollName = cms.string("generalCascadeCandidatesNew"),
        v0CollName       = cms.string("generalV0CandidatesNew"),
        v0IDName         = cms.string("Omega"),
        etaCutMin        = cms.double(-2.4),
        etaCutMax        = cms.double(2.4),
        )
