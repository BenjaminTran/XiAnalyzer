import FWCore.ParameterSet.Config as cms

selectGenCandidates = cms.EDProducer('GenSelector',
        gnCollection = cms.InputTag('genParticles'),
        v0IDName     = cms.string("Kshort")
)
