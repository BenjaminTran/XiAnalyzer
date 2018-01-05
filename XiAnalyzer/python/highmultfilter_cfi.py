import FWCore.ParameterSet.Config as cms

HighMultFilter = cms.EDFilter('HighMultFilter',
        genSrc = cms.InputTag('genParticles'),
        vertexSrc = cms.InputTag('offlinePrimaryVertices'),
        trackSrc = cms.InputTag('generalTracks'),
        multMax = cms.double(35),
        multMin = cms.double(0),
        etaMax = cms.double(2.4),
        etaMin = cms.double(-2.4),
        doGenParticle = cms.bool(False)
        )
