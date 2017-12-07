import FWCore.ParameterSet.Config as cms

HighMultFilter = cms.EDFilter('HighMultFilter',
        genSrc = cms.InputTag('genParticles'),
        vertexSrc = cms.InputTag('offlinePrimaryVertices'),
        trackSrc = cms.InputTag('generalTracks'),
        multMax = cms.int32(250),
        multMin = cms.int32(185),
        etaMax = cms.double(2.4),
        etaMin = cms.double(-2.4),
        doGenParticle = cms.bool(False)
        )
