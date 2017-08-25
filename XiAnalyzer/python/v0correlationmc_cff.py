import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.v0correlationmc_cfi import *

v0CorrelationMCRapidity = v0CorrelationMC.clone(
        gnCollection = cms.InputTag('genParticles'),
        trkSrc         = cms.InputTag('generalTracks'),
        vertexCollName = cms.InputTag('offlinePrimaryVertices'),
        doRap  = cms.untracked.bool(True),
        rapMin = cms.untracked.double(-1.0),
        rapMax = cms.untracked.double(1.0),
        ptcut_ks       = cms.untracked.vdouble(0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0,15.0,20.0,25.0,30.0),
        ptcut_la       = cms.untracked.vdouble(0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0,15.0,20.0,25.0,30.0),
        multHigh = cms.untracked.double(9999999),
        multLow = cms.untracked.double(-1)
        )
#hltHM = hltHM.clone()
