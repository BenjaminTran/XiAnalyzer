import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.v0correlationmc_cfi import *

v0CorrelationMCRapidity = v0CorrelationMC.clone(
        gnCollection   = cms.InputTag('genParticles'),
        trkSrc         = cms.InputTag('generalTracks'),
        vertexCollName = cms.InputTag('offlinePrimaryVertices'),
        doRap          = cms.untracked.bool(True),
        doKs           = cms.untracked.bool(True),
        doLa           = cms.untracked.bool(True),
        doXi           = cms.untracked.bool(False),
        doReco         = cms.untracked.bool(False),
        rapMin         = cms.untracked.double(-1.0),
        rapMax         = cms.untracked.double(1.0),
        etaMax_trg     = cms.untracked.double(999),
        etaMin_trg     = cms.untracked.double(-999),
        ptcut_ks       = cms.untracked.vdouble(0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0,15.0,20.0,25.0,30.0),
        ptcut_la       = cms.untracked.vdouble(0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0,15.0,20.0,25.0,30.0),
        ptcut_xi       = cms.untracked.vdouble(1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10, 20),
        multHigh       = cms.untracked.double(999999),
        multLow        = cms.untracked.double(-1)
        )
