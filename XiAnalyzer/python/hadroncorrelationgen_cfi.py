import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.hadroncorrelationgen_cfi import *

HadronCorrelation = cms.EDAnalyzer('HadronCorrelationGen',
        vertexCollName = cms.InputTag('offlinePrimaryVertices'),
        trkSrc = cms.InputTag('generalTracks'),
        gnCollection = cms.InputTag('genParticles'),
        etaMin_trg = cms.double(-2.4),
        etaMax_trg = cms.double(2.4),
        etaMin_ass = cms.double(-2.4),
        etaMax_ass = cms.double(2.4),
        ptMin_ass = cms.double(0.3),
        ptMax_ass = cms.double(3.0),
        bkgFactor = cms.double(20),
        rapMax = cms.double(1.0),
        rapMin = cms.double(-1.0),
        multMax = cms.double(9999999),
        multMin = cms.double(-1),
        ptcut = cms.untracked.vdouble(0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0,15.0,20.0,25.0,30.0)
        )
