import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.massptproducer_cfi import *

XiMassPt = MassPt.clone(
    ks = cms.untracked.bool(False),
    la = cms.untracked.bool(False)
)

KslaMassPt = MassPt.clone(
    xi = cms.untracked.bool(False)
)

KsMassPt = MassPt.clone(
    xi = cms.untracked.bool(False),
    la = cms.untracked.bool(False)
)

LaMassPt = MassPt.clone(
    xi = cms.untracked.bool(False),
    ks = cms.untracked.bool(False)
)

MassPtRapidity = MassPt.clone(
        multHigh = cms.double(250),
        ksCollection = cms.InputTag('selectV0CandidatesNewkshortRapidity:Kshort'),
        laCollection = cms.InputTag('selectV0CandidatesNewlambdaRapidity:Lambda'),
        xiCollection = cms.InputTag('selectV0CandidatesLowXiRapidity:Xi')
)
#hltHM = hltHM.clone();
