import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.masspt_cfi import *

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
#hltHM = hltHM.clone();
