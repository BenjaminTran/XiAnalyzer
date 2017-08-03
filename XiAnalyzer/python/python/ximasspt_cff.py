import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.ximasspt_cfi import *

MassPt = xiMassPt.clone()

XiMassPt = xiMassPt.clone(
    ks = cms.untracked.bool(False),
    la = cms.untracked.bool(False)
)

KslaMassPt = xiMassPt.clone(
    xi = cms.untracked.bool(False)
)

KsMassPt = xiMassPt.clone(
    xi = cms.untracked.bool(False),
    la = cms.untracked.bool(False)
)

LaMassPt = xiMassPt.clone(
    xi = cms.untracked.bool(False),
    ks = cms.untracked.bool(False)
)
#hltHM = hltHM.clone();
