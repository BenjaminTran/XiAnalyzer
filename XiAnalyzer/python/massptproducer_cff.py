import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.massptproducer_cfi import *

#XiMassPt = MassPt.clone(
#    ks = cms.untracked.bool(False),
#    la = cms.untracked.bool(False),
#    om = cms.untracked.bool(False)
#)
#
#OmMassPt = MassPt.clone(
#        ks = cms.untracked.bool(False),
#        la = cms.untracked.bool(False),
#        xi = cms.untracked.bool(False)
#        )
#
#KslaMassPt = MassPt.clone(
#    xi = cms.untracked.bool(False)
#)
#
#KsMassPt = MassPt.clone(
#    xi = cms.untracked.bool(False),
#    la = cms.untracked.bool(False)
#)
#
#LaMassPt = MassPt.clone(
#    xi = cms.untracked.bool(False),
#    ks = cms.untracked.bool(False)
#)

MassPtRapidity = MassPt.clone(
        multHigh = cms.double(250),
        gnCollection = cms.InputTag('genParticles'),
        ksCollection = cms.InputTag('selectV0CandidatesNewkshortRapidity:Kshort'),
        laCollection = cms.InputTag('selectV0CandidatesNewlambdaRapidity:Lambda'),
        xiCollection = cms.InputTag('selectV0CandidatesLowXiRapidity:Xi'),
        omCollection = cms.InputTag('selectOmegaCandidatesNewRapidity:Omega')
)

MassPtRapidityMC = MassPt.clone(
        MC = cms.untracked.bool(True),
        om = cms.untracked.bool(True),
        multLow = cms.double(0),
        multHigh = cms.double(999999),
        gnCollection = cms.InputTag('genParticles'),
        ksCollection = cms.InputTag('selectV0CandidatesNewkshortRapidity:Kshort'),
        laCollection = cms.InputTag('selectV0CandidatesNewlambdaRapidity:Lambda'),
        xiCollection = cms.InputTag('selectV0CandidatesLowXiRapidity:Xi')
)
#hltHM = hltHM.clone();
