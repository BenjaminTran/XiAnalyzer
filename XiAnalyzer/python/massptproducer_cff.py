import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.massptproducer_cfi import *

MassPtRapidity = MassPt.clone(
        multHigh = cms.double(250),
        gnCollection = cms.InputTag('genParticles'),
        ksCollection = cms.InputTag('selectV0CandidatesNewkshortRapidity:Kshort'),
        laCollection = cms.InputTag('selectV0CandidatesNewlambdaRapidity:Lambda'),
        xiCollection = cms.InputTag('selectV0CandidatesLowXiRapidity:Xi'),
        omCollection = cms.InputTag('selectOmegaCandidatesNewRapidity:Omega')
)

MassPtRapidityMB = MassPtRapidity.clone(
        multHigh = cms.double(35),
        multLow = cms.double(0)
        )

MassPtRapidityMC = MassPt.clone(
        MC = cms.untracked.bool(True),
        multLow = cms.double(0),
        multHigh = cms.double(999999),
        gnCollection = cms.InputTag('genParticles'),
        ksCollection = cms.InputTag('selectV0CandidatesNewkshortRapidity:Kshort'),
        laCollection = cms.InputTag('selectV0CandidatesNewlambdaRapidity:Lambda'),
        xiCollection = cms.InputTag('selectV0CandidatesLowXiRapidity:Xi'),
        omCollection = cms.InputTag('selectOmegaCandidatesNewRapidity:Omega')
)
#hltHM = hltHM.clone();
