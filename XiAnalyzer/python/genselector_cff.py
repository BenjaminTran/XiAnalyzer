import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.genselector_cfi import *

selectGenCandidatesKshort = selectGenCandidates.clone()

selectGenCandidatesLambda = selectGenCandidates.clone(
        v0IDName = cms.string('Lambda')
        )
