import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.v0xiomttreeproducer_cfi import *

XiTreeProducerRapidity = TreeProducer.Clone(
        v0CollName = cms.string('generalCascadeCandidatesNew'),
        v0IDName = cms.string('Xi'),
        doRap = cms.bool(True),
        doXi = cms.bool(True),
        multmax = cms.double(250)
        )

KsTreeProducerRapidity = TreeProducer.Clone(
        v0CollName = cms.string('generalV0CandidatesNew'),
        v0IDName = cms.string('Kshort'),
        doRap = cms.bool(True),
        doKs = cms.bool(True),
        multmax = cms.double(250)
        )

LaTreeProducerRapidity = TreeProducer.Clone(
        v0CollName = cms.string('generalV0CandidatesNew'),
        v0IDName = cms.string('Lambda'),
        doRap = cms.bool(True),
        doLa = cms.bool(True),
        multmax = cms.double(250)
        )
