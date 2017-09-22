import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.v0xiomttreeproducer_cfi import *

KsTreeProducer = TreeProducer.clone(
        v0CollName = cms.string('generalV0CandidatesNew'),
        v0IDName   = cms.string('Kshort'),
        multmin = cms.double(185),
        multmax    = cms.double(250)
        )

LaTreeProducer = TreeProducer.clone(
        v0CollName   = cms.string('generalV0CandidatesNew'),
        v0IDName     = cms.string('Lambda'),
        multmin = cms.double(185),
        multmax      = cms.double(250),
        misIDMassCut = cms.double(0.020)
        )

XiTreeProducer = TreeProducer.clone(
        v0CollName = cms.string('generalCascadeCandidatesNew'),
        v0IDName   = cms.string('Xi'),
        multmin = cms.double(185),
        multmax    = cms.double(250)
        )

OmTreeProducer = TreeProducer.clone(
        v0CollName   = cms.string('generalCascadeCandidatesNew'),
        v0IDName     = cms.string('Omega'),
        multmin = cms.double(185),
        multmax      = cms.double(250),
        misIDMassCut = cms.double(0.020)
        )

OmTreeProducerRapidity = OmTreeProducer.clone(
        doRap = cms.bool(True)
        )

XiTreeProducerRapidity = XiTreeProducer.clone(
        doRap = cms.bool(True)
        )

KsTreeProducerRapidity = KsTreeProducer.clone(
        doRap = cms.bool(True)
        )

LaTreeProducerRapidity = LaTreeProducer.clone(
        doRap = cms.bool(True)
        )

KsTreeProducerRapidityPbPb = KsTreeProducerRapidity.clone(
        multmin = cms.double(0),
        useCentrality = cms.bool(True)
        )

LaTreeProducerRapidityPbPb = LaTreeProducerRapidity.clone(
        multmin = cms.double(0),
        useCentrality = cms.bool(True)
        )

XiTreeProducerRapidityPbPb = XiTreeProducerRapidity.clone(
        multmin = cms.double(0),
        useCentrality = cms.bool(True)
        )

OmTreeProducerRapidityPbPb = OmTreeProducerRapidity.clone(
        multmin = cms.double(0),
        useCentrality = cms.bool(True)
        )
