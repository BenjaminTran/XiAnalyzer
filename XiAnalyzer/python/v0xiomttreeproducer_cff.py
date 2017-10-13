import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.v0xiomttreeproducer_cfi import *

KsTreeProducer = TreeProducer.clone(
        v0IDName = cms.string('Kshort'),
        multLow = cms.double(185),
        multHigh    = cms.double(250),
        misIDMassCut = cms.double(0.010)
        )

LaTreeProducer = TreeProducer.clone(
        v0IDName = cms.string('Lambda'),
        multLow = cms.double(185),
        multHigh      = cms.double(250),
        misIDMassCut = cms.double(0.020)
        )

XiTreeProducer = TreeProducer.clone(
        v0CollName = cms.string('generalCascadeCandidatesNew'),
        v0IDName = cms.string('Xi'),
        multLow = cms.double(185),
        multHigh    = cms.double(250)
        )

OmTreeProducer = TreeProducer.clone(
        v0CollName = cms.string('generalCascadeCandidatesNew'),
        v0IDName = cms.string('Omega'),
        multLow = cms.double(185),
        multHigh      = cms.double(250),
        misIDMassCut = cms.double(0.015)
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
        multLow = cms.double(0),
        multHigh = cms.double(999999),
        useCentrality = cms.bool(True)
        )

LaTreeProducerRapidityPbPb = LaTreeProducerRapidity.clone(
        multLow = cms.double(0),
        multHigh = cms.double(999999),
        useCentrality = cms.bool(True)
        )

XiTreeProducerRapidityPbPb = XiTreeProducerRapidity.clone(
        multLow = cms.double(0),
        multHigh = cms.double(999999),
        useCentrality = cms.bool(True)
        )

OmTreeProducerRapidityPbPb = OmTreeProducerRapidity.clone(
        multLow = cms.double(0),
        multHigh = cms.double(999999),
        useCentrality = cms.bool(True)
        )
