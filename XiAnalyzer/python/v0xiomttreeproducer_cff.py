import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.v0xiomttreeproducer_cfi import *

KsTreeProducer = TreeProducer.clone(
        multmin = cms.double(185),
        multmax    = cms.double(250),
        doKs = cms.bool(True)
        )

LaTreeProducer = TreeProducer.clone(
        multmin = cms.double(185),
        multmax      = cms.double(250),
        doLa = cms.bool(True)
        )

XiTreeProducer = TreeProducer.clone(
        multmin = cms.double(185),
        multmax    = cms.double(250),
        doXi = cms.bool(True)
        )

OmTreeProducer = TreeProducer.clone(
        multmin = cms.double(185),
        multmax      = cms.double(250),
        misIDMassCut = cms.double(0.020),
        doOm = cms.bool(True)
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
