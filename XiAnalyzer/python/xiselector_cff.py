import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.xiselector_cfi import *
selectV0CandidatesLowXi = selectV0CandidatesLow.clone()
selectV0CandidatesLowOmega = selectV0CandidatesLow.clone(
        v0IDName           = cms.string("Omega"),
        cas3DIpSigValue     = cms.double(3.0),
        casBat3DIpSigValue   = cms.double(4), #Kaon
        VTrkPi3DIpSigValue = cms.double(3),
        VTrkP3DIpSigValue  = cms.double(2),
        casFlightSigValue   = cms.double(2),
        distanceSigValue   = cms.double(10)
        )

selectV0CandidatesLowXiRapidity = selectV0CandidatesLow.clone(
    doRap = cms.bool(True)
        )

selectV0CandidatesLowOmegaRapidity = selectV0CandidatesLowOmega.clone(
    doRap = cms.bool(True)
        )

selectV0CandidatesXiRapidityLoose = selectV0CandidatesLowXiRapidity.clone(
    cas3DIpSigValue        = cms.double(3.0),
    casBat3DIpSigValue      = cms.double(4.5),
    VTrkPi3DIpSigValue    = cms.double(3.5),
    VTrkP3DIpSigValue     = cms.double(2.5),
    casFlightSigValue      = cms.double(2.5),
    distanceSigValue      = cms.double(10.0)
        )

selectV0CandidatesXiRapidityTight = selectV0CandidatesLowXiRapidity.clone(
    cas3DIpSigValue        = cms.double(2.0),
    casBat3DIpSigValue      = cms.double(5.5),
    VTrkPi3DIpSigValue    = cms.double(4.5),
    VTrkP3DIpSigValue     = cms.double(3.5),
    casFlightSigValue      = cms.double(3.5),
    distanceSigValue      = cms.double(14.0)
        )
