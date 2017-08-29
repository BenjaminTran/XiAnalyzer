import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.xiselector_cfi import *
selectV0CandidatesLowXi = selectV0CandidatesLow.clone()

selectV0CandidatesLowXiRapidity = selectV0CandidatesLow.clone(
    doRap = cms.bool(True)
        )

selectV0CandidatesXiRapidityLoose = selectV0CandidatesLowXiRapidity.clone(
    xi3DIpSigValue        = cms.double(3.0),
    xiPi3DIpSigValue      = cms.double(4.5),
    VTrkPi3DIpSigValue    = cms.double(3.5),
    VTrkP3DIpSigValue     = cms.double(2.5),
    xiFlightSigValue      = cms.double(2.5),
    distanceSigValue      = cms.double(10.0)
        )

selectV0CandidatesXiRapidityTight = selectV0CandidatesLowXiRapidity.clone(
    xi3DIpSigValue        = cms.double(2.0),
    xiPi3DIpSigValue      = cms.double(5.5),
    VTrkPi3DIpSigValue    = cms.double(4.5),
    VTrkP3DIpSigValue     = cms.double(3.5),
    xiFlightSigValue      = cms.double(3.5),
    distanceSigValue      = cms.double(14.0)
        )
