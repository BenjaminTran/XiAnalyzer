import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.omegaselector_cfi import *
selectOmCandidatesNew = selectOmCandidates.clone()

selectOmCandidatesLowRapidity = selectOmCandidates.clone(
    dorap = cms.bool(True)
        )

selectOmCandidatesLoose = selectOmCandidates.clone(
    xi3DIpSigValue        = cms.double(2.5),
    xiPi3DIpSigValue      = cms.double(5.0),
    VTrkPi3DIpSigValue    = cms.double(4.0),
    VTrkP3DIpSigValue     = cms.double(3.0),
    xiFlightSigValue      = cms.double(3.0),
    distanceSigValue      = cms.double(12.0)
        )
