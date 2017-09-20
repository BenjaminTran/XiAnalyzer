import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.v0selector_cfi import *
selectV0CandidatesNewkshort = selectV0CandidatesNew.clone()
selectV0CandidatesNewlambda = selectV0CandidatesNew.clone(
  v0IDName = cms.string("Lambda"),
)

selectV0CandidatesNewkshortRapidity = selectV0CandidatesNew.clone(
        dorap = cms.bool(True)
)

selectV0CandidatesNewlambdaRapidity = selectV0CandidatesNew.clone(
  v0IDName     = cms.string("Lambda"),
  dorap        = cms.bool(True)
)

selectV0CandidatesNewkshortRapidityMC = selectV0CandidatesNewkshortRapidity.clone(
        v0CollName = cms.string("generalV0Candidates")
)

selectV0CandidatesNewlambdaRapidityMC = selectV0CandidatesNewlambdaRapidity.clone(
        v0CollName = cms.string("generalV0Candidates")
)

selectV0CandidatesNewlambdatight = selectV0CandidatesNewlambda.clone(
  dxySigCut1      = cms.double(1.25),
  dxySigCut2      = cms.double(1.25),
  dzSigCut1       = cms.double(1.25),
  dzSigCut2       = cms.double(1.25),
  vtxChi2Cut      = cms.double(10000.0),
  cosThetaCut     = cms.double(0.9995),
  decayLSigCut    = cms.double(7.0),
)

selectV0CandidatesNewlambdaloose = selectV0CandidatesNewlambda.clone(
  dxySigCut1      = cms.double(1.),
  dxySigCut2      = cms.double(1.),
  dzSigCut1       = cms.double(1.),
  dzSigCut2       = cms.double(1.),
  vtxChi2Cut      = cms.double(10000.0),
  cosThetaCut     = cms.double(0.995),
  decayLSigCut    = cms.double(3.0),
)

selectV0CandidatesNewlambdalooseRapidity = selectV0CandidatesNewlambdaloose.clone(
        dorap        = cms.bool(True),
        decayLSigCut = cms.double(2.5),
        cosThetaCut  = cms.double(0.99)
        )

selectV0CandidatesNewlambdatightRapidity = selectV0CandidatesNewlambdatight.clone(
        dorap        = cms.bool(True),
        decayLSigCut = cms.double(7.5),
        cosThetaCut  = cms.double(0.9999)
        )

selectV0CandidatesNewkshorttightRapidity = selectV0CandidatesNewkshortRapidity.clone(
        dxySigCut1      = cms.double(1.25),
        dxySigCut2      = cms.double(1.25),
        dzSigCut1       = cms.double(1.25),
        dzSigCut2       = cms.double(1.25),
        decayLSigCut = cms.double(7.5),
        cosThetaCut  = cms.double(0.9999)
        )

selectV0CandidatesNewkshortlooseRapidity = selectV0CandidatesNewkshortRapidity.clone(
        decayLSigCut = cms.double(2.5),
        cosThetaCut  = cms.double(0.99)
        )

selectV0CandidatesNewd0 = selectV0CandidatesNew.clone(
  v0IDName = cms.string("D0"),
  ptCut1          = cms.double(0.5),
  ptCut2          = cms.double(0.5),
  nHitCut1        = cms.int32(10),
  nHitCut2        = cms.int32(10),
  dxySigCut1      = cms.double(0.5),
  dxySigCut2      = cms.double(0.5),
  dzSigCut1       = cms.double(0.0),
  dzSigCut2       = cms.double(0.0),
  vtxChi2Cut      = cms.double(100000.0),
  cosThetaCut     = cms.double(0.98),
  decayLSigCut    = cms.double(4.0),
  misIDMassCut   = cms.double(0.000),
  misIDMassCutEE = cms.double(0.000)
)

