import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.v0selector_cfi import *
selectV0CandidatesNewkshort = selectV0CandidatesNew.clone()
selectV0CandidatesNewlambda = selectV0CandidatesNew.clone(
  v0IDName = cms.string("Lambda"),
  misIDMassCut   = cms.double(0.020)
)

selectV0CandidatesNewkshortRapidity = selectV0CandidatesNew.clone(
        dorap = cms.bool(True)
)

selectV0CandidatesNewlambdaRapidity = selectV0CandidatesNew.clone(
  v0IDName     = cms.string("Lambda"),
  misIDMassCut = cms.double(0.020),
  dorap        = cms.bool(True)
)

selectV0CandidatesNewkshortRapidityPbPb = selectV0CandidatesNew.clone(
        dorap = cms.bool(True),
)

selectV0CandidatesNewlambdaRapidityPbPb = selectV0CandidatesNewlambda.clone(
        dorap        = cms.bool(True),
        cosThetaCut0 = cms.double(0.9999), #First bin
        cosThetaCut = cms.double(0.9998) #Other bins
)

#### SYSTEMATICS

selectV0CandidatesNewlambdatight = selectV0CandidatesNewlambda.clone(
  dxySigCut1      = cms.double(1.25),
  dxySigCut2      = cms.double(1.25),
  dzSigCut1       = cms.double(1.25),
  dzSigCut2       = cms.double(1.25),
  vtxChi2Cut      = cms.double(10000.0),
  cosThetaCut     = cms.double(0.9999),
  decayLSigCut    = cms.double(7.5),
)

selectV0CandidatesNewlambdaloose = selectV0CandidatesNewlambda.clone(
  dxySigCut1      = cms.double(1.),
  dxySigCut2      = cms.double(1.),
  dzSigCut1       = cms.double(1.),
  dzSigCut2       = cms.double(1.),
  vtxChi2Cut      = cms.double(10000.0),
  cosThetaCut     = cms.double(0.99),
  decayLSigCut    = cms.double(2.5),
)

selectV0CandidatesNewlambdalooseRapidity = selectV0CandidatesNewlambdaloose.clone(
        dorap        = cms.bool(True),
        )

selectV0CandidatesNewlambdatightRapidity = selectV0CandidatesNewlambdatight.clone(
        dorap        = cms.bool(True),
        )

selectV0CandidatesNewkshorttightRapidity = selectV0CandidatesNewlambdaloose.clone(
        v0IDName     = cms.string("Kshort"),
        dorap = cms.bool(True)
        )

selectV0CandidatesNewkshortlooseRapidity = selectV0CandidatesNewlambdatight.clone(
        v0IDName     = cms.string("Kshort"),
        dorap = cms.bool(True)
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

