import FWCore.ParameterSet.Config as cms

selectV0CandidatesLow  = cms.EDProducer('XiSelector',
    vertexCollName     = cms.InputTag('offlinePrimaryVertices'),
    v0CollName         = cms.string("generalCascadeCandidatesNew"),
    v0IDName           = cms.string("Xi"),
    etaCutMin          = cms.double(-2.4),
    etaCutMax          = cms.double(2.4),
    zVertexLow         = cms.double(-15.0),
    zVertexHigh        = cms.double(15.0),
    ptCut1             = cms.double(0.0),
    ptCut2             = cms.double(0.0),
    nHitCut1           = cms.int32(3),
    xi3DIpSigValue     = cms.double(2.5),
    xiPi3DIpSigValue   = cms.double(5),
    VTrkPi3DIpSigValue = cms.double(4),
    VTrkP3DIpSigValue  = cms.double(3),
    xiFlightSigValue   = cms.double(3),
    distanceSigValue   = cms.double(12),
    vtxChi2Cut         = cms.double(10000.0),
    cosThetaCut        = cms.double(0.999),
    misIDMassCut       = cms.double(0.010),
    misIDMassCutEE     = cms.double(0.015)
    )
