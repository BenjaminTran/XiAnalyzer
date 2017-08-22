import FWCore.ParameterSet.Config as cms

selectOmCandidatesNew  = cms.EDProducer('OmSelector',
    vertexCollName     = cms.InputTag('offlinePrimaryVertices'),
    v0CollName         = cms.string("generalCascadeCandidatesNew"),
    v0IDName           = cms.string("Omega"),
    etaCutMin          = cms.double(-2.4),
    etaCutMax          = cms.double(2.4),
    zVertexLow         = cms.double(-15.0),
    zVertexHigh        = cms.double(15.0),
    ptCut1             = cms.double(0.0),
    ptCut2             = cms.double(0.0),
    nHitCut1           = cms.int32(3),
    xi3DIpSigValue     = cms.double(3.0),
    xiPi3DIpSigValue   = cms.double(4), #Kaon
    VTrkPi3DIpSigValue = cms.double(3),
    VTrkP3DIpSigValue  = cms.double(2),
    xiFlightSigValue   = cms.double(2),
    distanceSigValue   = cms.double(10),
    dorap              = cms.bool(False),
    rapMax             = cms.double(1.0),
    rapMin             = cms.double(-1.0)
    )
