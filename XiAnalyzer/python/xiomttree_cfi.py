import FWCore.ParameterSet.Config as cms

xiTree = cms.EDAnalyzer('XiOmTTree',
     vertexCollName  = cms.InputTag('offlinePrimaryVertices'),
     trkSrc = cms.InputTag('generalTracks'),
     #v0CollName      = cms.string('generalCascadeCandidatesNew'),
     v0CollName      = cms.string('generalV0CandidatesNew'), #For Test dataset
     v0IDName = cms.string('Xi')
)
