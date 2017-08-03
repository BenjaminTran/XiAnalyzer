import FWCore.ParameterSet.Config as cms

#import HLTrigger.HLTfilters.hltHighLevel_cfi
#hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
#hltHM.HLTPaths = [
#                'HLT_PAPixelTracks_Multiplicity100_v*',
#                'HLT_PAPixelTracks_Multiplicity130_v*',
#                'HLT_PAPixelTracks_Multiplicity160_v*'
#                #'HLT_PAPixelTracks_Multiplicity190_v',
#                #'HLT_PAPixelTracks_Multiplicity220_v'
#            ]

xiMassPt = cms.EDAnalyzer('XiMassPt',
        trkSrc = cms.InputTag('generalTracks'),
        xiCollection   = cms.InputTag('selectV0CandidatesLowXi:Xi'),
        laCollection = cms.InputTag('selectV0CandidatesNewlambdatight:Lambda'),
        ksCollection = cms.InputTag('selectV0CandidatesNewkshort:Kshort'),
        vertexCollName = cms.InputTag('offlinePrimaryVertices'),
        zVtxHigh       = cms.double(15.0),
        zVtxLow        = cms.double(-15.0),
        multHigh       = cms.double(220),
        multLow        = cms.double(185),
        #Leave these as true, select variations from cff cloned processes
        xi = cms.untracked.bool(True),
        ks = cms.untracked.bool(True),
        la = cms.untracked.bool(True)

)
