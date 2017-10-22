import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
hltHM.HLTPaths = [
                #'HLT_PAPixelTracks_Multiplicity100_v*',
                #'HLT_PAPixelTracks_Multiplicity130_v*',
                #'HLT_PAPixelTracks_Multiplicity160_v*'
                #'HLT_PAPixelTracks_Multiplicity190_v',
                #'HLT_PAPixelTracks_Multiplicity220_v'
                #For 8 TeV
                'HLT_PAFullTracks_Multiplicity185*_v*'
            ]

hltHM.andOr = cms.bool(True)
hltHM.throw = cms.bool(False)

v0CorrelationMC = cms.EDAnalyzer('V0CorrelationMC',
        trkSrc         = cms.InputTag('generalTracks'),
        vertexCollName = cms.InputTag('offlinePrimaryVertices'),
        gnCollection = cms.InputTag('genParticles'),
        multHigh     = cms.untracked.double(250),
        multLow      = cms.untracked.double(185),
        etaMax_ass   = cms.untracked.double(2.4),
        etaMin_ass   = cms.untracked.double(-2.4),
        etaMax_trg   = cms.untracked.double(2.4),
        etaMin_trg   = cms.untracked.double(-2.4),
        ptMin_ass    = cms.untracked.double(0.3),
        ptMax_ass    = cms.untracked.double(3.0),
        bkgnum       = cms.untracked.int32(20),
        rapMin       = cms.untracked.double(0),
        rapMax       = cms.untracked.double(0),
        doKs         = cms.untracked.bool(True),
        doLa         = cms.untracked.bool(True),
        doXi         = cms.untracked.bool(True),
        doRap        = cms.untracked.bool(False),
        doReco = cms.untracked.bool(False)
        )
