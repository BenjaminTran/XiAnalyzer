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

xiCorrelation          = cms.EDAnalyzer('XiCorrelation',
        #Input tags
        trkSrc         = cms.InputTag('generalTracks'),
        vertexCollName = cms.InputTag('offlinePrimaryVertices'),
        #casCollection   = cms.InputTag('selectV0CandidatesLowXi:Xi'),
        casCollection   = cms.InputTag('selectV0CandidatesLowXi:Xi'),
        #Parameters
        bkgnum         = cms.double(20.0),
        etaMax_ass     = cms.double(2.4),
        etaMin_ass     = cms.double(-2.4),
        multHigh       = cms.double(250),
        multLow        = cms.double(185),
        peakFactor     = cms.int32(2),
        ptBin          = cms.vdouble(1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0),
        ptMax_ass      = cms.double(3.0),
        ptMin_ass      = cms.double(0.3),
        sideFactor     = cms.int32(3),
        zVtxHigh       = cms.double(15.0),
        zVtxLow        = cms.double(-15.0),
        PtBinNum       = cms.int32(7),
        xiMassMean     = cms.vdouble(1.32296, 1.32254, 1.3225, 1.32246, 1.32242, 1.32226, 1.32214),
        xiMassSigma    = cms.vdouble(0.00558245, 0.00520491, 0.00518877, 0.00514213, 0.00512142, 0.00512142, 0.00490732),
        doRap          = cms.bool(False)
        )
