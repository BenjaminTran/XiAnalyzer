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
        trkSrc         = cms.InputTag('generalTracks'),
        xiCollection   = cms.InputTag('selectV0CandidatesLowXi:Xi'),
        vertexCollName = cms.InputTag('offlinePrimaryVertices'),
        zVtxHigh       = cms.double(15.0),
        zVtxLow        = cms.double(-15.0),
        ptMax_ass      = cms.double(3.0),
        ptMin_ass      = cms.double(0.3),
        etaMax_ass     = cms.double(2.4),
        etaMin_ass     = cms.double(-2.4),
        multHigh       = cms.int32(220),
        multLow        = cms.int32(185),
        bkgnum         = cms.double(20.0),
        ptBin          = cms.vdouble(1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0),
        # Root fit values
        #xiMassMean     = cms.vdouble(1.32245, 1.3226, 1.32259, 1.32247, 1.32239, 1.32224, 1.32212),
        #xiMassSigma    = cms.vdouble(0.00344097, 0.00544116, 0.00424282, 0.00539586, 0.00588911, 0.00537808, 0.00281928),
        # RooFit values
        xiMassMean     = cms.vdouble(1.32296, 1.32254, 1.3225, 1.32246, 1.32242, 1.32226, 1.32214),
        xiMassSigma    = cms.vdouble(0.00558245, 0.00520491, 0.00518877, 0.00514213, 0.00512142, 0.00512142, 0.00490732),
        peakFactor     = cms.int32(2),
        sideFactor     = cms.int32(3),
        PtBinNum       = cms.int32(7)
        )
