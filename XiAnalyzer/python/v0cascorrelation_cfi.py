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

v0CasCorrelation          = cms.EDAnalyzer('V0CasCorrelation',
        trkSrc         = cms.InputTag('generalTracks'),
        towerSrc      = cms.InputTag('towerMaker'),
        laCollection   = cms.untracked.InputTag('selectV0CandidatesNewlambda:Lambda'),
        ksCollection   = cms.untracked.InputTag('selectV0CandidatesNewkshort:Kshort'),
        xiCollection = cms.untracked.InputTag('selectV0CandidatesLowXi:Xi'),
        omCollection = cms.untracked.InputTag('selectOmegaCandidatesNew:Omega'),
        vertexCollName = cms.InputTag('offlinePrimaryVertices'),
        zVtxHigh       = cms.untracked.double(15.0),
        zVtxLow        = cms.untracked.double(-15.0),
        multHigh       = cms.untracked.double(220),
        multLow        = cms.untracked.double(185),
        etaMax_ass     = cms.untracked.double(2.4),
        etaMin_ass     = cms.untracked.double(-2.4),
        etaMax_trg     = cms.untracked.double(2.4),
        etaMin_trg     = cms.untracked.double(-2.4),
        ptMin_ass      = cms.untracked.double(0.3),
        ptMax_ass      = cms.untracked.double(3.0),
        bkgnum         = cms.untracked.int32(20),
        ptbin_n        = cms.untracked.int32(13),
        ptbin_n_cas_        = cms.untracked.int32(13),
        peakFactor     = cms.untracked.double(2.0),
        sideFactor     = cms.untracked.double(3.0),
        mis_ks_range   = cms.untracked.double(0.020),
        mis_la_range   = cms.untracked.double(0.010),
        mis_ph_range   = cms.untracked.double(0.015),
        sigma_ks       = cms.untracked.vdouble(0.00782525, 0.00775315 ,0.00742869 ,0.00710546 ,0.00683902 ,0.00664499 ,0.0066186 ,0.00670261 ,0.00680736 ,0.00700286 ,0.00728129 ,0.00742705 ,0.00774735),
        mean_ks        = cms.untracked.vdouble(0.497166 ,0.497335 ,0.497627 ,0.497753 ,0.49775 ,0.497681 ,0.497606 ,0.497546 ,0.497517 ,0.497492 ,0.497507 ,0.497516 ,0.497598),
        sigma_la       = cms.untracked.vdouble(0.00309807 ,0.00302995 ,0.00304411 ,0.00309178 ,0.00295109 ,0.00288246 ,0.00285832 ,0.00289492 ,0.00326894 ,0.00275142),
        mean_la        = cms.untracked.vdouble(1.11636 ,1.11605 ,1.11598 ,1.116 ,1.11597 ,1.1159 ,1.11585 ,1.11581 ,1.11582 ,1.11584),
        mean_xi     = cms.untracked.vdouble(1.32296, 1.32254, 1.3225, 1.32246, 1.32242, 1.32226, 1.32214),
        sigma_xi    = cms.untracked.vdouble(0.00558245, 0.00520491, 0.00518877, 0.00514213, 0.00512142, 0.00512142, 0.00490732),
        mean_om    = cms.untracked.vdouble(1.67306 ,1.67284 ,1.6728 ,1.67268 ,1.67272 ,1.6727 ,1.67262 ,1.67269 ,1.67279),
        sigma_om   = cms.untracked.vdouble(0.00762585 ,0.00910153 ,0.00762863 ,0.00588798 ,0.00553094 ,0.00401328 ,0.00372463 ,0.00354882 ,0.00479006),
        ptcut_ks       = cms.untracked.vdouble(0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,9.0,12.0),
        ptcut_la       = cms.untracked.vdouble(0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,9.0,12.0),
        ptcut_xi         = cms.untracked.vdouble(1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10),
        ptcut_om         = cms.untracked.vdouble(1.0, 1.5, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10.0),
        useCentrality = cms.untracked.bool(False),
        cent_bin_low = cms.untracked.int32(60),
        cent_bin_high = cms.untracked.int32(100),
        rejectDaughter = cms.untracked.bool(True),
        doRap = cms.untracked.bool(False),
        doDauEff = cms.untracked.bool(False),
        doGenRef = cms.untracked.bool(False),
        doKs = cms.untracked.bool(False),
        doLa = cms.untracked.bool(False),
        doXi = cms.untracked.bool(True),
        doOm = cms.untracked.bool(False)
        )
