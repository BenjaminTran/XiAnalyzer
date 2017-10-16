import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.v0cascorrelation_cfi import *

v0CasCorrelationRapidity = v0CasCorrelation.clone(
        multHigh       = cms.untracked.double(250),
        etaMin_trg = cms.untracked.double(-999999),
        etaMax_trg = cms.untracked.double(999999),
        ksCollection   = cms.untracked.InputTag('selectV0CandidatesNewkshortRapidity:Kshort'),
        laCollection   = cms.untracked.InputTag('selectV0CandidatesNewlambdaRapidity:Lambda'),
        xiCollection = cms.untracked.InputTag('selectV0CandidatesLowXiRapidity:Xi'),
        omCollection = cms.untracked.InputTag('selectOmegaCandidatesNewRapidity:Omega'),
        ptbin_n        = cms.untracked.int32(18),
        ptbin_n_cas        = cms.untracked.int32(13),
        sigma_ks       = cms.untracked.vdouble(0.00618431 ,0.0069063 ,0.00670385 ,0.00629998 ,0.0060034 ,0.0059126 ,0.00635576 ,0.00647386 ,0.00607519 ,0.00678218 ,0.00665193 ,0.00728848 ,0.00751979 ,0.0075107 ,0.00795764 ,0.00913178 ,0.00661404 ,0.00643239),
        mean_ks        = cms.untracked.vdouble(0.498204 ,0.497685 ,0.497711 ,0.497739 ,0.497693 ,0.497636 ,0.497599 ,0.497574 ,0.497576 ,0.497598 ,0.49755 ,0.497652 ,0.497674 ,0.49772 ,0.497758 ,0.497828 ,0.497605 ,0.498467),
        sigma_la       = cms.untracked.vdouble(0.00228409 ,0.00191783 ,0.00181896 ,0.00188397 ,0.00194954 ,0.00200742 ,0.00207498 ,0.0021983 ,0.00212885 ,0.00236541 ,0.00269782 ,0.00289265 ,0.0057154 ,0.00261388 ,0.00728782),
        mean_la        = cms.untracked.vdouble(1.11668 ,1.11618 ,1.11597 ,1.11593 ,1.1159 ,1.11584 ,1.11582 ,1.11582 ,1.11583 ,1.11585 ,1.11586 ,1.11593 ,1.11587 ,1.11616 ,1.11543),
        mean_xi    = cms.untracked.vdouble(1.32339 ,1.3227 ,1.32243 ,1.32234 ,1.32224 ,1.32218 ,1.32214 ,1.3221 ,1.32208),
        sigma_xi   = cms.untracked.vdouble(0.00546234 ,0.00477726 ,0.00445783 ,0.00400124 ,0.00382872 ,0.00379973 ,0.00392485 ,0.00394356 ,0.00484037),
        mean_om    = cms.untracked.vdouble(1.67306 ,1.67284 ,1.6728 ,1.67268 ,1.67272 ,1.6727 ,1.67262 ,1.67269 ,1.67279),
        sigma_om   = cms.untracked.vdouble(0.00762585 ,0.00910153 ,0.00762863 ,0.00588798 ,0.00553094 ,0.00401328 ,0.00372463 ,0.00354882 ,0.00479006),
        ptcut_ks       = cms.untracked.vdouble(0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0,15.0,20.0,25.0,30.0),
        ptcut_la       = cms.untracked.vdouble(0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0,15.0,20.0,25.0,30.0),
        ptcut_xi         = cms.untracked.vdouble(1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10),
        ptcut_om         = cms.untracked.vdouble(1.0, 1.5, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10.0),
        doRap = cms.untracked.bool(True)
        )

v0CasCorrelationRapidityPbPb = v0CasCorrelationRapidity.clone(
        useCentrality = cms.bool(True),
        ksCollection   = cms.untracked.InputTag('selectV0CandidatesNewkshortRapidityPbPb:Kshort'),
        laCollection   = cms.untracked.InputTag('selectV0CandidatesNewlambdaRapidityPbPb:Lambda'),
        xiCollection = cms.untracked.InputTag('selectV0CandidatesLowXiRapidityPbPb:Xi'),
        omCollection = cms.untracked.InputTag('selectOmegaCandidatesNewRapidityPbPb:Omega'),
        ptbin_n_cas = cms.untracked.int32(9),
        sigma_ks       = cms.untracked.vdouble(0.0066,0.0068,0.0063,0.0056,0.0049,0.0048,0.0047,0.0049,0.0053,0.0056,0.0057,0.0062,0.0068,0.0077),
        mean_ks        = cms.untracked.vdouble(0.49816 ,0.49754 ,0.49763 ,0.4977 ,0.49768 ,0.49764 ,0.49762 ,0.49757 ,0.49759 ,0.49761 ,0.49763 ,0.49763 ,0.4977 ,0.49769),
        sigma_la       = cms.untracked.vdouble(0.0024 ,0.0025 ,0.002 ,0.00188397 ,0.00194954 ,0.00200742 ,0.00207498 ,0.0021983 ,0.00212885 ,0.00236541),
        mean_la        = cms.untracked.vdouble(1.11668 ,1.1161 ,1.1159 ,1.11593 ,1.1159 ,1.11584 ,1.11582 ,1.11582 ,1.11583 ,1.11585),
        mean_xi    = cms.untracked.vdouble(1.3232 ,1.3226 ,1.3224 ,1.3223 ,1.3222 ,1.3222 ,1.3219, 1.3221),
        sigma_xi   = cms.untracked.vdouble(0.0043 ,0.0044,0.004 ,0.0035 ,0.0034 ,0.0034 ,0.0033, 0.0035),
        mean_om    = cms.untracked.vdouble(1.6728 ,1.6726 ,1.6727 ,1.6727 ,1.6728 ,1.6727),
        sigma_om   = cms.untracked.vdouble(0.011,0.0083 ,0.0061 ,0.0043 ,0.0049 ,0.003),
        )

v0CasCorrelationlooseRapidity = v0CasCorrelationRapidity.clone(
        laCollection   = cms.untracked.InputTag('selectV0CandidatesNewlambdalooseRapidity:Lambda'),
        ksCollection   = cms.untracked.InputTag('selectV0CandidatesNewkshortlooseRapidity:Kshort')
        )

v0CasCorrelationtightRapidity = v0CasCorrelationRapidity.clone(
        laCollection   = cms.untracked.InputTag('selectV0CandidatesNewlambdatightRapidity:Lambda'),
        ksCollection   = cms.untracked.InputTag('selectV0CandidatesNewkshorttightRapidity:Kshort')
        )

v0CasCorrelationRapidityMC = v0CasCorrelationRapidity.clone(
        multHigh = cms.untracked.double(999999),
        multLow = cms.untracked.double(-1)
        )

v0CasCorrelationRapidityMatchMC = v0CasCorrelationRapidity.clone(
        gnCollection   = cms.InputTag('genParticles'),
        multHigh = cms.untracked.double(999999),
        multLow = cms.untracked.double(-1),
        laCollection   = cms.untracked.InputTag('MatchCandidatesLambda:Lambda'),
        ksCollection   = cms.untracked.InputTag('MatchCandidatesKshort:Kshort'),
        doGenRef = cms.untracked.bool(True)
        #TO USE THIS MODULE ALSO COMMENT OUT ALL KINEMATIC CUTS IN SRC
        )
#hltHM = hltHM.clone()

v0CasCorrelationRapidityPeriSub = v0CasCorrelationRapidity.clone(
        multHigh = cms.untracked.double(20),
        multLow = cms.untracked.double(0)
        )
