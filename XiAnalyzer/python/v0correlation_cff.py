import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.v0correlation_cfi import *

v0CorrelationRapidity = v0Correlation.clone(
        multHigh       = cms.untracked.double(250),
        etaMin_trg = cms.untracked.double(-9999999),
        etaMax_trg = cms.untracked.double(9999999),
        laCollection   = cms.untracked.InputTag('selectV0CandidatesNewlambdaRapidity:Lambda'),
        ksCollection   = cms.untracked.InputTag('selectV0CandidatesNewkshortRapidity:Kshort'),
        ptbin_n        = cms.untracked.int32(18),
        sigma_ks       = cms.untracked.vdouble(0.00618431 ,0.0069063 ,0.00670385 ,0.00629998 ,0.0060034 ,0.0059126 ,0.00635576 ,0.00647386 ,0.00607519 ,0.00678218 ,0.00665193 ,0.00728848 ,0.00751979 ,0.0075107 ,0.00795764 ,0.00913178 ,0.00661404 ,0.00643239),
        mean_ks        = cms.untracked.vdouble(0.498204 ,0.497685 ,0.497711 ,0.497739 ,0.497693 ,0.497636 ,0.497599 ,0.497574 ,0.497576 ,0.497598 ,0.49755 ,0.497652 ,0.497674 ,0.49772 ,0.497758 ,0.497828 ,0.497605 ,0.498467),
        sigma_la       = cms.untracked.vdouble(0.00228409 ,0.00191783 ,0.00181896 ,0.00188397 ,0.00194954 ,0.00200742 ,0.00207498 ,0.0021983 ,0.00212885 ,0.00236541 ,0.00269782 ,0.00289265 ,0.0057154 ,0.00261388 ,0.00728782),
        mean_la        = cms.untracked.vdouble(1.11668 ,1.11618 ,1.11597 ,1.11593 ,1.1159 ,1.11584 ,1.11582 ,1.11582 ,1.11583 ,1.11585 ,1.11586 ,1.11593 ,1.11587 ,1.11616 ,1.11543),
        ptcut_ks       = cms.untracked.vdouble(0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0,15.0,20.0,25.0,30.0),
        ptcut_la       = cms.untracked.vdouble(0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0,15.0,20.0,25.0,30.0)
        )

v0CorrelationlooseRapidity = v0CorrelationRapidity.clone(
        laCollection   = cms.untracked.InputTag('selectV0CandidatesNewlambdalooseRapidity:Lambda'),
        ksCollection   = cms.untracked.InputTag('selectV0CandidatesNewkshortlooseRapidity:Kshort'),
        )

v0CorrelationtightRapidity = v0CorrelationRapidity.clone(
        laCollection   = cms.untracked.InputTag('selectV0CandidatesNewlambdatightRapidity:Lambda'),
        ksCollection   = cms.untracked.InputTag('selectV0CandidatesNewkshorttightRapidity:Kshort'),
        )

v0CorrelationRapidityMC = v0CorrelationRapidity.clone(
        multHigh = cms.untracked.double(99999999),
        multLow = cms.untracked.double(-1)
        )
#hltHM = hltHM.clone()
