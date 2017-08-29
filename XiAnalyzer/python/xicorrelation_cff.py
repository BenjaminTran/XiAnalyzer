import FWCore.ParameterSet.Config as cms

#from RiceHig.V0Analysis.xicorrelation_cfi import *
# For interactive
from XiAnalyzer.XiAnalyzer.xicorrelation_cfi import *

hltHM = hltHM.clone()


xiCorrelationRapidity = xiCorrelation.clone(
        xiCollection = cms.InputTag('selectV0CandidatesLowXiRapidity:Xi'),
        ptBin         = cms.vdouble(1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10, 20),
        PtBinNum      = cms.int32(10),
        xiMassMean    = cms.vdouble(1.32338 ,1.3227 ,1.32243 ,1.32235 ,1.32224 ,1.32218 ,1.32214 ,1.32211 ,1.32208 ,1.32236),
        xiMassSigma   = cms.vdouble(0.00512355 ,0.00460879 ,0.0042615 ,0.0038439 ,0.00360748 ,0.00352663 ,0.00364015 ,0.00369722 ,0.00377408 ,0.00425549)
        )

xiCorrelationRapidityLoose = xiCorrelationRapidity.clone(
        xiCollection = cms.InputTag('selectV0CandidatesXiRapidityLoose:Xi')
        )

xiCorrelationRapidityTight = xiCorrelationRapidity.clone(
        xiCollection = cms.InputTag('selectV0CandidatesXiRapidityTight:Xi')
        )
