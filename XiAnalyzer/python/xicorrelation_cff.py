import FWCore.ParameterSet.Config as cms

#from RiceHig.V0Analysis.xicorrelation_cfi import *
# For interactive
from XiAnalyzer.XiAnalyzer.xicorrelation_cfi import *

hltHM = hltHM.clone()

omCorrelation = xiCorrelation.clone(
        casCollection = cms.InputTag('selectOmCandidates:Omega'),
        ptBin = cms.vdouble(1.0, 1.5, 1.9, 2.3, 2.7, 3.3, 4.1, 5.0, 6.0, 8.0),
        PtBinNum = cms.int32(9),
        xiMassMean = cms.vdouble(1.67316 ,1.67273 ,1.67282 ,1.67275 ,1.67275 ,1.67285 ,1.67273 ,1.67271 ,1.6726),
        xiMassSigma = cms.vdouble(0.00635809 ,0.00531942 ,0.00476514 ,0.00492592 ,0.00484927 ,0.00462567 ,0.00438588 ,0.00457689 ,0.00457444)
        )


xiCorrelationRapidity = xiCorrelation.clone(
        casCollection = cms.InputTag('selectV0CandidatesLowXiRapidity:Xi'),
        ptBin         = cms.vdouble(1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10, 20),
        PtBinNum      = cms.int32(10),
        xiMassMean    = cms.vdouble(1.32338 ,1.3227 ,1.32243 ,1.32235 ,1.32224 ,1.32218 ,1.32214 ,1.32211 ,1.32208 ,1.32236),
        xiMassSigma   = cms.vdouble(0.00512355 ,0.00460879 ,0.0042615 ,0.0038439 ,0.00360748 ,0.00352663 ,0.00364015 ,0.00369722 ,0.00377408 ,0.00425549),
        doRap = cms.bool(True)
        )

xiCorrelationRapidityLoose = xiCorrelationRapidity.clone(
        casCollection = cms.InputTag('selectV0CandidatesXiRapidityLoose:Xi'),
        doRap = cms.bool(True)
        )

xiCorrelationRapidityTight = xiCorrelationRapidity.clone(
        casCollection = cms.InputTag('selectV0CandidatesXiRapidityTight:Xi'),
        doRap = cms.bool(True)
        )

xiCorrelationRapidityMC = xiCorrelationRapidity.clone(
        multHigh = cms.double(999999),
        multLow = cms.double(-1),
        doRap = cms.bool(True)
        )
