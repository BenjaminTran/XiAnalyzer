import FWCore.ParameterSet.Config as cms

#from RiceHig.V0Analysis.xicorrelation_cfi import *
# For interactive
from XiAnalyzer.XiAnalyzer.xicorrelation_cfi import *

hltHM = hltHM.clone()

omCorrelation = xiCorrelation.clone(
        multHigh      = cms.double(220),
        casCollection = cms.InputTag('selectOmegaCandidatesNew:Omega'),
        ptBin         = cms.vdouble(1.0, 1.5, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10.0),
        PtBinNum      = cms.int32(9),
        xiMassMean    = cms.vdouble(1.67316 ,1.6728 ,1.67278 ,1.67273 ,1.67283 ,1.67283 ,1.67269 ,1.67251 ,1.673),
        xiMassSigma   = cms.vdouble(0.00635809 ,0.0052194 ,0.00463206 ,0.0049998 ,0.00469137 ,0.00451139 ,0.00444609 ,0.00430358 ,0.00556587)
        )

omCorrelationRapidity = xiCorrelation.clone(
        casCollection = cms.InputTag('selectOmegaCandidatesNewRapidity:Omega'),
        ptBin         = cms.vdouble(1.0, 1.5, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10.0),
        PtBinNum      = cms.int32(9),
        xiMassMean    = cms.vdouble(1.67306 ,1.67284 ,1.6728 ,1.67268 ,1.67272 ,1.6727 ,1.67262 ,1.67269 ,1.67279),
        xiMassSigma   = cms.vdouble(0.00762585 ,0.00910153 ,0.00762863 ,0.00588798 ,0.00553094 ,0.00401328 ,0.00372463 ,0.00354882 ,0.00479006),
        doRap         = cms.bool(True)
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
