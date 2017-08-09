import FWCore.ParameterSet.Config as cms

#from RiceHig.V0Analysis.xicorrelation_cfi import *
# For interactive
from XiAnalyzer.XiAnalyzer.xicorrelation_cfi import *



xiCorrelationRapidity = xiCorrelation.clone(
        ptBin         = cms.vdouble(1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 8.5, 10, 20),
        PtBinNum      = cms.int32(11),
        xiMassMean    = cms.vdouble(1.32338 ,1.3227 ,1.32243 ,1.32235 ,1.32224 ,1.32218 ,1.32214 ,1.32211 ,1.32207 ,1.32208 ,1.32236),
        xiMassSigma   = cms.vdouble(0.00559704, 0.00538643, 0.00517754, 0.00468905, 0.00448518, 0.00413725, 0.00434247, 0.00382751, 0.00463059, 0.00316714, 0.00352118)
        )
hltHM = hltHM.clone()
