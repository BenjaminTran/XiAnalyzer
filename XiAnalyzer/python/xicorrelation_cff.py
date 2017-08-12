import FWCore.ParameterSet.Config as cms

#from RiceHig.V0Analysis.xicorrelation_cfi import *
# For interactive
from XiAnalyzer.XiAnalyzer.xicorrelation_cfi import *



xiCorrelationRapidity = xiCorrelation.clone(
        ptBin         = cms.vdouble(1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 10, 20),
        PtBinNum      = cms.int32(9),
        xiMassMean    = cms.vdouble(1.32338 ,1.3227 ,1.32243 ,1.32235 ,1.32224 ,1.32218 ,1.32214 ,1.3221 ,1.32236),
        xiMassSigma   = cms.vdouble(0.00512355 ,0.00460879 ,0.0042615 ,0.0038439 ,0.00360748 ,0.00352663 ,0.00364015 ,0.00370284 ,0.00425549)
        )
hltHM = hltHM.clone()
