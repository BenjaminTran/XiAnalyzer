import FWCore.ParameterSet.Config as cms

from XiAnalyzer.XiAnalyzer.highmultfilter_cfi import *

HighMultFilterPeriSub = HighMultFilter.clone(
        multMax = cms.double(35),
        multMin = cms.double(0)
)
