#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "XiAnalyzer/XiAnalyzer/interface/XiSelector.h"
#include "XiAnalyzer/XiAnalyzer/interface/V0Selector.h"
#include "XiAnalyzer/XiAnalyzer/interface/OmegaSelector.h"
#include "XiAnalyzer/XiAnalyzer/interface/XiCorrelation.h"
#include "XiAnalyzer/XiAnalyzer/interface/V0Correlation.h"
#include "XiAnalyzer/XiAnalyzer/interface/MassPtProducer.h"
#include "XiAnalyzer/XiAnalyzer/interface/V0XiOmTTreeProducer.h"

DEFINE_FWK_MODULE(XiSelector);
DEFINE_FWK_MODULE(V0Selector);
DEFINE_FWK_MODULE(OmegaSelector);
DEFINE_FWK_MODULE(XiCorrelation);
DEFINE_FWK_MODULE(V0Correlation);
DEFINE_FWK_MODULE(MassPtProducer);
DEFINE_FWK_MODULE(V0XiOmTTreeProducer);
