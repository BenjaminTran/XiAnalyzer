#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "XiAnalyzer/XiAnalyzer/interface/XiSelector.h"
#include "XiAnalyzer/XiAnalyzer/interface/V0Selector.h"
#include "XiAnalyzer/XiAnalyzer/interface/XiCorrelation.h"
#include "XiAnalyzer/XiAnalyzer/interface/V0Correlation.h"
#include "XiAnalyzer/XiAnalyzer/interface/XiMassPt.h"
#include "XiAnalyzer/XiAnalyzer/interface/XiOmTTree.h"

DEFINE_FWK_MODULE( XiSelector );
DEFINE_FWK_MODULE( V0Selector );
DEFINE_FWK_MODULE( XiCorrelation );
DEFINE_FWK_MODULE( V0Correlation );
DEFINE_FWK_MODULE( XiMassPt );
DEFINE_FWK_MODULE( XiOmTTree );
