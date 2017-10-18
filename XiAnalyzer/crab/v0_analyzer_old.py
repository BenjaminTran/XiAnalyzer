import os
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = os.path.expandvars('$CMSSW_BASE/src/XiAnalyzer/XiAnalyzer/test/v0casanalysis_cfg.py')
#config.General.requestName = 'HLT185_250FlowCascadev22016PDCorrelationJL2'
config.General.workArea = 'crab_dir/HLT185_250FlowCascadev22016PDRap'

config.section_("Data")
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
#config.Data.primaryDataset = 'MinBias_TuneZ2star_7TeV_pythia6'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.inputDataset = '/HIMinimumBias5/davidlw-RecoSkim2015_pprereco_V0Cascade_Golden_v2-a2a36526d6b050b4e6f00846a47a9f83/USER'
#config.Data.totalUnits = 300
config.Data.unitsPerJob = 90000
config.Data.outLFNDirBase = '/store/group/phys_heavyions/btran/'
config.Data.useParent = True
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'RapidityCut_Cascade_Rereco_HM185_250Flow'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_US_MIT']
config.Site.whitelist = ['T2_US_Vanderbilt']
