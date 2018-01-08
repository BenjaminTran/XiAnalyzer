import os
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = os.path.expandvars('$CMSSW_BASE/src/XiAnalyzer/XiAnalyzer/test/v0casanalysis_cfg.py')

config.section_("Data")
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
#config.Data.primaryDataset = 'MinBias_TuneZ2star_7TeV_pythia6'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.totalUnits = 3500000
#######MB
#config.Data.unitsPerJob = 375000 #V0s 0-35
#config.Data.unitsPerJob = 375000 #Omega 0-35 Seems like could have gone higher
#config.Data.unitsPerJob = 550000 #Xi 0-35
#config.Data.unitsPerJob = 150000 #V0s + Xi 0-35
#######HM
#config.Data.unitsPerJob = 80000 #Xi
config.Data.unitsPerJob = 80000 #Xi
config.Data.outLFNDirBase = '/store/group/phys_heavyions/btran/'
config.Data.useParent = True
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'RapidityCut_V0Cas_Rereco_HM185_250Flow'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist = ['T2_US_MIT','T2_US_Vanderbilt']
#config.Site.whitelist = ['T2_US_Vanderbilt']
