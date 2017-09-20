import os
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = os.path.expandvars('/afs/cern.ch/user/b/btran/work/CMSSW_8_0_24/src/XiAnalyzer/XiAnalyzer/test/v0analysis_cfg.py')

config.section_("Data")
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
#config.Data.primaryDataset = 'MinBias_TuneZ2star_7TeV_pythia6'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 2
config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = 2000000
config.Data.unitsPerJob = 200000
config.Data.outLFNDirBase = '/store/group/phys_heavyions/btran/'
config.Data.useParent = True
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'RapidityCut_V0_Rereco_HM185_250Flow'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist = ['T2_US_MIT']
#config.Site.whitelist = ['T2_US_Vanderbilt']
