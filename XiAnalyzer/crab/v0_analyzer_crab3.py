import os
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'HLT185_220FlowKslav2ppb2016PD5CorrelationJL8'
config.General.workArea = 'HLT185_220FlowKslav2ppb2016PD5'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = os.path.expandvars('/afs/cern.ch/user/b/btran/work/CMSSW_8_0_24/src/XiAnalyzer/XiAnalyzer/test/v0analysis_cfg.py')

config.section_("Data")
#config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.userInputFiles = list(open('HMMC90.txt'))
config.Data.inputDBS = 'phys03'
#config.Data.primaryDataset = 'MinBias_TuneZ2star_7TeV_pythia6'
config.Data.splitting = 'FileBased'
#config.Data.totalUnits = 300
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/group/phys_heavyions/btran/'
config.Data.useParent = True
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'pPb2016_Cascade_Rereco_HM185_220Flow'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist = ['T2_US_MIT']

