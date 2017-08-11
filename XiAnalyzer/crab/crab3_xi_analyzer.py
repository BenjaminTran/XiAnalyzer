import os
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
# config.General.requestName = 'HLT185_250FlowCascadev2ppb2016PD6CorrelationJL6'
# config.General.workArea = 'HLT185_250FlowCascadev2ppb2016PD6Rap'
config.General.requestName = 'HLT185_250FlowCascadev2pbp2016PD3CorrelationJL9'
config.General.workArea = 'HLT185_250FlowCascadev2pbp2016PD3Rap'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = os.path.expandvars('$CMSSW_BASE/src/XiAnalyzer/XiAnalyzer/test/xianalysis_cfg.py')

config.section_("Data")
#config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#
#config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.userInputFiles = list(open('HMMC90.txt'))
config.Data.inputDBS = 'phys03'
#config.Data.primaryDataset = 'MinBias_TuneZ2star_7TeV_pythia6'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.totalUnits = 300
config.Data.unitsPerJob = 40000
config.Data.outLFNDirBase = '/store/group/phys_heavyions/btran/'
config.Data.useParent = True
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'Pbp2016_Cascade_Rereco_HM185_250Flow'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist = ['T2_US_MIT']

