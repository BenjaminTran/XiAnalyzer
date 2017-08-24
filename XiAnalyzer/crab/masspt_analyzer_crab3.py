import os
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'HLT185_250FlowMassPtppb2016MCJL28'
config.General.workArea = 'crab_dir/HLT185_250FlowMassPtv2ppb2016MCPD1'
#config.General.requestName = 'HLT185_250FlowMassPtpbp2016OmegaPD1JL28'
#config.General.workArea = 'crab_dir/HLT185_250FlowMassPtv2pbp2016OmegaPD2'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = os.path.expandvars('/afs/cern.ch/user/b/btran/work/CMSSW_8_0_24/src/XiAnalyzer/XiAnalyzer/test/masspt_cfg.py')

config.section_("Data")
config.Data.inputDataset = '/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/davidlw-RecoSkim2016_pPb_V0_v1-2fc6918bc3c19ca88eae36cad5440243/USER'
#config.Data.inputDataset = '/ReggeGribovPartonMC_EposLHC_PbP_4080_4080_DataBS/davidlw-RecoSkim2016_Pbp_V0_v1-2fc6918bc3c19ca88eae36cad5440243/USER'
#config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
#config.Data.userInputFiles = list(open('HMMC90.txt'))
config.Data.inputDBS = 'phys03'
#config.Data.primaryDataset = 'MinBias_TuneZ2star_7TeV_pythia6'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
#config.Data.splitting = 'EventAwareLumiBased'
#config.Data.unitsPerJob = 60000
#config.Data.totalUnits = 120000
config.Data.outLFNDirBase = '/store/group/phys_heavyions/btran/'
config.Data.useParent = True
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'pPb2016_MassPt_MC_HM185_250Flow'
#config.Data.outputDatasetTag = 'Pbp2016_MassPt_Omega_HM185_220Flow'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist = ['T2_US_MIT']

