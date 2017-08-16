import os
from WMCore.Configuration import Configuration
config = Configuration()

collID = 'pPb'
#collID = 'Pbp'
count = 22

config.section_("General")
#config.General.requestName = 'HLT185_250FlowCascadev2ppb2016PD1CorrelationJL22'
#config.General.requestName = 'HLT185_250FlowCascadev2pbp2016PD6CorrelationJL15'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = os.path.expandvars('$CMSSW_BASE/src/XiAnalyzer/XiAnalyzer/test/xianalysis_cfg.py')

config.section_("Data")
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.totalUnits = 300
config.Data.unitsPerJob = 50000
config.Data.outLFNDirBase = '/store/group/phys_heavyions/btran/'
config.Data.useParent = True
config.Data.publication = False
config.Data.publishDBS = 'phys03'
#config.Data.outputDatasetTag = 'Pbp2016_Cascade_Rereco_HM185_250Flow'
config.Data.outputDatasetTag = 'RapidityCut_Cascade_Rereco_HM185_250Flow'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_US_MIT','T2_US_Vanderbilt']
config.Site.whitelist = ['T2_US_Vanderbilt']

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)


    for num in range(1,6):
        counter = count + num
        if collID == 'pPb':
            DataSet = ['/PAHighMultiplicity1/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity2/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity3/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity4/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity5/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity6/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER']
            print 'Input Dataset is %r ' % (Dataset[num])
            config.Data.inputDataset = Dataset[num]
            config.General.workArea = 'HLT185_250FlowCascadev2ppb2016PD' + str(num+1) + 'Rap'
            config.General.requestName = 'HLT185_250FlowCascadev2ppb2016PD' + str(num+1) + 'CorrelationJL' + str(counter)
        else:
            DataSet = ['/PAHighMultiplicity1/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity2/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity3/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity4/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity5/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity6/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER']
            print 'Input Dataset is %r ' % (Dataset[num])
            config.Data.inputDataset = Dataset[num]
            config.General.workArea = 'HLT185_250FlowCascadev2pbp2016PD' + str(num+1) + 'Rap'
            config.General.requestName = 'HLT185_250FlowCascadev2pbp2016PD' + str(num+1) + 'CorrelationJL' + str(counter)
        submit(config)


