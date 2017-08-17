import os
from WMCore.Configuration import Configuration
config = Configuration()

collID = 'pPb'
#collID = 'Pbp'

config.section_("General")

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = os.path.expandvars('/afs/cern.ch/user/b/btran/work/CMSSW_8_0_24/src/XiAnalyzer/XiAnalyzer/test/v0analysis_cfg.py')

config.section_("Data")
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 60000
#config.Data.totalUnits = 300
config.Data.outLFNDirBase = '/store/group/phys_heavyions/btran/'
config.Data.useParent = True
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'RapidityCut_V0_Rereco_HM185_250Flow'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist = ['T2_US_MIT','T2_US_Vanderbilt']

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


    for num in range(0,6):
        try:
            with open( 'V0VarStore.dat', 'r' ) as fle:
                counter = int( fle.readline() )
        except FileNotFoundError:
            counter = 0

        if collID == 'pPb':
            DataSet = ['/PAHighMultiplicity1/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity2/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity3/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity4/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity5/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity6/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER']
            print 'Input Dataset is %r ' % (DataSet[num])
            config.Data.inputDataset = DataSet[num]
            config.General.workArea = 'crab_dir/HLT185_250FlowV0v2ppb2016PD' + str(num+1) + 'Rap'
            config.General.requestName = 'HLT185_250FlowV0v2ppb2016PD' + str(num+1) + 'CorrelationJL' + str(counter)
        else:
            DataSet = ['/PAHighMultiplicity1/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity2/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity3/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity4/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity5/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity6/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER']
            print 'Input Dataset is %r ' % (DataSet[num])
            config.Data.inputDataset = DataSet[num]
            config.General.workArea = 'crab_dir/HLT185_250FlowV0v2pbp2016PD' + str(num+1) + 'Rap'
            config.General.requestName = 'HLT185_250FlowV0v2pbp2016PD' + str(num+1) + 'CorrelationJL' + str(counter)
        with open( 'V0VarStore.dat', 'w' ) as fle:
            counter = counter + 1
            fle.write( str(counter) )
        submit(config)


