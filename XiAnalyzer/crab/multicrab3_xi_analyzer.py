import xianalyzerheader as xi

collID = 'pPb'
#collID = 'Pbp'

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
        try:
            with open( 'XiVarStore.dat', 'r' ) as fle:
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
            xi.config.Data.inputDataset = DataSet[num]
            xi.config.General.workArea = 'crab_dir/HLT185_250FlowCascadev2ppb2016PD' + str(num+1) + 'Rap'
            xi.config.General.requestName = 'HLT185_250Flow2016CorrelationpPbCascadePD' + str(num+1) + 'JL' + str(counter)
        else:
            DataSet = ['/PAHighMultiplicity1/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity2/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity3/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity4/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity5/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                       '/PAHighMultiplicity6/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER']
            print 'Input Dataset is %r ' % (DataSet[num])
            xi.config.Data.inputDataset = DataSet[num]
            xi.config.General.workArea = 'crab_dir/HLT185_250FlowCascadev2pbp2016PD' + str(num+1) + 'Rap'
            xi.config.General.requestName = 'HLT185_250Flow2016CorrelationPbpCascadePD' + str(num+1) + 'JL' + str(counter)
        with open( 'XiVarStore.dat', 'w' ) as fle:
            counter = counter + 1
            fle.write( str(counter) )
        submit(xi.config)


