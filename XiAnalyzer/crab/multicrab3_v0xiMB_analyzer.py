import v0xianalyzerheader as v0xi

#collID = 'pPb'
#collID = 'Pbp'
#collID = 'pPbMB'
#collID = 'PbpMB'

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

    for num in range(0,1):
        try:
            with open( 'V0XiVarStore.dat', 'r' ) as fle:
                counter = int( fle.readline() )
        except FileNotFoundError:
            counter = 0

        DataSet = ['/PAMinimumBias1/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias2/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias3/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias4/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias5/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias6/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias7/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias8/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias1/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias2/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias3/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias4/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias5/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias6/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias7/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias8/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias9/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias10/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias11/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias12/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias13/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias14/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias15/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias16/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias17/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias18/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias19/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER',
                '/PAMinimumBias20/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER']
        print 'Input Dataset is %r ' % (DataSet[num])
        v0xi.config.Data.inputDataset = DataSet[num]
        v0xi.config.General.workArea = 'crab_dir/HLT185_250FlowCombinedv2ppbMB2016PD' + str(num+1) + 'Rap'
        v0xi.config.General.requestName = 'HLT185_250Flow2016CorrelationpPbMBCombinedPD' + str(num+1) + 'JL' + str(counter)
        with open( 'V0XiVarStore.dat', 'w' ) as fle:
            counter = counter + 1
            fle.write( str(counter) )
        submit(v0xi.config)


