import xianalyzerheader as xi

try:
    with open( 'XiVarStore.dat', 'r' ) as fle:
        counter = int( fle.readline() )
except FileNotFoundError:
    counter = 0

#To submit edit counter, CollID, inputDataset

PDcounter = 1
collID = 'ppb'
#collID = 'pbp'

xi.config.General.requestName = 'HLT185_250FlowCascadev2' + collID + '2016PD' + str(PDcounter) + 'CorrelationJL' + str(counter)
xi.config.General.workArea = 'crab_dir/HLT185_250FlowCascadev2' + collID +'2016PD' + str(PDcounter) + 'Rap'

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



    #xi.config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #xi.config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #xi.config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #xi.config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #xi.config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #xi.config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #
    #xi.config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #xi.config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #xi.config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #xi.config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #xi.config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #xi.config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #MC
    xi.config.Data.inputDataset = '/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/davidlw-RecoSkim2016_pPb_V0_v1-2fc6918bc3c19ca88eae36cad5440243/USER'
    #xi.config.Data.inputDataset = '/ReggeGribovPartonMC_EposLHC_PbP_4080_4080_DataBS/davidlw-RecoSkim2016_Pbp_V0_v1-2fc6918bc3c19ca88eae36cad5440243/USER'

    with open( 'XiVarStore.dat', 'w' ) as fle:
        counter = counter + 1
        fle.write( str(counter) )
    submit(xi.config)
