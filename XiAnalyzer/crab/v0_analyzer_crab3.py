import v0analyzerheader as v0

#To submit edit PDcounter, CollID, inputDataset
try:
    with open( 'V0VarStore.dat', 'r' ) as fle:
        counter = int( fle.readline() )
except FileNotFoundError:
    counter = 0

PDcounter = 2
#collID = 'ppb'
collID = 'Pbp'

v0.config.General.workArea = 'crab_dir/HLT185_250FlowV0v2' + collID + '2016PD' + str(PDcounter) + 'Rap'
v0.config.General.requestName = 'HLT185_250Flow2016Correlation' + collID + 'V0PD' + str(PDcounter) + 'JL' + str(counter)

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


    #v0.config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #v0.config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #v0.config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #v0.config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #v0.config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #v0.config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #
    #v0.config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #v0.config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #v0.config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #v0.config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #v0.config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #v0.config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #v0.config.Data.userInputFiles = list(open('HMMC90.txt'))
    # MC
    v0.config.Data.inputDataset = '/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/davidlw-RecoSkim2016_pPb_V0_v1-2fc6918bc3c19ca88eae36cad5440243/USER'
    #v0.config.Data.inputDataset = '/ReggeGribovPartonMC_EposLHC_PbP_4080_4080_DataBS/davidlw-RecoSkim2016_Pbp_V0_v1-2fc6918bc3c19ca88eae36cad5440243/USER'

    with open( 'V0VarStore.dat', 'w' ) as fle:
        counter = counter + 1
        fle.write( str(counter) )
    submit(v0.config)
