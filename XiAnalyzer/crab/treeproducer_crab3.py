import treeproducerheader as tp

try:
    with open( 'TreeVarStore.dat', 'r' ) as fle:
        counter = int( fle.readline() )
except FileNotFoundError:
    counter = 0

#To submit edit counter, CollID, inputDataset

PDcounter = 2
#collID = 'ppb'
collID = 'pbp'

tp.config.General.requestName = 'HLT185_250FlowTree' + collID + '2016PD' + str(PDcounter) + 'JL' + str(counter)
tp.config.General.workArea = 'crab_dir/HLT185_250FlowTree' + collID +'2016PD' + str(PDcounter) + 'Rap'

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



    #tp.config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #tp.config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #tp.config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #tp.config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #tp.config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #tp.config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-RecoSkim2016_pPb_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #
    #tp.config.Data.inputDataset = '/PAHighMultiplicity1/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    tp.config.Data.inputDataset = '/PAHighMultiplicity2/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #tp.config.Data.inputDataset = '/PAHighMultiplicity3/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #tp.config.Data.inputDataset = '/PAHighMultiplicity4/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #tp.config.Data.inputDataset = '/PAHighMultiplicity5/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #tp.config.Data.inputDataset = '/PAHighMultiplicity6/davidlw-RecoSkim2016_Pbp_V0Cascade_v1-97be9aa52ea60cba5455e64649c12464/USER'
    #tp.config.Data.userInputFiles = list(open('HMMC90.txt'))

    with open( 'TreeVarStore.dat', 'w' ) as fle:
        counter = counter + 1
        fle.write( str(counter) )
    submit(tp.config)
