import FWCore.ParameterSet.Config as cms
process = cms.Process("XiAna")

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )

### standard includes
process.load("Configuration.StandardSequences.Digi_cff")
#process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

### To use the TransientTrackBuilder
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_PostLS1_cff")

### GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.GlobalTag.globaltag = 'GR_P_V43F::All'
#process.GlobalTag.globaltag = 'GR_R_53_V21A::All'
process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v15'

#process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
#process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("XiAnalyzer.XiAnalyzer.xiselector_cff")
process.load("XiAnalyzer.XiAnalyzer.v0selector_cff")
process.load("XiAnalyzer.XiAnalyzer.xicorrelation_cff")
process.load("XiAnalyzer.XiAnalyzer.v0correlation_cff")
process.load("XiAnalyzer.XiAnalyzer.massptproducer_cff")
process.load("XiAnalyzer.XiAnalyzer.v0xiomttreeproducer_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #pPb
        # 'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_pPb_V0Cascade_FullSkim_v4/170803_222621/0000/pPb_HM_109.root'
        #PbPb
        'root://cmsxrootd.fnal.gov//store/user/davidlw/HIMinimumBias5/RecoSkim2015_pprereco_V0Cascade_Golden_v2/170302_083114/0000/PbPb_MB_109.root'
        ),
    secondaryFileNames = cms.untracked.vstring(
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2015/HIMinimumBias5/AOD/02May2016-v1/10000/52651877-FC24-E611-8885-003048F316DC.root'
        )
)

# Additional output definition
process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string(
                                        'V0CasTreePbPb.root'
                                        )
                                  )
# CORRELATION
# MinBias
#process.XiAnalysis = cms.Sequence(process.selectV0CandidatesLowXi*process.xiCorrelation)

#Tree producer
process.TreeProd = cms.Sequence(process.OmTreeProducerRapidityPbPb*process.XiTreeProducerRapidityPbPb*process.KsTreeProducerRapidityPbPb*process.LaTreeProducerRapidityPbPb)

process.LaOmTreeProd = cms.Sequence(process.LaTreeProducer + process.OmTreeProducer)

process.p = cms.Path(process.TreeProd)

process.schedule = cms.Schedule(process.p)


