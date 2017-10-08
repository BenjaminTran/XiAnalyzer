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
process.load("XiAnalyzer.XiAnalyzer.genselector_cff")
process.load("XiAnalyzer.XiAnalyzer.matchselector_cff")
process.load("XiAnalyzer.XiAnalyzer.xiselector_cff")
process.load("XiAnalyzer.XiAnalyzer.v0selector_cff")
process.load("XiAnalyzer.XiAnalyzer.xicorrelation_cff")
process.load("XiAnalyzer.XiAnalyzer.v0correlation_cff")
process.load("XiAnalyzer.XiAnalyzer.v0correlationmc_cff")
process.load("XiAnalyzer.XiAnalyzer.hadroncorrelationgen_cff")
process.load("XiAnalyzer.XiAnalyzer.massptproducer_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
   #'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_pPb_V0Cascade_FullSkim_v4/170803_222621/0000/pPb_HM_109.root'
   #For MC
   #'root://cmsxrootd.fnal.gov//store/user/davidlw/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/RecoSkim2016_pPb_V0_v1/170817_174330/0000/pPb_HM_1.root'
   #),
   # secondaryFileNames = cms.untracked.vstring(
   #     'root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/AODSIM/MB_80X_mcRun2_pA_v4-v2/120000/06A08BA8-4F10-E711-B5C0-00266CF3E3C4.root',
   #     'root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/AODSIM/MB_80X_mcRun2_pA_v4-v2/120000/06A49EE1-1A09-E711-B8FB-842B2B5C2299.root',
   #     'root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/AODSIM/MB_80X_mcRun2_pA_v4-v2/120000/06FE3920-490B-E711-A5AE-0CC47A01CAEA.root',
   #     'root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/AODSIM/MB_80X_mcRun2_pA_v4-v2/120000/080C92C2-7D08-E711-87B7-3417EBE6470E.root',
   #     'root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/AODSIM/MB_80X_mcRun2_pA_v4-v2/120000/0A000C53-F008-E711-A26F-001E67504D15.root',
   #     'root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/AODSIM/MB_80X_mcRun2_pA_v4-v2/120000/0A1413D5-2009-E711-8527-0CC47A7EEE76.root',
   #     'root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/AODSIM/MB_80X_mcRun2_pA_v4-v2/120000/0A193AD4-430B-E711-8BD6-02163E011842.root',
   #     'root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/AODSIM/MB_80X_mcRun2_pA_v4-v2/120000/003FC5B1-2009-E711-8AD2-0025904C6414.root',
   #     'root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/AODSIM/MB_80X_mcRun2_pA_v4-v2/120000/02BAB77A-6F0E-E711-9C60-E41D2D08DF30.root',
   #     'root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/AODSIM/MB_80X_mcRun2_pA_v4-v2/120000/061E97EF-C40B-E711-8781-549F3525A200.root'
   #)
   #For Peripheral Subtraction
   'root://cmsxrootd.fnal.gov//store/user/davidlw/PAMinimumBias1/RecoSkim2016_pPb_V0Cascade_v1/170302_094853/0000/pPb_HM_1.root'
   ),
    secondaryFileNames = cms.untracked.vstring(
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/479/00000/4E546F7B-DBAE-E611-B49E-FA163E63A392.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/480/00000/126AB4C4-2CAF-E611-93B9-FA163E41A46B.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/480/00000/2CB5DFE9-2BAF-E611-A0A9-FA163EA4BCD2.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/480/00000/52EDED08-28AF-E611-82AC-02163E0124B5.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/480/00000/9CDBDF3A-30AF-E611-8365-02163E01252E.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/480/00000/C8B078A7-24AF-E611-ADF7-FA163EDC366E.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/480/00000/F42B4F03-32AF-E611-9CC6-02163E011CAE.root'
   )
)

# Additional output definition
process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string(
                                    #'V0CorrelationCorrectMultB.root'
                                    #'V0CorrelationRapidityClosureRecoCollection.root'
                                    #'MatchV0ClosureB.root'
                                    #'MatchV0ClosureGenRef.root'
                                    'V0CorrelationRapidityPeriSubTestRAM.root'
				                    )
                                  )
# CORRELATION
# MinBias
#process.XiAnalysis = cms.Sequence(process.selectV0CandidatesLowXi*process.xiCorrelation)

# HM
#process.XiAnalysis = cms.Sequence(process.hltHM*process.selectV0CandidatesLowXi*process.xiCorrelation)
process.V0CorrAnalysis = cms.Sequence(process.selectV0CandidatesNewlambda*process.selectV0CandidatesNewkshort*process.v0Correlation)

process.V0CorrAnalysisRapidity = cms.Sequence(process.selectV0CandidatesNewlambdaRapidity*process.selectV0CandidatesNewkshortRapidity*process.v0CorrelationRapidity)

process.V0CorrAnalysisRapidityPeriSub = cms.Sequence(process.selectV0CandidatesNewlambdaRapidity*process.selectV0CandidatesNewkshortRapidity*process.v0CorrelationRapidityPeriSub)

process.V0CorrAnalysisRapidityMC = cms.Sequence(process.selectV0CandidatesNewlambdaRapidity*process.selectV0CandidatesNewkshortRapidity*process.v0CorrelationRapidityMC)

process.V0CorrAnalysisRapidityLoose = cms.Sequence(process.selectV0CandidatesNewlambdalooseRapidity*process.selectV0CandidatesNewkshortlooseRapidity*process.v0CorrelationlooseRapidity)

process.V0CorrAnalysisRapidityTight = cms.Sequence(process.selectV0CandidatesNewlambdatightRapidity*process.selectV0CandidatesNewkshorttightRapidity*process.v0CorrelationtightRapidity)

process.V0CorrAnalysisRapidityMCGen = cms.Sequence(process.v0CorrelationMCRapidity)
process.HadCorrAnalysisRapidityMCGen = cms.Sequence(process.HadronCorrelation)

process.genSelector = cms.Sequence(process.selectGenCandidatesKshort*process.selectGenCandidatesLambda*process.MatchCandidatesKshort*process.MatchCandidatesLambda*process.v0CorrelationRapidityMatchMC)

process.p = cms.Path(process.V0CorrAnalysisRapidityPeriSub)

process.schedule = cms.Schedule(process.p)


