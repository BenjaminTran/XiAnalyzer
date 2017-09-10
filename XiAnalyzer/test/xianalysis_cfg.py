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
process.load("XiAnalyzer.XiAnalyzer.omegaselector_cff")
process.load("XiAnalyzer.XiAnalyzer.xicorrelation_cff")
process.load("XiAnalyzer.XiAnalyzer.v0correlation_cff")
process.load("XiAnalyzer.XiAnalyzer.v0correlationmc_cff")
process.load("XiAnalyzer.XiAnalyzer.massptproducer_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5000))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_pPb_V0Cascade_FullSkim_v4/170803_222621/0000/pPb_HM_102.root'
    #'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_Pbp_V0Cascade_v1/170301_202152/0000/pPb_HM_105.root'),
    #secondaryFileNames = cms.untracked.vstring(
    #    '/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/3A1D69BC-01B7-E611-AD9A-FA163E98E135.root',
    #    '/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/3A8F14C0-F9B6-E611-9FDF-02163E0133B5.root',
    #    '/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/5878662E-F6B6-E611-A11A-FA163E02A339.root',
    #    '/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/AED413F1-F6B6-E611-AEBE-02163E0125FC.root',
    #    '/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/C408E634-F6B6-E611-AB11-FA163E592268.root',
    #    '/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/DECFA2B3-F7B6-E611-B733-FA163EB68B63.root'
    #)
    #'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_Pbp_V0Cascade_v1/170301_202152/0000/pPb_HM_101.root'
    #),
    #secondaryFileNames = cms.untracked.vstring(
    #    'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/001A0095-CFB6-E611-9975-02163E011DBF.root',
    #    'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/00F209E6-E2B6-E611-9416-FA163E59B88F.root',
    #    'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/4EA26194-CFB6-E611-9CB7-02163E014606.root',
    #    'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/6C66149B-CFB6-E611-A534-02163E0136E8.root',
    #    'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/72BD830F-D0B6-E611-A435-02163E0145A8.root',
    #    'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/8438FB8D-CFB6-E611-9DFF-FA163ED95C35.root',
    #    'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/B891E58E-CFB6-E611-A5A6-02163E0145A1.root',
    #    'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/CAF852B2-CFB6-E611-B077-02163E01356E.root',
    #    'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/CE1A628C-CFB6-E611-93DD-02163E0126F7.root',
    #    'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/CE41AA8B-CFB6-E611-B15F-FA163EA17EEB.root',
    #    'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/D47BEDDA-CFB6-E611-8E37-02163E0141C6.root'
    #)
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
   )
)

# Additional output definition
process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string(
                                    #'V0Correlation.root'
                                    #'XiCorrelationRapidityPbpLoose.root'
                                    #'XiCorrelationRapidityRECOClosure.root'
                                    'XiCorrelationRapidityCorrectionA.root'
				    )
                                  )
# CORRELATION
# MinBias
#process.XiAnalysis = cms.Sequence(process.selectV0CandidatesLowXi*process.xiCorrelation)

# HM
#process.XiAnalysis = cms.Sequence(process.hltHM*process.selectV0CandidatesLowXi*process.xiCorrelation)
process.XiCorrAnalysis = cms.Sequence(process.selectV0CandidatesLowXi*process.xiCorrelation)

process.XiCorrAnalysisRapidity = cms.Sequence(process.selectV0CandidatesLowXiRapidity*process.xiCorrelationRapidity)

process.XiCorrAnalysisRapidityLoose = cms.Sequence(process.selectV0CandidatesXiRapidityLoose*process.xiCorrelationRapidityLoose)

process.XiCorrAnalysisRapidityTight = cms.Sequence(process.selectV0CandidatesXiRapidityTight*process.xiCorrelationRapidityTight)

process.XiCorrAnalysisRapidityMC = cms.Sequence(process.selectV0CandidatesLowXiRapidity*process.xiCorrelationRapidityMC)

# process.RapidityAnalysis = cms.Sequence(process.selectV0CandidatesNewlambdaRapidity*process.selectV0CandidatesNewkshortRapidity*process.selectV0CandidatesLowXiRapidity*process.v0CorrelationRapidity*process.xiCorrelationRapidity)

process.OmCorrAnalysis = cms.Sequence(process.selectOmegaCandidatesNew)

process.GenCorrAnalysis = cms.Sequence(process.v0CorrelationMCRapidity)

process.p = cms.Path(process.XiCorrAnalysisRapidity)

process.schedule = cms.Schedule(process.p)


