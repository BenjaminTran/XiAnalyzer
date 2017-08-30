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
process.load("XiAnalyzer.XiAnalyzer.massptproducer_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5000))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    # MinBias File 5 TeV
    #'/store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1000_1_Bgt.root',
    #'/store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1001_1_NY3.root',
    #'/store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1002_1_Qga.root',
    #'/store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1003_1_oEd.root',
    #'/store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1004_1_d8Z.root',
    #'/store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1005_1_ssw.root'
    # HM file 5 TeV
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Rereco_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1000_1_YQB.root',
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Rereco_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1001_1_Ab9.root',
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Rereco_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1002_1_mXl.root'
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Rereco_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1003_1_Kq0.root',
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Rereco_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1004_1_M1B.root',
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Rereco_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1005_1_c98.root',
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Rereco_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1006_1_gVF.root',
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Rereco_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1007_1_rlc.root',
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Rereco_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1008_1_my1.root',
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Rereco_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1009_1_jim.root',
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Rereco_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_100_1_s12.root ',
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Reverse_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1000_1_YyR.root',
    #'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1000_1_BPd.root'
    # HM file 8 TeV
    #'/store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_Pbp_V0Cascade_FullSkim_v1/170301_205443/0000/pPb_HM_1.root'
    #'/store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_Pbp_V0Cascade_FullSkim_v1/170301_205443/0000/pPb_HM_100.root'
    #'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_Pbp_V0Cascade_FullSkim_v1/170301_205443/0000/pPb_HM_90.root'
    #'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_pPb_V0Cascade_FullSkim_v1/170301_205416/0000/pPb_HM_100.root'
    'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_pPb_V0Cascade_FullSkim_v4/170803_222621/0000/pPb_HM_101.root'
   )
)

# Additional output definition
process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string(
                                    #'V0Correlation.root'
                                    'XiCorrelationRapidityPbpLoose.root'
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

# process.RapidityAnalysis = cms.Sequence(process.selectV0CandidatesNewlambdaRapidity*process.selectV0CandidatesNewkshortRapidity*process.selectV0CandidatesLowXiRapidity*process.v0CorrelationRapidity*process.xiCorrelationRapidity)

process.OmCorrAnalysis = cms.Sequence(process.selectOmegaCandidatesNew)

process.p = cms.Path(process.XiCorrAnalysisRapidityLoose)

process.schedule = cms.Schedule(process.p)


