import FWCore.ParameterSet.Config as cms

process = cms.Process("XiAna")

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )

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
process.load("XiAnalyzer.XiAnalyzer.highmultfilter_cfi")
process.load("XiAnalyzer.XiAnalyzer.v0selector_cff")
process.load("XiAnalyzer.XiAnalyzer.omegaselector_cff")
process.load("XiAnalyzer.XiAnalyzer.v0cascorrelation_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #pPb MB
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
   #PbPb
   #'root://cmsxrootd.fnal.gov//store/user/davidlw/HIMinimumBias5/RecoSkim2015_pprereco_V0Cascade_Golden_v2/170302_083114/0000/PbPb_MB_107.root'
   #),
   # secondaryFileNames = cms.untracked.vstring(
   #     'root://cmsxrootd.fnal.gov//store/hidata/HIRun2015/HIMinimumBias5/AOD/02May2016-v1/10000/5A15CD8E-BE24-E611-BE6B-842B2B6AECDD.root'
   #)
)

# Additional output definition
process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string(
                                    #'V0CasCorrelationPbPbTest.root'
                                    #'V0CasCorrelationMBTEST.root'
                                    #'V0CasCorrelationPbPb.root'
                                    #'V0CasCorrelationPeriSubOmegaTest.root'
                                    #'MBomCorrelation_0_35_Rebin_v1.root'
                                    #'MBV0Correlation_0_35_V0DifferenceTest.root'
                                    'MBV0Correlation_0_35_Pt8p5_10.root'
                                    #'MBXiCorrelation_0_35_wMultFilter.root'
				    )
                                  )
# CORRELATION

# HM
#process.XiAnalysis = cms.Sequence(process.hltHM*process.selectV0CandidatesLowXi*process.xiCorrelation)
process.RapidityAnalysis = cms.Sequence(process.selectV0CandidatesNewlambdaRapidity*process.selectV0CandidatesNewkshortRapidity*process.selectV0CandidatesLowXiRapidity*process.selectOmegaCandidatesNewRapidity*process.v0CasCorrelationRapidity)

# MB
process.RapidityAnalysisPeriSub = cms.Sequence(process.selectV0CandidatesNewlambdaRapidity*process.selectV0CandidatesNewkshortRapidity*process.selectV0CandidatesLowXiRapidity*process.selectOmegaCandidatesNewRapidity*process.v0CasCorrelationRapidityPeriSub)

process.RapidityAnalysisPeriSubV0 = cms.Sequence(process.HighMultFilter*process.selectV0CandidatesNewlambdaRapidity*process.selectV0CandidatesNewkshortRapidity*process.v0CasCorrelationRapidityPeriSub)

process.RapidityAnalysisPeriSubV0Xi = cms.Sequence(process.selectV0CandidatesNewlambdaRapidity*process.selectV0CandidatesNewkshortRapidity*process.selectV0CandidatesLowXiRapidity*process.v0CasCorrelationRapidityPeriSub)

process.RapidityAnalysisPeriSubXi = cms.Sequence(process.HighMultFilter*process.selectV0CandidatesLowXiRapidity*process.v0CasCorrelationRapidityPeriSub)

process.RapidityAnalysisPeriSubXiOmega = cms.Sequence(process.selectV0CandidatesLowXiRapidity*process.selectOmegaCandidatesNewRapidity*process.v0CasCorrelationRapidityPeriSub)

process.RapidityAnalysisPeriSubOmega = cms.Sequence(process.selectOmegaCandidatesNewRapidity*process.v0CasCorrelationRapidityPeriSub)

# PbPb
process.RapidityAnalysisPbPb = cms.Sequence(process.selectV0CandidatesNewlambdaRapidityPbPb*process.selectV0CandidatesNewkshortRapidityPbPb*process.selectV0CandidatesLowXiRapidityPbPb*process.selectV0CandidatesLowOmegaRapidityPbPb*process.v0CasCorrelationRapidityPbPb)

# Reco Cuts pPb
process.RapidityAnalysisLoose = cms.Sequence(process.HighMultFilter*process.selectV0CandidatesLowXiRapidityLoose*process.v0CasCorrelationlooseRapidity)

process.RapidityAnalysisTight = cms.Sequence(process.HighMultFilter*process.selectV0CandidatesLowXiRapidityTight*process.v0CasCorrelationtightRapidity)

#process.p = cms.Path(process.selectOmegaCandidatesNewRapidity*process.v0CasCorrelationRapidity)
process.p = cms.Path(process.RapidityAnalysisLoose)

process.schedule = cms.Schedule(process.p)


