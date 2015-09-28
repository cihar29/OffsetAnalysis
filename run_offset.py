# PYTHON configuration file for class: OffsetAnalysis
# Author: C. Harrington
# Date:  19 - January - 2015

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = readFiles)
readFiles.extend( [

   '/store/mc/RunIISpring15DR74/SingleNeutrino/AODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/004C5D84-9708-E511-AB2E-00259073E384.root'

  # '/store/mc/RunIISpring15DR74/SingleNeutrino/AODSIM/Asympt50nsRaw_MCRUN2_74_V9A-v2/10000/00350BD2-FD08-E511-97F6-00074305CE2D.root'

  # '/store/data/Run2015B/ZeroBias/RECO/PromptReco-v1/000/251/251/00000/02A7FB5A-7427-E511-ABE1-02163E012284.root'

  # '/store/data/Run2015C/ZeroBias/RECO/PromptReco-v1/000/254/790/00000/0A17F7BE-064A-E511-AD50-02163E0141C4.root'

  # '/store/data/Run2015D/ZeroBias/RECO/PromptReco-v3/000/256/729/00000/1E7C4AFE-995F-E511-9D6B-02163E011ABA.root'

] );

isMC = cms.bool(True)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
OutputName = "_MC"

if isMC == False:

  process.load( "Configuration.Geometry.GeometryIdeal_cff" )
  process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
  process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
  from Configuration.AlCa.GlobalTag import GlobalTag
  # process.GlobalTag = GlobalTag( process.GlobalTag, '74X_dataRun2_Prompt_v1' ) # Run2015 B/C
  process.GlobalTag = GlobalTag( process.GlobalTag, '74X_dataRun2_Prompt_v2' ) # Run2015 D

  # ZeroBias Trigger
  process.HLTZeroBias =cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = cms.vstring('HLT_ZeroBias_*'),
    eventSetupPathsKey = cms.string(''),
    andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw = cms.bool(False)
  )

  #Beam Halo
  process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')

  #PV Filter
  process.primaryVertexFilter = cms.EDFilter( "GoodVertexFilter",
                                              vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                              minimumNDOF = cms.uint32(4),
                                              maxAbsZ = cms.double(24),
                                              maxd0 = cms.double(2) )

  #HCAL HBHE
  process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
  process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
  process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
    # inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'), # this is for the 50ns
    # inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResultRun2Loose'),  # this is for the 25ns
     inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResultRun2Tight'),  # this is for the 25ns
    reverseDecision = cms.bool(False)
  )

  #Bad EE Supercrystal filter
  process.load('RecoMET.METFilters.eeBadScFilter_cfi')

  OutputName = "_Data"

params = cms.PSet(
    PUstep = cms.double(1.),
    rhoStep = cms.double(1.),
    maxNPV = cms.int32(50),
    maxNPU = cms.int32(50),
    maxRho = cms.int32(50),
    isMC = isMC,
    reweight = cms.bool(True),
    pvTag = cms.InputTag("offlinePrimaryVertices"),
    puTag = cms.InputTag("addPileupInfo"),
    pfTag = cms.InputTag("particleFlow"),
    rhoTag = cms.InputTag("fixedGridRhoAll")
)

process.pfR4 = cms.EDAnalyzer("OffsetAnalysis",
    params,
    RootFileName = cms.string("Offset" + OutputName + "_R4.root"),
    coneDR = cms.double(0.4),
)

process.pfR5 = cms.EDAnalyzer("OffsetAnalysis",
    params,
    RootFileName = cms.string("Offset" + OutputName + "_R5.root"),
    coneDR = cms.double(0.5),
)

process.pfR7 = cms.EDAnalyzer("OffsetAnalysis",
    params,
    RootFileName = cms.string("Offset" + OutputName + "_R7.root"),
    coneDR = cms.double(0.7),
)

process.pfR8 = cms.EDAnalyzer("OffsetAnalysis",
    params,
    RootFileName = cms.string("Offset" + OutputName + "_R8.root"),
    coneDR = cms.double(0.8),
)

process.myseq = cms.Sequence( process.pfR4 )
                              #process.pfR5 *
                              #process.pfR7 *
                              #process.pfR8 )

if isMC :
  process.p = cms.Path( process.myseq )
else:
  process.p = cms.Path( process.HLTZeroBias * 
                        process.CSCTightHaloFilter *
                        process.HBHENoiseFilterResultProducer *
                        process.ApplyBaselineHBHENoiseFilter *
                        process.primaryVertexFilter *
                        process.eeBadScFilter *
                        process.myseq )
