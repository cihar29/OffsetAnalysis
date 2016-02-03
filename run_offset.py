# PYTHON configuration file for class: OffsetAnalysis
# Author: C. Harrington
# Date:  19 - January - 2015

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = readFiles)
readFiles.extend( [

   '/store/mc/RunIIFall15DR76/SingleNeutrino/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/004252D9-C094-E511-928F-AC162DA8C2B0.root'

  # '/store/data/Run2015D/ZeroBias1/AOD/16Dec2015-v1/10000/180AD4E6-D0B0-E511-B655-0CC47A4D7662.root'  

] );

isMC = cms.bool(True)
OutputName = "_MC"

if isMC == False:

  process.load( "Configuration.Geometry.GeometryIdeal_cff" )
  process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
  process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
  from Configuration.AlCa.GlobalTag import GlobalTag
  process.GlobalTag = GlobalTag( process.GlobalTag, '76X_dataRun2_v15' )

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
  #process.primaryVertexFilter = cms.EDFilter( "GoodVertexFilter",
  #                                            vertexCollection = cms.InputTag('offlinePrimaryVertices'),
  #                                            minimumNDOF = cms.uint32(4),
  #                                            maxAbsZ = cms.double(24),
  #                                            maxd0 = cms.double(2) )

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
    numSkip = cms.int32(89),
    maxNPV = cms.int32(50),
    maxNPU = cms.int32(50),
    maxRho = cms.int32(50),
    isMC = isMC,
    reweight = cms.bool(False),
    pvTag = cms.InputTag("offlinePrimaryVertices"),
    puTag = cms.InputTag("addPileupInfo"),
    pfTag = cms.InputTag("particleFlow"),
    rhoTag = cms.InputTag("fixedGridRhoFastjetAll")
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
                        #process.primaryVertexFilter *
                        process.eeBadScFilter *
                        process.myseq )
