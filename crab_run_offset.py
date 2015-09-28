from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'offset_analysis'
config.General.workArea = 'crab_projects/RunII_MCD_2509_R4'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_offset.py'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.inputFiles = ["pileup_25ns_2509.txt"]
config.JobType.outputFiles = ["Offset_MC_R4.root"]

config.section_("Data")
config.Data.inputDataset = '/SingleNeutrino/RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v2/AODSIM'
  #'/ZeroBias/Run2015D-PromptReco-v3/RECO'
  #'/ZeroBias/Run2015B-PromptReco-v1/RECO'
  #'/ZeroBias/Run2015C-PromptReco-v1/RECO'
  #'/SingleNeutrino/RunIISpring15DR74-Asympt50nsRaw_MCRUN2_74_V9A-v2/AODSIM'
  #'/SingleNeutrino/RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v2/AODSIM'
config.Data.splitting = 'FileBased' #'LumiBased'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-256869_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
config.Data.unitsPerJob = 10 #30 for data 10 for MC
config.Data.outLFNDirBase = '/store/user/charring/OffsetAnalysis/RunII_MCD_2509_R4'
config.Data.publication = False
#config.Data.ignoreLocality = True
#config.Data.publishDataName = 'offset_analysis'

config.section_("Site")
#config.Site.blacklist = ['T1_US_FNAL']
#config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = "T3_US_FNALLPC"

# source /cvmfs/cms.cern.ch/crab3/crab.csh
