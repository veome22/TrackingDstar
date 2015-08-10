from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DStar_MinBias_TuneMBR-July31'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'trkeffanalyzer_MC_cfg.py'

#config.Data.inputDataset = '/SingleNeutrino/RunIISpring15DR74-Asympt50nsRaw_MCRUN2_74_V9A-v2/GEN-SIM-RAW'
config.Data.inputDataset = '/MinBias_TuneMBR_13TeV-pythia8/RunIISpring15DR74-NoPURawReco_castor_MCRUN2_74_V8B-v1/GEN-SIM-RECO'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY_Run2015B.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'
#config.Data.lumiMask = 'Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'
#config.Data.lumiMask = 'Cert_251721_13TeV_lowPURun_JSON.txt'
#config.Data.runRange = '251000-252000'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.publishDataName = 'DStar_MinBias_TuneMBR-July31'

config.Site.storageSite = 'T2_CH_CERN' 
