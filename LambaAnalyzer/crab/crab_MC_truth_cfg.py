from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#config.General.requestName = 'LambdaAnalysis_Feb28_QCD_GENSIMRECO_wPXB1RemovedVertexAttempt_Full'
#config.General.requestName = 'lambdaanalysis_apr7_singleneutrinogun_1000e'
config.General.requestName = 'lambdaanalysis_relval11'

#config.General.requestName = 'LambdaAnalysis_Dec7_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_GENSIMRECO'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'trkeffanalyzer_MC_GeneralTracks_cfg.py'

#config.Data.inputDataset = '/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/RunIIAutumn18DRPremix-102X_upgrade2018_realistic_v15_ext1-v2/GEN-SIM-RECO'
#config.Data.inputDataset = '/SingleNeutrino/RunIIAutumn18DRPremix-forRECO_102X_upgrade2018_realistic_v15_ext1-v1/GEN-SIM-RECO'
config.Data.inputDataset = '/RelValMinBias_13/CMSSW_10_2_0-102X_upgrade2018_realistic_v9_gcc7-v1/GEN-SIM-RECO'
##original relval config.Data.inputDataset = '/RelValQCD_FlatPt_15_3000HS_13/CMSSW_10_2_3-PU25ns_102X_upgrade2018_realistic_v11-v1/GEN-SIM-RECO'
#10mil custom MC sample (needs to get on disk): config.Data.inputDataset = '/SingleNeutrinoGun/RunIIAutumn18DR-102X_upgrade2018_realistic_v15_ext1-v2/GEN-SIM-RECO'
#config.Data.userInputFiles = ['/store/relval/CMSSW_10_2_3/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PU25ns_102X_upgrade2018_realistic_v11-v1/20000/FAC4072D-E9C0-7D47-B58C-7C9C8580E39A.root']
config.Data.useParent = True
config.Data.inputDBS = 'global'
#config.Data.splitting = 'LumiBased'
#config.Data.unitsPerJob = 20
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY_Run2015B.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'
#config.Data.lumiMask = 'Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'
#config.Data.lumiMask = 'Cert_251721_13TeV_lowPURun_JSON.txt'
#config.Data.lumiMask = 'Cert_314472-319851_13TeV_PromptReco_Collisions18_JSON.txt'
#config.Data.lumiMask = 'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
#config.Data.runRange = '251000-252000'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'DStar_ZeroBias_Run251721-lowPU-July27'

config.Site.storageSite = 'T3_US_FNALLPC' 
