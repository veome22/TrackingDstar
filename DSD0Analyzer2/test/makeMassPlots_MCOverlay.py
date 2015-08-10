## make basic plots from crab output
# Author: Stephane Cooperstein (Princeton)

import ROOT
import numpy
import tdrstyle, CMS_lumi
from math import sqrt

tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_13TeV = "40.02 pb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"

presel = "cosAlphaKpi>0.99"

plotlist = {}
#plotlist["nPV"]         = ["nPV", 30,0,30,"Kpi", presel]
#plotlist["D0MassK3pi1"] = ["D0MassK3pi1",40,1.5,2.1,"K3pi", presel]
#plotlist["DSMassK3pi1"] = ["DSMassK3pi1",40,1.7,2.5,"K3pi", presel]
#plotlist["D0MassKpi"] = ["D0MassKpi",40,1.65,2.05,"Kpi", presel]
#plotlist["DSMassKpi"] = ["DSMassKpi",50,1.7,2.5,"Kpi", presel]
#plotlist["NK3piCand"] = ["NK3piCand",50,0,50,"K3pi", presel]
#plotlist["NKpiCand"] = ["NKpiCand",50,0,50,"Kpi", presel]
#plotlist["dM_K3pi"] = ["DSMassK3pi1 - D0MassK3pi1", 40, 0.14, 0.16, "K3pi", presel]
#plotlist["dM_Kpi"] = ["DSMassKpi - D0MassKpi", 40, 0.14, 0.16, "Kpi", presel]
#plotlist["KpiTrkKpt"] = ["KpiTrkKpt", 40, 3, 10, "Kpi", presel]
#plotlist["KpiTrkKeta"] = ["KpiTrkKeta", 40, -5, 5, "Kpi", presel]
#plotlist["KpiTrkKphi"] = ["KpiTrkKphi", 40, -3.2, 3.2, "Kpi", presel]
plotlist["KpiTrkKdz"] = ["KpiTrkKdz", 40, 0, 1, "Kpi", presel]
#plotlist["KpiTrkSpt"] = ["KpiTrkSpt", 40, 0, 2, "Kpi", presel]
#plotlist["KpiTrkSeta"] = ["KpiTrkSeta", 40, -5, 5, "Kpi", presel]
#plotlist["KpiTrkSphi"] = ["KpiTrkSphi", 40, -3.2, 3.2, "Kpi", presel]
plotlist["KpiTrkSdz"] = ["KpiTrkSdz", 40, 0, 1, "Kpi", presel]
#plotlist["KpiTrkSnhits"] = ["KpiTrkSnhits", 30, 0, 30, "Kpi", presel]
#plotlist["KpiTrkSchi2"] = ["KpiTrkSchi2", 50, 0, 4, "Kpi", presel]
#plotlist["KpiTrkSdxy"]  = ["KpiTrkSdxy", 40, -0.1, 0.1, "Kpi", presel]
#plotlist["KpiTrkpipt"] = ["KpiTrkpipt", 40, 3, 10, "Kpi", presel]
#plotlist["KpiTrkpieta"] = ["KpiTrkpieta", 40, -5, 5, "Kpi", presel]
#plotlist["KpiTrkpiphi"] = ["KpiTrkpiphi", 40, -3.2, 3.2, "Kpi", presel]
plotlist["KpiTrkpidz"] = ["KpiTrkpidz", 40, 0, 1, "Kpi", presel]
#plotlist["K3piTrkKpt"] = ["K3piTrkKpt", 40, 3, 10, "K3pi", presel]
#plotlist["K3piTrkKeta"] = ["K3piTrkKeta", 40, -5, 5, "K3pi", presel]
#plotlist["K3piTrkKphi"] = ["K3piTrkKphi", 40, -3.2, 3.2, "K3pi", presel]
#plotlist["K3piTrkKdz"] = ["K3piTrkKdz", 40, 0, 1, "K3pi", presel]
#plotlist["K3piTrkSpt"] = ["K3piTrkSpt", 40, 0, 2, "K3pi", presel]
#plotlist["K3piTrkSeta"] = ["K3piTrkSeta", 40, -5, 5, "K3pi", presel]
#plotlist["K3piTrkSphi"] = ["K3piTrkSphi", 40, -3.2, 3.2, "K3pi", presel]
#plotlist["K3piTrkSdz"] = ["K3piTrkSdz", 40, 0, 1, "K3pi", presel]
#plotlist["K3piTrk1pipt"] = ["K3piTrk1pipt", 40, 3, 10, "K3pi", presel]
#plotlist["K3piTrk1pieta"] = ["K3piTrk1pieta", 40, -5, 5, "K3pi", presel]
#plotlist["K3piTrk1piphi"] = ["K3piTrk1piphi", 40, -3.2, 3.2, "K3pi", presel]
#plotlist["K3piTrk1pidz"] = ["K3piTrk1pidz", 40, 0, 1, "K3pi", presel]
#plotlist["K3piTrk2pipt"] = ["K3piTrk2pipt", 40, 3, 10, "K3pi", presel]
#plotlist["K3piTrk2pieta"] = ["K3piTrk2pieta", 40, -5, 5, "K3pi", presel]
#plotlist["K3piTrk2piphi"] = ["K3piTrk2piphi", 40, -3.2, 3.2, "K3pi", presel]
#plotlist["K3piTrk2pidz"] = ["K3piTrk2pidz", 40, 0, 3, "K3pi", presel]
#plotlist["K3piTrk3pipt"] = ["K3piTrk3pipt", 40, 3, 10, "K3pi", presel]
#plotlist["K3piTrk3pieta"] = ["K3piTrk3pieta", 40, -5, 5, "K3pi", presel]
#plotlist["K3piTrk3piphi"] = ["K3piTrk3piphi", 40, -3.2, 3.2, "K3pi", presel]
#plotlist["K3piTrk3pidz"] = ["K3piTrk3pidz", 40, 0, 3, "K3pi", presel]
#plotlist["D0VtxProb"]  = ["D0VtxProb",40,0,1.,"Kpi", presel]
#plotlist["D0VtxProb3"]  = ["D0VtxProb3",40,0,1.,"K3pi", presel]
#plotlist["D0etaKpi"]     = ["D0etaKpi",40,-5,5,"Kpi", presel]
#plotlist["D0etaK3pi"]     = ["D0etaK3pi",40,-5,5,"K3pi", presel]
#plotlist["D0phiKpi"]    = ["D0phiKpi",40,-3.2,3.2,"Kpi", presel]
#plotlist["D0phiK3pi"]    = ["D0phiK3pi",40,-3.2,3.2,"K3pi", presel]
#plotlist["D0PtKpi"]    = ["D0PtKpi",40,5,12,"Kpi", presel]
#plotlist["D0PtK3pi"]    = ["D0PtK3pi",40,5,12,"K3pi", presel]
#plotlist["DSPtKpi"]    = ["DSPtKpi",40,5,12,"Kpi", presel]
#plotlist["DSPtK3pi"]    = ["DSPtK3pi",40,5,12,"K3pi", presel]
#plotlist["DSetaKpi"]     = ["DSetaKpi",40,-5,5,"Kpi", presel]
#plotlist["DSetaK3pi"]     = ["DSetaK3pi",40,-5,5,"K3pi", presel]
#plotlist["DSphiKpi"]    = ["DSphiKpi",40,-3.2,3.2,"Kpi", presel]
#plotlist["DSphiK3pi"]    = ["DSphiK3pi",40,-3.2,3.2,"K3pi", presel]

#plotlist["cosAlphaKpi"] = ["cosAlphaKpi", 40, -1, 1, "Kpi", presel]
#plotlist["flightLengthKpi"] = ["flightLengthKpi", 40, -0.1, 0.1, "Kpi", presel]
#plotlist["cosAlphaK3pi"] = ["cosAlphaK3pi", 40, -1, 1, "K3pi", presel]
#plotlist["flightLengthK3pi"] = ["flightLengthK3pi", 40, -0.1, 0.1, "K3pi", presel]

#plotlist["MCDsDeltaRKpi"] = ["MCDsDeltaR", 100, 0, 3., "Kpi", presel]
#plotlist["MCDsDeltaRK3pi"] = ["MCDsDeltaR", 100, 0, 3 ,"K3pi", presel]

plots = plotlist.keys()
plots.sort()

#output_dir = "mount/cms/store/user/scoopers/ZeroBias1/DStar_ZeroBiasData-July15-v1/150715_111957/0000/" 
#output_dir = "mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBiasData-July15-v4/150715_213727/0000/"
#output_dir = "./"
#output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias_Run251721-lowPU-July17-v3/150717_174537/0000/"
#output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias_Run251721-lowPU-July17-v1/150717_155848/0000/"
#output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July17-v1/150717_154349/0000/"
#output_dir = " ~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBiasData-July15-v5/150716_100614/0000/"
#output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias_Run251721-lowPU-July17-v2/150717_160219/0000/"
#output_dir = "~/mount/cms/store/user/scoopers/JetHT/DStar_JetHT-July17-v1/150717_142543/0000/"
#output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July20/150720_140439/0000/"
#output_dir = " ~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias_Run251721-lowPU-July20-v2/150720_144629/0000/"
#output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July22/150722_114017/0000/"
output_dir = "mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July30/150730_133011/0000/"

chain_Kpi = ROOT.TChain("analyzer/tree1")
chain_K3pi = ROOT.TChain("analyzer/tree2")
mc_chain_Kpi = ROOT.TChain("analyzer/tree1")
mc_chain_K3pi = ROOT.TChain("analyzer/tree2")

chain_Kpi.Add("%s/*.root" % output_dir)
chain_K3pi.Add("%s/*.root" % output_dir )

mc_chain_Kpi.Add("mount/cms/store/user/scoopers/MinBias_TuneZ2star_13TeV-pythia6/DStar_MinBias_TuneZ2star-July21-v2-wHLT/150721_153230/0000/TrkAnalysis_MC_*.root")
mc_chain_Kpi.Add("mount/cms/store/user/scoopers/MinBias_TuneEE5C_13TeV-herwigpp/DStar_MinBias_TuneEE5C-July21/150721_133708/0000/TrkAnalysis_MC_*.root")
mc_chain_Kpi.Add("mount/cms/store/user/scoopers/MinBias_TuneCUETP8M1_13TeV-pythia8/DStar_MinBias_TuneCUETP8M1-July21-wHLT/150721_161541/0000/TrkAnalysis_MC_*.root")
mc_chain_Kpi.Add("mount/cms/store/user/scoopers/MinBias_TuneMBR_13TeV-pythia8/DStar_MinBias_TuneMBR-July21/150721_133938/0000/TrkAnalysis_MC_*.root")

mc_chain_K3pi.Add("mount/cms/store/user/scoopers/MinBias_TuneZ2star_13TeV-pythia6/DStar_MinBias_TuneZ2star-July21-v2-wHLT/150721_153230/0000/TrkAnalysis_MC_*.root")
mc_chain_K3pi.Add("mount/cms/store/user/scoopers/MinBias_TuneEE5C_13TeV-herwigpp/DStar_MinBias_TuneEE5C-July21/150721_133708/0000/TrkAnalysis_MC_*.root")
mc_chain_K3pi.Add("mount/cms/store/user/scoopers/MinBias_TuneCUETP8M1_13TeV-pythia8/DStar_MinBias_TuneCUETP8M1-July21-wHLT/150721_161541/0000/TrkAnalysis_MC_*.root")
mc_chain_K3pi.Add("mount/cms/store/user/scoopers/MinBias_TuneMBR_13TeV-pythia8/DStar_MinBias_TuneMBR-July21/150721_133938/0000/TrkAnalysis_MC_*.root")
#chain_Kpi.Print()

print chain_Kpi.GetEntries()
print chain_K3pi.GetEntries() 
for plot in plots:
    tree = ROOT.TTree("tree","tree")
    mc_tree = ROOT.TTree("mc_tree","mc_tree")
    if (plotlist[plot][4] == "K3pi"):
        continue
        tree = chain_K3pi
        mc_tree = mc_chain_K3pi
    elif (plotlist[plot][4] == "Kpi"):
        #continue
        tree = chain_Kpi
        mc_tree = mc_chain_Kpi
    else:
        print "invalid category type %s, skipping variables %s" % (plotlist[plot][4], plot)

    canv = ROOT.TCanvas("canv","canv")
    hist = ROOT.TH1F("temp",plot, plotlist[plot][1], plotlist[plot][2], plotlist[plot][3])
    hist_mc = ROOT.TH1F("temp_mc","%s_mc" % plot, plotlist[plot][1], plotlist[plot][2], plotlist[plot][3])

    tree.Draw("(%s)>>temp" % plotlist[plot][0], plotlist[plot][5])
    mc_tree.Draw("(%s)>>temp_mc" % plotlist[plot][0], plotlist[plot][5])

    #hist.GetYaxis().SetTitleOffset(1)
    #hist.GetXaxis().SetTitle("Reconstructed D0 Mass (GeV)")
    #hist.GetYaxis().SetTitle("Events / 10 MeV")

    # scale MC yield to match total number of data events
    hist_mc.Scale(hist.Integral()/hist_mc.Integral())

    hist.Draw("ep")
    hist_mc.Draw("hist same")

    CMS_lumi.CMS_lumi(canv, 4, 11)
    canv.SaveAs("%s/%s.png" % (plotlist[plot][4],plot) )
    #canv.SaveAs("%s/%s.C" % (plotlist[plot][4], plot) )

