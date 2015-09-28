// create output file, from CRAB output, that can be used by the fitting script
// Author: Stephane Cooperstein (Princeton)

void makeOutputTreeTextFile_K3pi() {
//output_dir = "mount/cms/store/user/scoopers/ZeroBias1/DStar_ZeroBiasData-July15-v1/150715_111957/0000/" 
//output_dir = "mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBiasData-July15-v4/150715_213727/0000/"
//output_dir = "./"
//output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias_Run251721-lowPU-July17-v3/150717_174537/0000/"
//output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias_Run251721-lowPU-July17-v1/150717_155848/0000/"
//std::string output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July17-v1/150717_154349/0000/";
//std::string output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July20/150720_140439/0000/";
//std::string output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July22/150722_114017/0000/";
//std::string output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July27/150727_100208/0000/";
//std::string output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July30/150730_133011/0000/";
std::string output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias-Aug18-noPV0/150818_115928/0000/";
//output_dir = " ~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBiasData-July15-v5/150716_100614/0000/"
//output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias_Run251721-lowPU-July17-v2/150717_160219/0000/"
//output_dir = "~/mount/cms/store/user/scoopers/JetHT/DStar_JetHT-July17-v1/150717_142543/0000/"

//TChain *chain_Kpi = new TChain("analyzer/tree1");
TChain *chain_K3pi = new TChain("analyzer/tree2");

//chain_Kpi->Add(Form("%s/*.root", output_dir.c_str()));
//chain_K3pi->Add(Form("%s/*.root", output_dir.c_str()) );
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/d0Plusk4pi_GEN_SIM_1/crab_DStar_EnrichedMC_K3pi-Sep16-v3/150918_144531/0000/*.root");
chain_K3pi->Add("150818_115928/0000/*.root");

//chain_Kpi->Add("~/mount/cms/store/user/scoopers/MinBias_TuneZ2star_13TeV-pythia6/DStar_MinBias_TuneZ2star-July21-v2-wHLT/150721_153230/0000/TrkAnalysis_MC_*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/MinBias_TuneEE5C_13TeV-herwigpp/DStar_MinBias_TuneEE5C-July21/150721_133708/0000/TrkAnalysis_MC_*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/MinBias_TuneCUETP8M1_13TeV-pythia8/DStar_MinBias_TuneCUETP8M1-July21-wHLT/150721_161541/0000/TrkAnalysis_MC_*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/MinBias_TuneMBR_13TeV-pythia8/DStar_MinBias_TuneMBR-July21/150721_133938/0000/TrkAnalysis_MC_*.root");

//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July17-v1/150717_154349/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July17-v1/150717_154349/0000/TrkAnalysis_generalTracks_9.root");
//chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/MinBias_TuneCUETP8M1_13TeV-pythia8/DStar_MinBias_TuneCUETP8M1-July21/150721_133237/0000/TrkAnalysis_MC_8.root");

//chain_Kpi->Print();

//TFile *ofile = new TFile(Form("%s/outputTree.dat",output_dir.c_str()),"RECREATE");
//TFile *ofile = new TFile("outputTree.dat","RECREATE");

FILE *ofile1 = fopen("D0Mass_K3pi.dat","w");
FILE *ofile2 = fopen("dM_K3pi.dat","w");

std::cout<<"copying tree"<<std::endl;
//TTree *otree = chain_K3pi->CopyTree("cosAlphaK3pi>0.99&&K3piTrkKpt>0.3&&K3piTrk1pipt>0.3&&K3piTrk2pipt>0.3&&K3piTrk3pipt>0.3&&K3piTrkKchi2<2.5&&K3piTrk1pichi2<2.5&&K3piTrk2pichi2<2.5&&K3piTrk3pichi2<2.5&&K3piTrkKnhits>=5&&K3piTrk1pinhits>=5&&K3piTrk2pinhits>=5&&K3piTrk3pinhits>=5&&K3piTrkKdz<1&&K3piTrk1pidz<1&&K3piTrk2pidz<1&&K3piTrk3pidz<1&&K3piTrkKdxy<0.1&&K3piTrk1pidxy<0.1&&K3piTrk2pidxy<0.1&&K3piTrk3pidxy<0.1&&abs(D0MassK3pi-1.864)<0.010&&D0VtxProb3>0.01&&K3piTrkSpt>0.25&&K3piTrkSchi2<3&&K3piTrkSnhits>2&&K3piTrkSdz<1&&K3piTrkSdxy<0.1&&DSPtK3pi>5.5");
//TTree *otree = chain_K3pi->CopyTree("cosAlphaK3pi>0.99&&K3piTrkKpt>0.5&&K3piTrk1pipt>0.5&&K3piTrk2pipt>0.5&&K3piTrk3pipt>0.5&&K3piTrkKchi2<2.5&&K3piTrk1pichi2<2.5&&K3piTrk2pichi2<2.5&&K3piTrk3pichi2<2.5&&K3piTrkKnhits>=5&&K3piTrk1pinhits>=5&&K3piTrk2pinhits>=5&&K3piTrk3pinhits>=5&&K3piTrkKdz<1&&K3piTrk1pidz<1&&K3piTrk2pidz<1&&K3piTrk3pidz<1&&K3piTrkKdxy<0.1&&K3piTrk1pidxy<0.1&&K3piTrk2pidxy<0.1&&K3piTrk3pidxy<0.1&&abs(D0MassK3pi-1.864)<0.012&&D0VtxProb3>0.03&&K3piTrkSpt>0.3&&K3piTrkSchi2<3&&K3piTrkSnhits>=2&&K3piTrkSdz<1&&K3piTrkSdxy<0.1&&DSPtK3pi>5.5");
TTree *otree = chain_K3pi->CloneTree();
std::cout<<"done copying tree"<<std::endl;
//otree->Print();
//std::vector<double> *cosAlphaKpi = 0;
//std::vector<double> *D0PtKpi = 0;
//std::vector<double> *D0etaKpi = 0;
//std::vector<double> *D0phiKpi = 0;
std::vector<double> *D0MassK3pi = 0;
std::vector<double> *DSMassK3pi = 0;
std::vector<double> *cosAlphaK3pi = 0;
std::vector<double> *K3piTrkKpt = 0;
std::vector<double> *K3piTrk1pipt = 0;
std::vector<double> *K3piTrk2pipt = 0;
std::vector<double> *K3piTrk3pipt = 0;
std::vector<double> *K3piTrkKchi2 = 0;
std::vector<double> *K3piTrk1pichi2 = 0;
std::vector<double> *K3piTrk2pichi2 = 0;
std::vector<double> *K3piTrk3pichi2 = 0;
std::vector<double> *K3piTrkKnhits = 0;
std::vector<double> *K3piTrk1pinhits = 0;
std::vector<double> *K3piTrk2pinhits = 0;
std::vector<double> *K3piTrk3pinhits = 0;
std::vector<double> *K3piTrkKdz = 0;
std::vector<double> *K3piTrk1pidz = 0;
std::vector<double> *K3piTrk2pidz = 0;
std::vector<double> *K3piTrk3pidz = 0;
std::vector<double> *K3piTrkKdxy = 0;
std::vector<double> *K3piTrk1pidxy = 0;
std::vector<double> *K3piTrk2pidxy = 0;
std::vector<double> *K3piTrk3pidxy = 0;
std::vector<double> *D0VtxProb3 = 0;
std::vector<double> *K3piTrkSpt = 0;
std::vector<double> *K3piTrkSchi2 = 0;
std::vector<double> *K3piTrkSnhits = 0;
std::vector<double> *K3piTrkSdz = 0;
std::vector<double> *K3piTrkSdxy = 0;
std::vector<double> *DSPtK3pi = 0;
//std::vector<double> *D0VtxPosx = 0;
//std::vector<double> *D0VtxPosy = 0;
//std::vector<double> *D0VtxPosz = 0;
int NK3piCand;
//double PVx;
//double PVy;
//double PVz;

//otree->SetBranchAddress("cosAlphaKpi",&cosAlphaKpi, &b_cosAlphaKpi);
//otree->SetBranchAddress("D0PtKpi", &D0PtKpi, &b_D0PtKpi);
//otree->SetBranchAddress("D0etaKpi", &D0etaKpi, &b_D0etaKpi);
//otree->SetBranchAddress("D0phiKpi", &D0phiKpi, &b_D0phiKpi);
otree->SetBranchAddress("D0MassK3pi", &D0MassK3pi);
otree->SetBranchAddress("DSMassK3pi", &DSMassK3pi);
otree->SetBranchAddress("cosAlphaK3pi", &cosAlphaK3pi);
otree->SetBranchAddress("K3piTrkKpt", &K3piTrkKpt);
otree->SetBranchAddress("K3piTrk1pipt", &K3piTrk1pipt);
otree->SetBranchAddress("K3piTrk2pipt", &K3piTrk2pipt);
otree->SetBranchAddress("K3piTrk3pipt", &K3piTrk3pipt);
otree->SetBranchAddress("K3piTrkKchi2", &K3piTrkKchi2);
otree->SetBranchAddress("K3piTrk1pichi2", &K3piTrk1pichi2);
otree->SetBranchAddress("K3piTrk2pichi2", &K3piTrk2pichi2);
otree->SetBranchAddress("K3piTrk3pichi2", &K3piTrk3pichi2);
otree->SetBranchAddress("K3piTrkKnhits", &K3piTrkKnhits);
otree->SetBranchAddress("K3piTrk1pinhits", &K3piTrk1pinhits);
otree->SetBranchAddress("K3piTrk2pinhits", &K3piTrk2pinhits);
otree->SetBranchAddress("K3piTrk3pinhits", &K3piTrk3pinhits);
otree->SetBranchAddress("K3piTrkKdz", &K3piTrkKdz);
otree->SetBranchAddress("K3piTrk1pidz", &K3piTrk1pidz);
otree->SetBranchAddress("K3piTrk2pidz", &K3piTrk2pidz);
otree->SetBranchAddress("K3piTrk3pidz", &K3piTrk3pidz);
otree->SetBranchAddress("K3piTrkKdxy", &K3piTrkKdxy);
otree->SetBranchAddress("K3piTrk1pidxy", &K3piTrk1pidxy);
otree->SetBranchAddress("K3piTrk2pidxy", &K3piTrk2pidxy);
otree->SetBranchAddress("K3piTrk3pidxy", &K3piTrk3pidxy);
otree->SetBranchAddress("D0VtxProb3", &D0VtxProb3);
otree->SetBranchAddress("K3piTrkSpt", &K3piTrkSpt);
otree->SetBranchAddress("K3piTrkSchi2", &K3piTrkSchi2);
otree->SetBranchAddress("K3piTrkSnhits", &K3piTrkSnhits);
otree->SetBranchAddress("K3piTrkSdz", &K3piTrkSdz);
otree->SetBranchAddress("K3piTrkSdxy", &K3piTrkSdxy);
otree->SetBranchAddress("DSPtK3pi", &DSPtK3pi);

//otree->SetBranchAddress("D0VtxPosx", &D0VtxPosx, &b_D0VtxPosx);
//otree->SetBranchAddress("D0VtxPosy", &D0VtxPosy, &b_D0VtxPosy);
//otree->SetBranchAddress("D0VtxPosz", &D0VtxPosz, &b_D0VtxPosz);

otree->SetBranchAddress("NK3piCand", &NK3piCand);
//otree->SetBranchAddress("PVx", &PVx);
//otree->SetBranchAddress("PVy", &PVy);
//otree->SetBranchAddress("PVz", &PVz);


double cosAlphaCut = 0.;
double TrkPtCut = 0.3;
double Trkchi2Cut = 2.5;
double TrknhitsCut = 5.;
double TrkdzCut = 1.;
double TrkdxyCut = 0.1;
double massWindowCut = 0.010;
double vtxProbCut = 0.01;
double TrkSptCut = 0.25;
double TrkSchi2Cut = 3.;
double TrkSnhitsCut = 2.;
double TrkSdzCut = 1.;
double TrkSdxyCut = 0.1;
double DSPtCut = 5.5;

int nentries = otree->GetEntries();
std::cout<<nentries<<std::endl;
for (int i=0; i<nentries; i++) {
    if ((i % 1000) == 0) { std::cout<<"processing entry: "<<i<<std::endl; }
    otree->GetEntry(i);
    for (int j=0; j<NK3piCand; j++) {
        if ( (*cosAlphaK3pi)[j] < cosAlphaCut)  continue;
        if ( (*K3piTrkKpt)[j]<TrkPtCut || (*K3piTrk1pipt)[j]<TrkPtCut || (*K3piTrk2pipt)[j]<TrkPtCut || (*K3piTrk3pipt)[j]<TrkPtCut) continue;
        if ( (*K3piTrkKchi2)[j]>Trkchi2Cut || (*K3piTrk1pichi2)[j]>Trkchi2Cut || (*K3piTrk2pichi2)[j]>Trkchi2Cut || (*K3piTrk3pichi2)[j]>Trkchi2Cut) continue;
        if ( (*K3piTrkKnhits)[j]<=TrknhitsCut || (*K3piTrk1pinhits)[j]<=TrknhitsCut || (*K3piTrk2pinhits)[j]<=TrknhitsCut || (*K3piTrk3pinhits)[j]<TrknhitsCut) continue;
        if ( (*K3piTrkKdz)[j]>TrkdzCut || (*K3piTrk1pidz)[j]>TrkdzCut || (*K3piTrk2pidz)[j]>TrkdzCut || (*K3piTrk3pidz)[j]>TrkdzCut) continue;
        if ( (*K3piTrkKdxy)[j]>TrkdxyCut || (*K3piTrk1pidxy)[j]>TrkdxyCut || (*K3piTrk2pidxy)[j]>TrkdxyCut || (*K3piTrk3pidxy)[j]>TrkdxyCut) continue;
        if ( fabs((*D0MassK3pi)[j] - 1.864) > massWindowCut ) continue;
        if ( (*D0VtxProb3)[j] < vtxProbCut) continue;
        if ( (*K3piTrkSpt)[j] < TrkSptCut) continue; 
        if ( (*K3piTrkSchi2)[j] > TrkSchi2Cut) continue;
        if ( (*K3piTrkSnhits)[j] <= TrkSnhitsCut) continue;
        if ( (*K3piTrkSdz)[j] > TrkSdzCut) continue;
        if ( (*K3piTrkSdxy)[j] > TrkSdxyCut) continue;
        if ( (*DSPtK3pi)[j] < DSPtCut) continue;
         //std::cout<<(*DSMassKpi)[j]<<std::endl;
        //std::cout<<"DM = "<<((*DSMassKpi)[j] - (*D0MassKpi)[j])<<std::endl;
        fprintf(ofile1,"%f\n",((*D0MassK3pi)[j]));
        fprintf(ofile2,"%f\n",((*DSMassK3pi)[j] - (*D0MassK3pi)[j]));
        
    } 
}

//ofile->cd();
//otree->Write();
//ofile->Write();
fclose(ofile1);
fclose(ofile2);

}
