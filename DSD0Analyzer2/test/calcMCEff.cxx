// Author: Stephane Cooperstein (Princeton)


void calcMCEff() {

int nDStar_Kpi_den = 2503478; // total number of real D*'s produced
int nDStar_K3pi_den = 2690337; 

TChain *chain_Kpi = new TChain("analyzer/tree1");
TChain *chain_K3pi = new TChain("analyzer/tree2");

//chain_Kpi->Add(Form("%s/*.root", output_dir.c_str()));
//chain_K3pi->Add(Form("%s/*.root", output_dir.c_str()) );

chain_Kpi->Add("~/mount/cms/store/user/scoopers/d0Plusk2pi_GEN_SIM/crab_DStar_EnrichedMC_Kpi-Sep14/150914_172350/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/d0Plusk4pi_GEN_SIM_1/crab_DStar_EnrichedMC_K3pi-Sep16-v3/150918_144531/0000/*.root"); // none processed yet

//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July17-v1/150717_154349/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July17-v1/150717_154349/0000/TrkAnalysis_generalTracks_9.root");

//chain_Kpi->Print();

//std::string presel_Kpi = "cosAlphaKpi>0&&KpiTrkKpt>0.3&&KpiTrkpipt>0.3&&KpiTrkKchi2<2.5&&KpiTrkpichi2<2.5&&KpiTrkKnhits>=5&&KpiTrkpinhits>=5&&KpiTrkKdz<1&&KpiTrkpidz<1&&KpiTrkKdxy<0.1&&KpiTrkpidxy<0.1&&abs(D0MassKpi-1.864)<0.025&&D0VtxProb>0.01&&KpiTrkSpt>0.25&&KpiTrkSchi2<3&&KpiTrkSnhits>2&&KpiTrkSdz<1&&KpiTrkSdxy<0.1&&DSPtKpi>5.5";

//std::string presel_K3pi = "cosAlphaK3pi>0&&K3piTrkKpt>0.3&&K3piTrk1pipt>0.3&&K3piTrk2pipt>0.3&&K3piTrk3pipt>0.3&&K3piTrkKchi2<2.5&&K3piTrk1pichi2<2.5&&K3piTrk2pichi2<2.5&&K3piTrk3pichi2<2.5&&K3piTrkKnhits>=5&&K3piTrk1pinhits>=5&&K3piTrk2pinhits>=5&&K3piTrk3pinhits>=5&&K3piTrkKdz<1&&K3piTrk1pidz<1&&K3piTrk2pidz<1&&K3piTrk3pidz<1&&K3piTrkKdxy<0.1&&K3piTrk1pidxy<0.1&&K3piTrk2pidxy<0.1&&K3piTrk3pidxy<0.1&&abs(D0MassK3pi-1.864)<0.010&&D0VtxProb3>0.01&&K3piTrkSpt>0.25&&K3piTrkSchi2<3&&K3piTrkSnhits>2&&K3piTrkSdz<1&&K3piTrkSdxy<0.1&&DSPtK3pi>5.5";

std::cout<<"Doing time-consuming CopyTree() call..."<<std::endl;
//TTree *tree_Kpi = chain_Kpi->CopyTree(presel_Kpi.c_str());
TTree *tree_Kpi = (TTree*) chain_Kpi;
std::cout<<"Done copying Kpi tree"<<std::endl;
//TTree *tree_K3pi = chain_K3pi->CopyTree(presel_K3pi.c_str());
TTree *tree_K3pi = (TTree*) chain_K3pi;
std::cout<<"Done copying K3pi tree"<<std::endl;

std::vector<double> *MCDsDeltaR = 0;
std::vector<double> *MCDsDeltaR3 = 0;
std::vector<double> *D0MassKpi = 0;
std::vector<double> *DSMassKpi = 0;
std::vector<double> *cosAlphaKpi = 0;
std::vector<double> *KpiTrkKpt = 0;
std::vector<double> *KpiTrkpipt = 0;
std::vector<double> *KpiTrkKchi2 = 0;
std::vector<double> *KpiTrkpichi2 = 0;
std::vector<double> *KpiTrkKnhits = 0;
std::vector<double> *KpiTrkpinhits = 0;
std::vector<double> *KpiTrkKdz = 0;
std::vector<double> *KpiTrkpidz = 0;
std::vector<double> *KpiTrkKdxy = 0;
std::vector<double> *KpiTrkpidxy = 0;
std::vector<double> *D0VtxProb = 0;
std::vector<double> *KpiTrkSpt = 0;
std::vector<double> *KpiTrkSchi2 = 0;
std::vector<double> *KpiTrkSnhits = 0;
std::vector<double> *KpiTrkSdz = 0;
std::vector<double> *KpiTrkSdxy = 0;
std::vector<double> *DSPtKpi = 0;
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
int NKpiCand;
int NK3piCand;

tree_Kpi->SetBranchAddress("MCDsDeltaR",&MCDsDeltaR);
tree_K3pi->SetBranchAddress("MCDsDeltaR3",&MCDsDeltaR3);
tree_Kpi->SetBranchAddress("NKpiCand", &NKpiCand);
tree_K3pi->SetBranchAddress("NK3piCand", &NK3piCand);
tree_Kpi->SetBranchAddress("D0MassKpi", &D0MassKpi);
tree_Kpi->SetBranchAddress("DSMassKpi", &DSMassKpi);
tree_Kpi->SetBranchAddress("cosAlphaKpi", &cosAlphaKpi);
tree_Kpi->SetBranchAddress("KpiTrkKpt", &KpiTrkKpt);
tree_Kpi->SetBranchAddress("KpiTrkpipt", &KpiTrkpipt);
tree_Kpi->SetBranchAddress("KpiTrkKchi2", &KpiTrkKchi2);
tree_Kpi->SetBranchAddress("KpiTrkpichi2", &KpiTrkpichi2);
tree_Kpi->SetBranchAddress("KpiTrkKnhits", &KpiTrkKnhits);
tree_Kpi->SetBranchAddress("KpiTrkpinhits", &KpiTrkpinhits);
tree_Kpi->SetBranchAddress("KpiTrkKdz", &KpiTrkKdz);
tree_Kpi->SetBranchAddress("KpiTrkpidz", &KpiTrkpidz);
tree_Kpi->SetBranchAddress("KpiTrkKdxy", &KpiTrkKdxy);
tree_Kpi->SetBranchAddress("KpiTrkpidxy", &KpiTrkpidxy);
tree_Kpi->SetBranchAddress("D0VtxProb", &D0VtxProb);
tree_Kpi->SetBranchAddress("KpiTrkSpt", &KpiTrkSpt);
tree_Kpi->SetBranchAddress("KpiTrkSchi2", &KpiTrkSchi2);
tree_Kpi->SetBranchAddress("KpiTrkSnhits", &KpiTrkSnhits);
tree_Kpi->SetBranchAddress("KpiTrkSdz", &KpiTrkSdz);
tree_Kpi->SetBranchAddress("KpiTrkSdxy", &KpiTrkSdxy);
tree_Kpi->SetBranchAddress("DSPtKpi", &DSPtKpi);

tree_K3pi->SetBranchAddress("D0MassK3pi", &D0MassK3pi);
tree_K3pi->SetBranchAddress("DSMassK3pi", &DSMassK3pi);
tree_K3pi->SetBranchAddress("cosAlphaK3pi", &cosAlphaK3pi);
tree_K3pi->SetBranchAddress("K3piTrkKpt", &K3piTrkKpt);
tree_K3pi->SetBranchAddress("K3piTrk1pipt", &K3piTrk1pipt);
tree_K3pi->SetBranchAddress("K3piTrk2pipt", &K3piTrk2pipt);
tree_K3pi->SetBranchAddress("K3piTrk3pipt", &K3piTrk3pipt);
tree_K3pi->SetBranchAddress("K3piTrkKchi2", &K3piTrkKchi2);
tree_K3pi->SetBranchAddress("K3piTrk1pichi2", &K3piTrk1pichi2);
tree_K3pi->SetBranchAddress("K3piTrk2pichi2", &K3piTrk2pichi2);
tree_K3pi->SetBranchAddress("K3piTrk3pichi2", &K3piTrk3pichi2);
tree_K3pi->SetBranchAddress("K3piTrkKnhits", &K3piTrkKnhits);
tree_K3pi->SetBranchAddress("K3piTrk1pinhits", &K3piTrk1pinhits);
tree_K3pi->SetBranchAddress("K3piTrk2pinhits", &K3piTrk2pinhits);
tree_K3pi->SetBranchAddress("K3piTrk3pinhits", &K3piTrk3pinhits);
tree_K3pi->SetBranchAddress("K3piTrkKdz", &K3piTrkKdz);
tree_K3pi->SetBranchAddress("K3piTrk1pidz", &K3piTrk1pidz);
tree_K3pi->SetBranchAddress("K3piTrk2pidz", &K3piTrk2pidz);
tree_K3pi->SetBranchAddress("K3piTrk3pidz", &K3piTrk3pidz);
tree_K3pi->SetBranchAddress("K3piTrkKdxy", &K3piTrkKdxy);
tree_K3pi->SetBranchAddress("K3piTrk1pidxy", &K3piTrk1pidxy);
tree_K3pi->SetBranchAddress("K3piTrk2pidxy", &K3piTrk2pidxy);
tree_K3pi->SetBranchAddress("K3piTrk3pidxy", &K3piTrk3pidxy);
tree_K3pi->SetBranchAddress("D0VtxProb3", &D0VtxProb3);
tree_K3pi->SetBranchAddress("K3piTrkSpt", &K3piTrkSpt);
tree_K3pi->SetBranchAddress("K3piTrkSchi2", &K3piTrkSchi2);
tree_K3pi->SetBranchAddress("K3piTrkSnhits", &K3piTrkSnhits);
tree_K3pi->SetBranchAddress("K3piTrkSdz", &K3piTrkSdz);
tree_K3pi->SetBranchAddress("K3piTrkSdxy", &K3piTrkSdxy);
tree_K3pi->SetBranchAddress("DSPtK3pi", &DSPtK3pi);

int nDStar_Kpi = 0;
int nDStar_K3pi = 0;

double cosAlphaCut = 0.99;
double TrkPtCut = 0.5;
double Trkchi2Cut = 2.5;
double TrknhitsCut = 5.;
double TrkdzCut = 1.;
double TrkdxyCut = 0.1;
double massWindowCutKpi = 0.025;
double massWindowCutK3pi = 0.010;
double vtxProbCut = 0.03;
double TrkSptCut = 0.3;
double TrkSchi2Cut = 3.;
double TrkSnhitsCut = 2.;
double TrkSdzCut = 1.;
double TrkSdxyCut = 0.1;
double DSPtCut = 5.5;

std::cout<<"Kpi"<<std::endl;
int nentries = tree_Kpi->GetEntries();
std::cout<<nentries<<std::endl;
for (int i=0; i<nentries; i++) {
    if ((i % 1000) == 0) { std::cout<<"processing entry: "<<i<<std::endl; }
    tree_Kpi->GetEntry(i);
    for (int j=0; j<NKpiCand; j++) {
        if ( (*cosAlphaKpi)[j] < cosAlphaCut)  continue;
        if ( (*KpiTrkKpt)[j]<TrkPtCut || (*KpiTrkpipt)[j]<TrkPtCut) continue;
        if ( (*KpiTrkKchi2)[j]>Trkchi2Cut || (*KpiTrkpichi2)[j]>Trkchi2Cut) continue;
        if ( (*KpiTrkKnhits)[j]<=TrknhitsCut || (*KpiTrkpinhits)[j]<=TrknhitsCut) continue;
        if ( (*KpiTrkKdz)[j]>TrkdzCut || (*KpiTrkpidz)[j]>TrkdzCut ) continue;
        if ( (*KpiTrkKdxy)[j]>TrkdxyCut || (*KpiTrkpidxy)[j]>TrkdxyCut ) continue;
        if ( fabs((*D0MassKpi)[j] - 1.864) > massWindowCutKpi ) continue;
        if ( (*D0VtxProb)[j] < vtxProbCut) continue;
        if ( (*KpiTrkSpt)[j] < TrkSptCut) continue;
        if ( (*KpiTrkSchi2)[j] > TrkSchi2Cut) continue;
        if ( (*KpiTrkSnhits)[j] <= TrkSnhitsCut) continue;
        if ( (*KpiTrkSdz)[j] > TrkSdzCut) continue;
        if ( (*KpiTrkSdxy)[j] > TrkSdxyCut) continue;
        if ( (*DSPtKpi)[j] < DSPtCut) continue;
        if ((*MCDsDeltaR)[j] < 0.1) {
             // cand. is matched to gen. D*
             nDStar_Kpi++;
             break; // don't want to count more than one per event
        }
    } 
}

std::cout<<"K3pi"<<std::endl;
nentries = tree_K3pi->GetEntries();
std::cout<<nentries<<std::endl;
for (int i=0; i<nentries; i++) {
    if ((i % 1000) == 0) { std::cout<<"processing entry: "<<i<<std::endl; }
    tree_K3pi->GetEntry(i);
    for (int j=0; j<NK3piCand; j++) {
        if ( (*cosAlphaK3pi)[j] < cosAlphaCut)  continue;
        if ( (*K3piTrkKpt)[j]<TrkPtCut || (*K3piTrk1pipt)[j]<TrkPtCut || (*K3piTrk2pipt)[j]<TrkPtCut || (*K3piTrk3pipt)[j]<TrkPtCut) continue;
        if ( (*K3piTrkKchi2)[j]>Trkchi2Cut || (*K3piTrk1pichi2)[j]>Trkchi2Cut || (*K3piTrk2pichi2)[j]>Trkchi2Cut || (*K3piTrk3pichi2)[j]>Trkchi2Cut) continue;
        if ( (*K3piTrkKnhits)[j]<=TrknhitsCut || (*K3piTrk1pinhits)[j]<=TrknhitsCut || (*K3piTrk2pinhits)[j]<=TrknhitsCut || (*K3piTrk3pinhits)[j]<TrknhitsCut) continue;
        if ( (*K3piTrkKdz)[j]>TrkdzCut || (*K3piTrk1pidz)[j]>TrkdzCut || (*K3piTrk2pidz)[j]>TrkdzCut || (*K3piTrk3pidz)[j]>TrkdzCut) continue;
        if ( (*K3piTrkKdxy)[j]>TrkdxyCut || (*K3piTrk1pidxy)[j]>TrkdxyCut || (*K3piTrk2pidxy)[j]>TrkdxyCut || (*K3piTrk3pidxy)[j]>TrkdxyCut) continue;
        if ( fabs((*D0MassK3pi)[j] - 1.864) > massWindowCutK3pi ) continue;
        if ( (*D0VtxProb3)[j] < vtxProbCut) continue;
        if ( (*K3piTrkSpt)[j] < TrkSptCut) continue;
        if ( (*K3piTrkSchi2)[j] > TrkSchi2Cut) continue;
        if ( (*K3piTrkSnhits)[j] <= TrkSnhitsCut) continue;
        if ( (*K3piTrkSdz)[j] > TrkSdzCut) continue;
        if ( (*K3piTrkSdxy)[j] > TrkSdxyCut) continue;
        if ( (*DSPtK3pi)[j] < DSPtCut) continue;
        if ((*MCDsDeltaR3)[j] < 0.1) {
             // cand. is matched to gen. D*
             nDStar_K3pi++;
             break; // don't want to count more than one per event
        }
    } 
}

std::cout<<"Total number of two-body D* candidates produced: "<<nDStar_Kpi_den<<std::endl;
std::cout<<"Number of two-body D* candidates reconstructed: "<<nDStar_Kpi<<std::endl;
std::cout<<"Efficiency: "<<((1.0*nDStar_Kpi)/nDStar_Kpi_den)<<std::endl;

std::cout<<"Total number of four-body D* candidates produced: "<<nDStar_K3pi_den<<std::endl;
std::cout<<"Number of four-body D* candidates reconstructed: "<<nDStar_K3pi<<std::endl;
std::cout<<"Efficiency: "<<((1.0*nDStar_K3pi)/nDStar_K3pi_den)<<std::endl;
}
