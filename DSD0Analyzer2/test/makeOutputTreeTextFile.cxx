// create output file, from CRAB output, that can be used by the fitting script
// Author: Stephane Cooperstein (Princeton)

void makeOutputTreeTextFile() {
//output_dir = "mount/cms/store/user/scoopers/ZeroBias1/DStar_ZeroBiasData-July15-v1/150715_111957/0000/" 
//output_dir = "mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBiasData-July15-v4/150715_213727/0000/"
//output_dir = "./"
//output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias_Run251721-lowPU-July17-v3/150717_174537/0000/"
//output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias_Run251721-lowPU-July17-v1/150717_155848/0000/"
//std::string output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July17-v1/150717_154349/0000/";
//std::string output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July20/150720_140439/0000/";
//std::string output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July22/150722_114017/0000/";
//std::string output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July27/150727_100208/0000/";
std::string output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July30/150730_133011/0000/";
//output_dir = " ~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBiasData-July15-v5/150716_100614/0000/"
//output_dir = "~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias_Run251721-lowPU-July17-v2/150717_160219/0000/"
//output_dir = "~/mount/cms/store/user/scoopers/JetHT/DStar_JetHT-July17-v1/150717_142543/0000/"

TChain *chain_Kpi = new TChain("analyzer/tree1");
//TChain *chain_K3pi = new TChain("analyzer/tree2");

//chain_Kpi->Add(Form("%s/*.root", output_dir.c_str()));
//chain_K3pi->Add(Form("%s/*.root", output_dir.c_str()) );

chain_Kpi->Add("~/mount/cms/store/user/scoopers/MinBias_TuneZ2star_13TeV-pythia6/DStar_MinBias_TuneZ2star-July21-v2-wHLT/150721_153230/0000/TrkAnalysis_MC_*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/MinBias_TuneEE5C_13TeV-herwigpp/DStar_MinBias_TuneEE5C-July21/150721_133708/0000/TrkAnalysis_MC_*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/MinBias_TuneCUETP8M1_13TeV-pythia8/DStar_MinBias_TuneCUETP8M1-July21-wHLT/150721_161541/0000/TrkAnalysis_MC_*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/MinBias_TuneMBR_13TeV-pythia8/DStar_MinBias_TuneMBR-July21/150721_133938/0000/TrkAnalysis_MC_*.root");

//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July17-v1/150717_154349/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July17-v1/150717_154349/0000/TrkAnalysis_generalTracks_9.root");
//chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/MinBias_TuneCUETP8M1_13TeV-pythia8/DStar_MinBias_TuneCUETP8M1-July21/150721_133237/0000/TrkAnalysis_MC_8.root");

//chain_Kpi->Print();

//TFile *ofile = new TFile(Form("%s/outputTree.dat",output_dir.c_str()),"RECREATE");
//TFile *ofile = new TFile("outputTree.dat","RECREATE");

FILE *ofile1 = fopen("D0Mass.dat","w");
FILE *ofile2 = fopen("dM.dat","w");

std::cout<<"copying tree"<<std::endl;
TTree *otree = chain_Kpi->CopyTree("cosAlphaKpi>0.99");
std::cout<<"done copying tree"<<std::endl;
//otree->Print();
//std::vector<double> *cosAlphaKpi = 0;
//std::vector<double> *D0PtKpi = 0;
//std::vector<double> *D0etaKpi = 0;
//std::vector<double> *D0phiKpi = 0;
std::vector<double> *D0MassKpi = 0;
std::vector<double> *DSMassKpi = 0;
//std::vector<double> *D0VtxPosx = 0;
//std::vector<double> *D0VtxPosy = 0;
//std::vector<double> *D0VtxPosz = 0;
int NKpiCand;
//double PVx;
//double PVy;
//double PVz;

//TBranch *b_D0PtKpi;
//TBranch *b_D0etaKpi;
//TBranch *b_D0phiKpi;
TBranch *b_D0MassKpi;
TBranch *b_DSMassKpi;
//TBranch *b_D0VtxPosx;
//TBranch *b_D0VtxPosy;
//TBranch *b_D0VtxPosz;
//TBranch *b_cosAlphaKpi;

//otree->SetBranchAddress("cosAlphaKpi",&cosAlphaKpi, &b_cosAlphaKpi);
//otree->SetBranchAddress("D0PtKpi", &D0PtKpi, &b_D0PtKpi);
//otree->SetBranchAddress("D0etaKpi", &D0etaKpi, &b_D0etaKpi);
//otree->SetBranchAddress("D0phiKpi", &D0phiKpi, &b_D0phiKpi);
otree->SetBranchAddress("D0MassKpi", &D0MassKpi, &b_D0MassKpi);
otree->SetBranchAddress("DSMassKpi", &DSMassKpi, &b_DSMassKpi);
//otree->SetBranchAddress("D0VtxPosx", &D0VtxPosx, &b_D0VtxPosx);
//otree->SetBranchAddress("D0VtxPosy", &D0VtxPosy, &b_D0VtxPosy);
//otree->SetBranchAddress("D0VtxPosz", &D0VtxPosz, &b_D0VtxPosz);

otree->SetBranchAddress("NKpiCand", &NKpiCand);
//otree->SetBranchAddress("PVx", &PVx);
//otree->SetBranchAddress("PVy", &PVy);
//otree->SetBranchAddress("PVz", &PVz);


// Add extra variables

int nentries = otree->GetEntries();
std::cout<<nentries<<std::endl;
for (int i=0; i<nentries; i++) {
    if ((i % 1000) == 0) { std::cout<<"processing entry: "<<i<<std::endl; }
    otree->GetEntry(i);
    for (int j=0; j<NKpiCand; j++) {
         //std::cout<<(*DSMassKpi)[j]<<std::endl;
        //std::cout<<"DM = "<<((*DSMassKpi)[j] - (*D0MassKpi)[j])<<std::endl;
        fprintf(ofile1,"%f\n",((*D0MassKpi)[j]));
        fprintf(ofile2,"%f\n",((*DSMassKpi)[j] - (*D0MassKpi)[j]));
        
    } 
}

//ofile->cd();
//otree->Write();
//ofile->Write();
fclose(ofile1);
fclose(ofile2);

}
