// create output file, from CRAB output, that can be used by the fitting script
// Author: Stephane Cooperstein (Princeton)
#include<cmath>
#include "TChain.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TFile.h"
#include<iostream>

using namespace std;

void makeOutputTreeTextFile(char* name, bool doData=true) {
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

TChain *chain_Kpi = new TChain("analyzer/tree1");
TChain *chain_K3pi = new TChain("analyzer/tree2");

//chain_Kpi->Add(Form("%s/*.root", output_dir.c_str()));
//chain_K3pi->Add(Form("%s/*.root", output_dir.c_str()) );
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/d0Plusk4pi_GEN_SIM_1/crab_DStar_EnrichedMC_K3pi-Sep16-v3/150918_144531/0000/*.root");
//chain_Kpi->Add("150818_115928/0000/*.root");

//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015D_Oct8/151008_150034/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015D_Oct8/151008_150209/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015D_Oct8/151008_150254/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015D_Oct8/151008_150335/0000/*.root");
/*chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015D_Oct13/151013_204514/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015D_Oct13/151013_204901/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015D_Oct13/151013_204942/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015D_Oct13/151013_205021/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015C_Oct8/151008_151934/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015C_Oct8/151008_152023/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015C_Oct8/151008_152103/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015C_Oct8/151008_152141/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias5/crab_DStar_ZeroBias5_Run2015C_Oct8/151008_152222/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias6/crab_DStar_ZeroBias6_Run2015C_Oct8/151008_152304/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias7/crab_DStar_ZeroBias7_Run2015C_Oct8/151008_152526/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias8/crab_DStar_ZeroBias8_Run2015C_Oct8/151008_152636/0000/*.root");


//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015D_Oct8/151008_150034/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015D_Oct8/151008_150209/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015D_Oct8/151008_150254/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015D_Oct8/151008_150335/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015D_Oct13/151013_204514/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015D_Oct13/151013_204901/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015D_Oct13/151013_204942/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015D_Oct13/151013_205021/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015C_Oct8/151008_151934/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015C_Oct8/151008_152023/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015C_Oct8/151008_152103/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015C_Oct8/151008_152141/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias5/crab_DStar_ZeroBias5_Run2015C_Oct8/151008_152222/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias6/crab_DStar_ZeroBias6_Run2015C_Oct8/151008_152304/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias7/crab_DStar_ZeroBias7_Run2015C_Oct8/151008_152526/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias8/crab_DStar_ZeroBias8_Run2015C_Oct8/151008_152636/0000/*.root");

//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015D_Oct5-FillByPV/151005_225806/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015D_Oct13/151013_204335/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015C_Oct5-FillByPV/151005_225536/0000/*.root");
chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015B_Oct5-FillByPV/151005_210907/0000/*.root");

//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015D_Oct5-FillByPV/151005_225806/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015D_Oct13/151013_204335/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015C_Oct5-FillByPV/151005_225536/0000/*.root");
chain_K3pi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015B_Oct5-FillByPV/151005_210907/0000/*.root");
*/

/*chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015B_Oct17-v2/151018_032218/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015C_Oct17/151018_031702/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015D_Oct17-v2/151018_025955/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015C_Oct17/151018_031805/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015D_Oct17-v2/151018_030348/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015C_Oct17/151018_031844/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015D_Oct17-v2/151018_030428/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015C_Oct17/151018_031922/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015D_Oct17-v2/151018_030557/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015C_Oct17/151018_032002/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015D_Oct17-v2/151018_030647/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias5/crab_DStar_ZeroBias5_Run2015C_Oct17-v2/151018_032454/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias6/crab_DStar_ZeroBias6_Run2015C_Oct17/151018_032534/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias7/crab_DStar_ZeroBias7_Run2015C_Oct17/151018_032613/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias8/crab_DStar_ZeroBias8_Run2015C_Oct17/151018_032655/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015B_Oct17-v2/151018_032218/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015C_Oct17/151018_031702/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015D_Oct17-v2/151018_025955/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015C_Oct17/151018_031805/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015D_Oct17-v2/151018_030348/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015C_Oct17/151018_031844/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015D_Oct17-v2/151018_030428/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015C_Oct17/151018_031922/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015D_Oct17-v2/151018_030557/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015C_Oct17/151018_032002/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015D_Oct17-v2/151018_030647/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias5/crab_DStar_ZeroBias5_Run2015C_Oct17-v2/151018_032454/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias6/crab_DStar_ZeroBias6_Run2015C_Oct17/151018_032534/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias7/crab_DStar_ZeroBias7_Run2015C_Oct17/151018_032613/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias8/crab_DStar_ZeroBias8_Run2015C_Oct17/151018_032655/0000/*.root");
*/

//chain_Kpi->Add("~/mount/cms/store/user/scoopers/d0Plusk2pi_GEN_SIM/crab_DStar_EnrichedMC_Kpi-Oct17/151018_215251/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/d0Plusk4pi_GEN_SIM_1/crab_DStar_EnrichedMC_K3pi-Oct17/151018_220323/0000/*.root");

/*chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015C_Nov12/151112_192953/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015C_Nov12/151112_193034/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015C_Nov12/151112_193123/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015C_Nov12/151112_193158/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015C_Nov12/151112_193239/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias5/crab_DStar_ZeroBias5_Run2015C_Nov12/151112_193324/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias6/crab_DStar_ZeroBias6_Run2015C_Nov12/151112_193356/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias7/crab_DStar_ZeroBias7_Run2015C_Nov12/151112_193430/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias8/crab_DStar_ZeroBias8_Run2015C_Nov12/151112_193508/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015C_Nov12/151112_192953/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015C_Nov12/151112_193034/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015C_Nov12/151112_193123/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015C_Nov12/151112_193158/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015C_Nov12/151112_193239/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias5/crab_DStar_ZeroBias5_Run2015C_Nov12/151112_193324/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias6/crab_DStar_ZeroBias6_Run2015C_Nov12/151112_193356/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias7/crab_DStar_ZeroBias7_Run2015C_Nov12/151112_193430/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias8/crab_DStar_ZeroBias8_Run2015C_Nov12/151112_193508/0000/*.root");
*/

//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015C_Nov16/151116_163029/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015C_Nov16/151116_163343/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015C_Nov16/151116_163723/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015C_Nov16/151116_164126/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015C_Nov16/151116_164417/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias5/crab_DStar_ZeroBias5_Run2015C_Nov16/151116_165857/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias6/crab_DStar_ZeroBias6_Run2015C_Nov16/151116_175815/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias7/crab_DStar_ZeroBias6_7un2015C_Nov16/151116_175859/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias8/crab_DStar_ZeroBias8_Run2015C_Nov16/151116_180100/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015C_Nov16/151116_163029/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015C_Nov16/151116_163343/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015C_Nov16/151116_163723/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015C_Nov16/151116_164126/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015C_Nov16/151116_164417/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias5/crab_DStar_ZeroBias5_Run2015C_Nov16/151116_165857/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias6/crab_DStar_ZeroBias6_Run2015C_Nov16/151116_175815/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias7/crab_DStar_ZeroBias6_7un2015C_Nov16/151116_175859/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias8/crab_DStar_ZeroBias8_Run2015C_Nov16/151116_180100/0000/*.root");

if (doData) {
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015C_Dec7/151207_222810/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015C_Dec7/151207_223906/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015C_Dec7/151207_223939/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015C_Dec7/151207_224055/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015C_Dec7/151208_041138/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias5/crab_DStar_ZeroBias5_Run2015C_Dec7/151208_041317/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias6/crab_DStar_ZeroBias6_Run2015C_Dec7/151208_041407/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias7/crab_DStar_ZeroBias7_Run2015C_Dec7/151208_041514/0000/*.root");
chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias8/crab_DStar_ZeroBias8_Run2015C_Dec7/151208_041601/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias/crab_DStar_ZeroBias_Run2015C_Dec7/151207_222810/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias1/crab_DStar_ZeroBias1_Run2015C_Dec7/151207_223906/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias2/crab_DStar_ZeroBias2_Run2015C_Dec7/151207_223939/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias3/crab_DStar_ZeroBias3_Run2015C_Dec7/151207_224055/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias4/crab_DStar_ZeroBias4_Run2015C_Dec7/151208_041138/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias5/crab_DStar_ZeroBias5_Run2015C_Dec7/151208_041317/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias6/crab_DStar_ZeroBias6_Run2015C_Dec7/151208_041407/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias7/crab_DStar_ZeroBias7_Run2015C_Dec7/151208_041514/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/ZeroBias8/crab_DStar_ZeroBias8_Run2015C_Dec7/151208_041601/0000/*.root");
}
else {
chain_Kpi->Add("~/mount/cms/store/user/scoopers/d0Plusk2pi_GEN_SIM/crab_DStar_EnrichedMC_Kpi-Decy7/151207_223230/0000/*.root");
chain_K3pi->Add("~/mount/cms/store/user/scoopers/d0Plusk4pi_GEN_SIM_1/crab_DStar_EnrichedMC_K3pi-Dec7/151207_223419/0000/*.root");
}
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/d0Plusk2pi_GEN_SIM/crab_DStar_EnrichedMC_Kpi-Oct14/151014_160055/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/d0Plusk4pi_GEN_SIM_1/crab_DStar_EnrichedMC_K3pi-Oct14/151014_155926/0000/*.root");

//chain_Kpi->Add("~/mount/cms/store/user/scoopers/d0Plusk2pi_GEN_SIM/crab_DStar_EnrichedMC_Kpi-Sep14/150914_172350/0000/*.root");
//chain_K3pi->Add("~/mount/cms/store/user/scoopers/d0Plusk4pi_GEN_SIM_1/crab_DStar_EnrichedMC_K3pi-Sep16-v3/150918_144531/0000/*.root");

//chain_Kpi->Add("~/mount/cms/store/user/scoopers/MinBias_TuneZ2star_13TeV-pythia6/DStar_MinBias_TuneZ2star-July21-v2-wHLT/150721_153230/0000/TrkAnalysis_MC_*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/MinBias_TuneEE5C_13TeV-herwigpp/DStar_MinBias_TuneEE5C-July21/150721_133708/0000/TrkAnalysis_MC_*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/MinBias_TuneCUETP8M1_13TeV-pythia8/DStar_MinBias_TuneCUETP8M1-July21-wHLT/150721_161541/0000/TrkAnalysis_MC_*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/MinBias_TuneMBR_13TeV-pythia8/DStar_MinBias_TuneMBR-July21/150721_133938/0000/TrkAnalysis_MC_*.root");

//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July17-v1/150717_154349/0000/*.root");
//chain_Kpi->Add("~/mount/cms/store/user/scoopers/ZeroBias/DStar_ZeroBias-July17-v1/150717_154349/0000/TrkAnalysis_generalTracks_9.root");
//chain_Kpi->Add("/afs/cern.ch/user/s/scoopers/mount/cms/store/user/scoopers/MinBias_TuneCUETP8M1_13TeV-pythia8/DStar_MinBias_TuneCUETP8M1-July21/150721_133237/0000/TrkAnalysis_MC_8.root");

//chain_Kpi->Add("AllData-Oct15-v2/ofile_tree_Job80.root");
//chain_K3pi->Add("AllData-Oct15-v2/ofile_tree_Job80.root");

//chain_Kpi->Print();

//TFile *ofile = new TFile(Form("%s/outputTree.dat",output_dir.c_str()),"RECREATE");
//TFile *ofile = new TFile("outputTree.dat","RECREATE");

//FILE *ofile1_Kpi = fopen("D0Mass_MC.dat","w");
FILE *ofile2_Kpi = fopen(Form("%s_dM.dat",name),"w");
//FILE *ofile1_K3pi = fopen("D0Mass_MC_K3pi.dat","w");
FILE *ofile2_K3pi = fopen(Form("%s_dM_K3pi.dat", name),"w");

TFile *ofile_tree = new TFile(Form("%s_ntuple.root",name),"RECREATE");
TDirectory *out_dir = ofile_tree->mkdir("analyzer");
out_dir->cd();
TTree *ofile_tree_Kpi;
TTree *ofile_tree_K3pi;


cout<<"copying tree"<<std::endl;
//TTree *otree = chain_K3pi->CopyTree("cosAlphaK3pi>0.99&&K3piTrkKpt>0.3&&K3piTrk1pipt>0.3&&K3piTrk2pipt>0.3&&K3piTrk3pipt>0.3&&K3piTrkKchi2<2.5&&K3piTrk1pichi2<2.5&&K3piTrk2pichi2<2.5&&K3piTrk3pichi2<2.5&&K3piTrkKnhits>=5&&K3piTrk1pinhits>=5&&K3piTrk2pinhits>=5&&K3piTrk3pinhits>=5&&K3piTrkKdz<1&&K3piTrk1pidz<1&&K3piTrk2pidz<1&&K3piTrk3pidz<1&&K3piTrkKdxy<0.1&&K3piTrk1pidxy<0.1&&K3piTrk2pidxy<0.1&&K3piTrk3pidxy<0.1&&abs(D0MassK3pi-1.864)<0.010&&D0VtxProb3>0.01&&K3piTrkSpt>0.25&&K3piTrkSchi2<3&&K3piTrkSnhits>2&&K3piTrkSdz<1&&K3piTrkSdxy<0.1&&DSPtK3pi>5.5");
//TTree *otree = chain_K3pi->CopyTree("cosAlphaK3pi>0.99&&K3piTrkKpt>0.5&&K3piTrk1pipt>0.5&&K3piTrk2pipt>0.5&&K3piTrk3pipt>0.5&&K3piTrkKchi2<2.5&&K3piTrk1pichi2<2.5&&K3piTrk2pichi2<2.5&&K3piTrk3pichi2<2.5&&K3piTrkKnhits>=5&&K3piTrk1pinhits>=5&&K3piTrk2pinhits>=5&&K3piTrk3pinhits>=5&&K3piTrkKdz<1&&K3piTrk1pidz<1&&K3piTrk2pidz<1&&K3piTrk3pidz<1&&K3piTrkKdxy<0.1&&K3piTrk1pidxy<0.1&&K3piTrk2pidxy<0.1&&K3piTrk3pidxy<0.1&&abs(D0MassK3pi-1.864)<0.012&&D0VtxProb3>0.03&&K3piTrkSpt>0.3&&K3piTrkSchi2<3&&K3piTrkSnhits>=2&&K3piTrkSdz<1&&K3piTrkSdxy<0.1&&DSPtK3pi>5.5");
TTree *otree_Kpi = (TTree*) chain_Kpi;
TTree *otree_K3pi = (TTree*) chain_K3pi;

cout<<"done copying tree"<<std::endl;
//otree->Print();
//std::vector<double> *cosAlphaKpi = 0;
//std::vector<double> *D0PtKpi = 0;
//std::vector<double> *D0etaKpi = 0;
//std::vector<double> *D0phiKpi = 0;
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
std::vector<double> *KpiTrkScharge = 0;
std::vector<double> *DSPtKpi = 0;
std::vector<double> *MCDsDeltaR = 0;
//std::vector<double> *D0VtxPosx = 0;
//std::vector<double> *D0VtxPosy = 0;
//std::vector<double> *D0VtxPosz = 0;
int NKpiCand;
double PVx;
double PVy;
double PVz;
double PVerrx;
double PVerry;
double PVerrz;
int PVOrder;

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
std::vector<double> *K3piTrkScharge = 0;
std::vector<double> *DSPtK3pi = 0;
std::vector<double> *flightLengthK3pi = 0;
std::vector<double> *MCDsDeltaR3 = 0;
int NK3piCand;

std::vector<double> *D0VtxPosx = 0;
std::vector<double> *D0VtxPosy = 0;
std::vector<double> *D0VtxPosz = 0;
std::vector<double> *D0Vtxerrx = 0;
std::vector<double> *D0Vtxerry = 0;
std::vector<double> *D0Vtxerrz = 0;
std::vector<double> *D0VtxPosx3 = 0;
std::vector<double> *D0VtxPosy3 = 0;
std::vector<double> *D0VtxPosz3 = 0;
std::vector<double> *D0Vtxerrx3 = 0;
std::vector<double> *D0Vtxerry3 = 0;
std::vector<double> *D0Vtxerrz3 = 0;

//otree->SetBranchAddress("cosAlphaKpi",&cosAlphaKpi, &b_cosAlphaKpi);
//otree->SetBranchAddress("D0PtKpi", &D0PtKpi, &b_D0PtKpi);
//otree->SetBranchAddress("D0etaKpi", &D0etaKpi, &b_D0etaKpi);
//otree->SetBranchAddress("D0phiKpi", &D0phiKpi, &b_D0phiKpi);
otree_Kpi->SetBranchAddress("D0MassKpi", &D0MassKpi);
otree_Kpi->SetBranchAddress("DSMassKpi", &DSMassKpi);
otree_Kpi->SetBranchAddress("cosAlphaKpi", &cosAlphaKpi);
otree_Kpi->SetBranchAddress("KpiTrkKpt", &KpiTrkKpt);
otree_Kpi->SetBranchAddress("KpiTrkpipt", &KpiTrkpipt);
otree_Kpi->SetBranchAddress("KpiTrkKchi2", &KpiTrkKchi2);
otree_Kpi->SetBranchAddress("KpiTrkpichi2", &KpiTrkpichi2);
otree_Kpi->SetBranchAddress("KpiTrkKnhits", &KpiTrkKnhits);
otree_Kpi->SetBranchAddress("KpiTrkpinhits", &KpiTrkpinhits);
otree_Kpi->SetBranchAddress("KpiTrkKdz", &KpiTrkKdz);
otree_Kpi->SetBranchAddress("KpiTrkpidz", &KpiTrkpidz);
otree_Kpi->SetBranchAddress("KpiTrkKdxy", &KpiTrkKdxy);
otree_Kpi->SetBranchAddress("KpiTrkpidxy", &KpiTrkpidxy);
otree_Kpi->SetBranchAddress("D0VtxProb", &D0VtxProb);
otree_Kpi->SetBranchAddress("KpiTrkSpt", &KpiTrkSpt);
otree_Kpi->SetBranchAddress("KpiTrkSchi2", &KpiTrkSchi2);
otree_Kpi->SetBranchAddress("KpiTrkSnhits", &KpiTrkSnhits);
otree_Kpi->SetBranchAddress("KpiTrkSdz", &KpiTrkSdz);
otree_Kpi->SetBranchAddress("KpiTrkSdxy", &KpiTrkSdxy);
otree_Kpi->SetBranchAddress("KpiTrkScharge", &KpiTrkScharge);
otree_Kpi->SetBranchAddress("DSPtKpi", &DSPtKpi);
otree_Kpi->SetBranchAddress("MCDsDeltaR", &MCDsDeltaR);
otree_Kpi->SetBranchAddress("PVx", &PVx);
otree_Kpi->SetBranchAddress("PVy", &PVy);
otree_Kpi->SetBranchAddress("PVz", &PVz);
otree_Kpi->SetBranchAddress("PVerrx", &PVerrx);
otree_Kpi->SetBranchAddress("PVerry", &PVerry);
otree_Kpi->SetBranchAddress("PVerrz", &PVerrz);
otree_Kpi->SetBranchAddress("D0VtxPosx", &D0VtxPosx);
otree_Kpi->SetBranchAddress("D0VtxPosy", &D0VtxPosy);
otree_Kpi->SetBranchAddress("D0VtxPosz", &D0VtxPosz);
otree_Kpi->SetBranchAddress("D0Vtxerrx", &D0Vtxerrx);
otree_Kpi->SetBranchAddress("D0Vtxerry", &D0Vtxerry);
otree_Kpi->SetBranchAddress("D0Vtxerrz", &D0Vtxerrz);
otree_Kpi->SetBranchAddress("PVOrder", &PVOrder);

//otree->SetBranchAddress("D0VtxPosx", &D0VtxPosx, &b_D0VtxPosx);
//otree->SetBranchAddress("D0VtxPosy", &D0VtxPosy, &b_D0VtxPosy);
//otree->SetBranchAddress("D0VtxPosz", &D0VtxPosz, &b_D0VtxPosz);

otree_Kpi->SetBranchAddress("NKpiCand", &NKpiCand);
//otree->SetBranchAddress("PVx", &PVx);
//otree->SetBranchAddress("PVy", &PVy);
//otree->SetBranchAddress("PVz", &PVz);

otree_K3pi->SetBranchAddress("D0MassK3pi", &D0MassK3pi);
otree_K3pi->SetBranchAddress("DSMassK3pi", &DSMassK3pi);
otree_K3pi->SetBranchAddress("cosAlphaK3pi", &cosAlphaK3pi);
otree_K3pi->SetBranchAddress("K3piTrkKpt", &K3piTrkKpt);
otree_K3pi->SetBranchAddress("K3piTrk1pipt", &K3piTrk1pipt);
otree_K3pi->SetBranchAddress("K3piTrk2pipt", &K3piTrk2pipt);
otree_K3pi->SetBranchAddress("K3piTrk3pipt", &K3piTrk3pipt);
otree_K3pi->SetBranchAddress("K3piTrkKchi2", &K3piTrkKchi2);
otree_K3pi->SetBranchAddress("K3piTrk1pichi2", &K3piTrk1pichi2);
otree_K3pi->SetBranchAddress("K3piTrk2pichi2", &K3piTrk2pichi2);
otree_K3pi->SetBranchAddress("K3piTrk3pichi2", &K3piTrk3pichi2);
otree_K3pi->SetBranchAddress("K3piTrkKnhits", &K3piTrkKnhits);
otree_K3pi->SetBranchAddress("K3piTrk1pinhits", &K3piTrk1pinhits);
otree_K3pi->SetBranchAddress("K3piTrk2pinhits", &K3piTrk2pinhits);
otree_K3pi->SetBranchAddress("K3piTrk3pinhits", &K3piTrk3pinhits);
otree_K3pi->SetBranchAddress("K3piTrkKdz", &K3piTrkKdz);
otree_K3pi->SetBranchAddress("K3piTrk1pidz", &K3piTrk1pidz);
otree_K3pi->SetBranchAddress("K3piTrk2pidz", &K3piTrk2pidz);
otree_K3pi->SetBranchAddress("K3piTrk3pidz", &K3piTrk3pidz);
otree_K3pi->SetBranchAddress("K3piTrkKdxy", &K3piTrkKdxy);
otree_K3pi->SetBranchAddress("K3piTrk1pidxy", &K3piTrk1pidxy);
otree_K3pi->SetBranchAddress("K3piTrk2pidxy", &K3piTrk2pidxy);
otree_K3pi->SetBranchAddress("K3piTrk3pidxy", &K3piTrk3pidxy);
otree_K3pi->SetBranchAddress("D0VtxProb3", &D0VtxProb3);
otree_K3pi->SetBranchAddress("K3piTrkSpt", &K3piTrkSpt);
otree_K3pi->SetBranchAddress("K3piTrkSchi2", &K3piTrkSchi2);
otree_K3pi->SetBranchAddress("K3piTrkSnhits", &K3piTrkSnhits);
otree_K3pi->SetBranchAddress("K3piTrkSdz", &K3piTrkSdz);
otree_K3pi->SetBranchAddress("K3piTrkSdxy", &K3piTrkSdxy);
otree_K3pi->SetBranchAddress("K3piTrkScharge", &K3piTrkScharge);
otree_K3pi->SetBranchAddress("DSPtK3pi", &DSPtK3pi);
otree_K3pi->SetBranchAddress("flightLengthK3pi", &flightLengthK3pi);
otree_K3pi->SetBranchAddress("NK3piCand", &NK3piCand);
otree_K3pi->SetBranchAddress("MCDsDeltaR3", &MCDsDeltaR3);
otree_K3pi->SetBranchAddress("PVx", &PVx);
otree_K3pi->SetBranchAddress("PVy", &PVy);
otree_K3pi->SetBranchAddress("PVz", &PVz);
otree_K3pi->SetBranchAddress("PVerrx", &PVerrx);
otree_K3pi->SetBranchAddress("PVerry", &PVerry);
otree_K3pi->SetBranchAddress("PVerrz", &PVerrz);
otree_K3pi->SetBranchAddress("D0VtxPosx3", &D0VtxPosx3);
otree_K3pi->SetBranchAddress("D0VtxPosy3", &D0VtxPosy3);
otree_K3pi->SetBranchAddress("D0VtxPosz3", &D0VtxPosz3);
otree_K3pi->SetBranchAddress("D0Vtxerrx3", &D0Vtxerrx3);
otree_K3pi->SetBranchAddress("D0Vtxerry3", &D0Vtxerry3);
otree_K3pi->SetBranchAddress("D0Vtxerrz3", &D0Vtxerrz3);
otree_K3pi->SetBranchAddress("PVOrder", &PVOrder);

ofile_tree_Kpi = otree_Kpi->CloneTree(0);
ofile_tree_K3pi = otree_K3pi->CloneTree(0);

double cosAlphaCut = 0.99;
double TrkPtCut = 0.5;
double Trkchi2Cut = 2.5;
double TrknhitsCut = 5.;
double TrkdzCut = 1.;
double TrkdxyCut = 0.1;
double massWindowCutKpi = 0.024;
double massWindowCutK3pi = 0.018;
double vtxProbCut = 0.03;
double TrkSptCut = 0.3;
double TrkSchi2Cut = 3.;
double TrkSnhitsCut = 2.;
double TrkSdzCut = 1.;
double TrkSdxyCut = 0.1;
double DSPtCut = 5.5;
double LSigCut = 3.0;

cout<<"Kpi"<<std::endl;
int nentries = otree_Kpi->GetEntries();
//int nentries = 10000;
//double divsize = (nentries*1.0)/nJobs;
cout<<nentries<<std::endl;
//int entry_first = (int) (divsize*jobNum);
//int entry_last = (int) entry_first + divsize; 
//cout<<"Running on entries "<<entry_first<<"through "<<entry_last<<std::endl;;
for (int i=0; i<nentries; i++) {
//for (int i=entry_first; i < entry_last; i++) {
    if ((i % 100000) == 0) { cout<<"processing entry: "<<i<<std::endl; }
    otree_Kpi->GetEntry(i);
    for (int j=0; j<NKpiCand; j++) {
        if ( (*cosAlphaKpi)[j] < cosAlphaCut)  continue;
        if ( (*KpiTrkKpt)[j]<TrkPtCut || (*KpiTrkpipt)[j]<TrkPtCut) continue;
        if ( (*KpiTrkKchi2)[j]>Trkchi2Cut || (*KpiTrkpichi2)[j]>Trkchi2Cut) continue;
        if ( (*KpiTrkKnhits)[j] < TrknhitsCut || (*KpiTrkpinhits)[j] < TrknhitsCut) continue;
        if ( (*KpiTrkKdz)[j]>TrkdzCut || (*KpiTrkpidz)[j]>TrkdzCut ) continue;
        if ( (*KpiTrkKdxy)[j]>TrkdxyCut || (*KpiTrkpidxy)[j]>TrkdxyCut ) continue;
        if ( fabs((*D0MassKpi)[j] - 1.864) > massWindowCutKpi ) continue;
        if ( (*D0VtxProb)[j] < vtxProbCut) continue;
        if ( (*KpiTrkSpt)[j] < TrkSptCut) continue; 
        if ( (*KpiTrkSchi2)[j] > TrkSchi2Cut) continue;
        if ( (*KpiTrkSnhits)[j] < TrkSnhitsCut) continue;
        if ( (*KpiTrkSdz)[j] > TrkSdzCut) continue;
        if ( (*KpiTrkSdxy)[j] > TrkSdxyCut) continue;
        if ( (*DSPtKpi)[j] < DSPtCut) continue;
        if (PVOrder != 1) continue;
        //if ( (*KpiTrkScharge)[j] <= 0) continue; // only take D*+ cand, to check directly with Valentina
        if ( ( 1.0 / ( sqrt( pow((*D0Vtxerrx)[j]/(*D0VtxPosx)[j],2) + pow((*D0Vtxerry)[j]/(*D0VtxPosy)[j],2) + pow(PVerrx/PVx,2) + pow(PVerry/PVy,2) ) ) ) < LSigCut) continue;
        //if ( (*MCDsDeltaR)[j] > 0.1) continue;
         //cout<<(*DSMassKpi)[j]<<std::endl;
        //cout<<"DM = "<<((*DSMassKpi)[j] - (*D0MassKpi)[j])<<std::endl;
        //fprintf(ofile1_Kpi,"%f\n",((*D0MassKpi)[j]));
        fprintf(ofile2_Kpi,"%f\n",((*DSMassKpi)[j] - (*D0MassKpi)[j]));
        ofile_tree_Kpi->Fill();
        //break;
        
    } 
}

cout<<"K3pi"<<std::endl;
nentries = otree_K3pi->GetEntries();
//divsize = (nentries*1.0)/nJobs;
//entry_first = (int) (divsize*jobNum);
//entry_last = (int) entry_first + divsize;
//nentries=10000;
cout<<nentries<<std::endl;
//cout<<"Running on entries "<<entry_first<<"through "<<entry_last<<std::endl;;
for (int i=0; i<nentries; i++) {
//for (int i=entry_first; i < entry_last; i++) {
    if ((i % 100000) == 0) { cout<<"processing entry: "<<i<<std::endl; }
    otree_K3pi->GetEntry(i);
    for (int j=0; j<NK3piCand; j++) {
        if ( (*cosAlphaK3pi)[j] < cosAlphaCut)  continue;
        if ( (*K3piTrkKpt)[j]<TrkPtCut || (*K3piTrk1pipt)[j]<TrkPtCut || (*K3piTrk2pipt)[j]<TrkPtCut || (*K3piTrk3pipt)[j]<TrkPtCut) continue;
        if ( (*K3piTrkKchi2)[j]>Trkchi2Cut || (*K3piTrk1pichi2)[j]>Trkchi2Cut || (*K3piTrk2pichi2)[j]>Trkchi2Cut || (*K3piTrk3pichi2)[j]>Trkchi2Cut) continue;
        if ( (*K3piTrkKnhits)[j]<TrknhitsCut || (*K3piTrk1pinhits)[j]<TrknhitsCut || (*K3piTrk2pinhits)[j]<TrknhitsCut || (*K3piTrk3pinhits)[j]<TrknhitsCut) continue;
        if ( (*K3piTrkKdz)[j]>TrkdzCut || (*K3piTrk1pidz)[j]>TrkdzCut || (*K3piTrk2pidz)[j]>TrkdzCut || (*K3piTrk3pidz)[j]>TrkdzCut) continue;
        if ( (*K3piTrkKdxy)[j]>TrkdxyCut || (*K3piTrk1pidxy)[j]>TrkdxyCut || (*K3piTrk2pidxy)[j]>TrkdxyCut || (*K3piTrk3pidxy)[j]>TrkdxyCut) continue;
        if ( fabs((*D0MassK3pi)[j] - 1.864) > massWindowCutK3pi ) continue;
        if ( (*D0VtxProb3)[j] < vtxProbCut) continue;
        if ( (*K3piTrkSpt)[j] < TrkSptCut) continue;
        if ( (*K3piTrkSchi2)[j] > TrkSchi2Cut) continue;
        if ( (*K3piTrkSnhits)[j] < TrkSnhitsCut) continue;
        if ( (*K3piTrkSdz)[j] > TrkSdzCut) continue;
        if ( (*K3piTrkSdxy)[j] > TrkSdxyCut) continue;
        if ( (*DSPtK3pi)[j] < DSPtCut) continue;
        if (PVOrder != 1) continue;
        //if ( (*K3piTrkScharge)[j] <= 0) continue; // only take D*+ cand, to check directly with Valentina
        if ( ( 1.0 / ( sqrt( pow((*D0Vtxerrx3)[j]/(*D0VtxPosx3)[j],2) + pow((*D0Vtxerry3)[j]/(*D0VtxPosy3)[j],2) + pow(PVerrx/PVx,2) + pow(PVerry/PVy,2) ) ) ) < LSigCut) continue;
        //if ( (*MCDsDeltaR3)[j] > 0.1) continue;
        //fprintf(ofile1_K3pi,"%f\n",((*D0MassK3pi)[j]));
        fprintf(ofile2_K3pi,"%f\n",((*DSMassK3pi)[j] - (*D0MassK3pi)[j])); 
        ofile_tree_K3pi->Fill();
        //break;
    }
}

ofile_tree->cd();
out_dir->cd();
ofile_tree_Kpi->Write("tree1");
ofile_tree_K3pi->Write("tree2");
ofile_tree->Close();
//fclose(ofile1_Kpi);
fclose(ofile2_Kpi);
//fclose(ofile1_K3pi);
fclose(ofile2_K3pi);

}
