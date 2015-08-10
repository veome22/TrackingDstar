#define final4_cxx
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooDstD0BG.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h" 
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooAbsReal.h"
#include "RooGlobalFunc.h"
#include "RooChebychev.h"
#include "RooGenericPdf.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooRealVar.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooGaussian.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooExponential.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooPolynomial.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooDstD0BG.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooCBShape.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooAddPdf.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooDataSet.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooGlobalFunc.h" 
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooPlot.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooFitResult.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooExtendPdf.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooAbsReal.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooGlobalFunc.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooChebychev.h"
//#include "/uscmst1/prod/sw/cms/slc5_ia32_gcc434/lcg/roofit/5.27.06-cms3/include/RooGenericPdf.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TROOT.h>
#include <fstream>
#include <iostream>
#include <algorithm>
//#include "/uscms/home/mrmooney/setTDRStyle.C"
#include "TGraphErrors.h"
#include "TPaveText.h"

using namespace RooFit;

class final4 {
public :

   final4();
   virtual ~final4();
   virtual double doFit(bool usePixel, bool isMC);
};

final4::final4()
{
}

final4::~final4()
{
}

double final4::doFit(bool usePixel, bool isMC)
{
  gROOT->SetBatch(kTRUE);
  gROOT->SetStyle("Plain");

  //setTDRStyle();

  RooRealVar x("","",0.140,0.160);
  x.SetTitle("M(K#pi#pi_{S}) - M(K#pi) [GeV/c^{ 2}]");

  RooRealVar mean("mean", "mean", 0.145556,0.144, 0.147);
  //RooRealVar mean("mean", "mean", 0.145556);
  RooRealVar sigma("sigma", "sigma", 0.00069, 0.0002, 0.0013);
  //RooRealVar sigma("sigma", "sigma", 0.00069);
  RooGaussian gauss("gauss","gaussian PDF", x, mean, sigma);

  RooRealVar alpha("alpha", "alpha", -1.0, -10.0, 10.0);
  RooRealVar power("power", "power", 3.0, 0.0, 50.0);
  RooCBShape cball("cball", "crystal ball PDF", x, mean, sigma, alpha, power);

  RooRealVar dm0("dm0", "dm0", 0.13957);
  dm0.setConstant(kTRUE);  
  RooRealVar shape("shape","shape",0.,-100.,100.);
  RooRealVar dstp1("p1","p1",0.,-500.,500.);
  RooRealVar dstp2("p2","p2",0.,-500.,500.);

  shape.setRange(0.000001,10.0);//was 0.02
  shape.setVal(0.0017);
  dstp1.setVal(0.45);
  dstp2.setVal(13.0);

  RooDstD0BG bkg("bkg","bkg",x,dm0,shape,dstp1,dstp2);

  RooRealVar c0("c0","c0",10.0,-10.0,11.0);
  RooRealVar c1("c1","c1",10.0,-10.0,11.0);
  RooRealVar c2("c2","c2",10.0,-10.0,11.0);
  RooRealVar c3("c3","c3",10.0,-10.0,11.0);
  RooRealVar c4("c4","c4",10.0,-10.0,11.0);
  RooRealVar c5("c5","c5",10.0,-10.0,11.0);
  RooRealVar c6("c6","c6",10.0,-10.0,11.0);
  RooRealVar c7("c7","c7",10.0,-10.0,11.0);
  RooRealVar c8("c8","c8",10.0,-10.0,11.0);
  RooGenericPdf cutoff("cutoff","cutoff","(@0 > @1)*(@2*abs(@0-@1) + @3*pow(abs(@0-@1),2) + @4*pow(abs(@0-@1),3) + @5*pow(abs(@0-@1),4) + @6*pow(abs(@0-@1),5) + @7*pow(abs(@0-@1),6) + @8*pow(abs(@0-@1),7))",RooArgSet(x,dm0,c0,c1,c2,c3,c4,c5,c6));

  RooRealVar poly1("poly1","poly1",1.0,-5000.0,5000.0);
  RooRealVar poly2("poly2","poly2",1.0,-5000.0,5000.0);
  RooRealVar poly3("poly3","poly3",1.0,-5000.0,5000.0);
  RooRealVar poly4("poly4","poly4",1.0,-5000.0,5000.0);

  RooPolynomial polybkg("polybkg","polybkg",x,RooArgSet(poly1,poly2,poly3,poly4));

  RooRealVar cheby0("cheby0","cheby0",1.0,-500.0,500.0);
  RooRealVar cheby1("cheby1","cheby1",1.0,-500.0,500.0);
  RooRealVar cheby2("cheby2","cheby2",1.0,-500.0,500.0);
  RooRealVar cheby3("cheby3","cheby3",1.0,-500.0,500.0);

  RooChebychev chebybkg("chebybkg","chebybkg",x,RooArgSet(cheby0,cheby1,cheby2,cheby3));

  RooRealVar mean2("mean2", "mean2", 0.14548,0.144, 0.147);
  RooRealVar sigma2("sigma2", "sigma2", 0.00065, 0.0002, 0.005);
  RooGaussian gauss2("gauss2","gaussian PDF 2", x, mean2, sigma2);

  RooRealVar S("S", "Signal Yield", 2000, 0, 300000);
  //RooRealVar S("S", "Signal Yield", 0, 0, 300000);
  RooRealVar SS("SS", "Signal Yield #2", 100, 0, 100000);
  RooRealVar S2("S2", "Signal2 Yield (MC only)", 0, 0, 200);
  RooRealVar B("B", "Background Yield", 40000, 0, 30000000);
  //RooRealVar B("B", "Background Yield", 0, 0, 30000000);

  RooAddPdf sum("sum", "gaussian plus threshold PDF",RooArgList(gauss, bkg), RooArgList(S, B));
  //RooAddPdf sum("sum", "background PDF",RooArgList(polybkg), RooArgList(B));
  //RooAddPdf sum("sum", "background PDF",RooArgList(chebybkg), RooArgList(B));
  //RooAddPdf sum("sum", "background PDF",RooArgList(cutoff), RooArgList(B));
  // RooAddPdf sum("sum", "gaussians plus threshold PDF",RooArgList(gauss, gauss2, bkg), RooArgList(S, SS, B));
  //RooAddPdf sum("sum", "crystal ball plus threshold PDF",RooArgList(cball, bkg), RooArgList(S, B));
  RooAddPdf sumMC("sumMC","double gaussian",RooArgList(gauss, gauss2), RooArgList(S, S2));

  fstream file;

  char filename[50];
  double cut=5.5;
  sprintf(filename,"dM.dat");
  
  RooDataSet* data = RooDataSet::read(filename,RooArgList(x));
  RooFitResult* fit = 0;
  if (isMC == 0)
  {
    fit = sum.fitTo(*data,RooFit::Extended(),PrintLevel(1),Save(true),RooFit::NumCPU(8),RooFit::Strategy(2));
    file << "cut: " << cut << "GeV" << endl;
    file << "status: " << fit->status() << endl;
    file << "covQual: " << fit->covQual() << endl;
    file << "edm: " << fit->edm() <<  endl;
    file << "Yield: " <<  S.getVal() << " " << S.getError() << endl;
    file << "Bkg: " << B.getVal() << " " << B.getError() << endl; 
    file << "sigma: " << sigma.getVal() <<  " " << sigma.getError() << endl;
    file << "mean: " << mean.getVal() << " " << mean.getError() << endl;
    file << "shape: " << shape.getVal() << " " << shape.getError() << endl;
    file << "dstp1: " << dstp1.getVal() << " " << dstp1.getError() << endl;
    file << "dstp2: " << dstp2.getVal() << " " << dstp2.getError() << endl;
    file << endl;
  }
  else
  {
    fit = sumMC.fitTo(*data,RooFit::Extended(),PrintLevel(1),Save(true),RooFit::NumCPU(8),RooFit::Strategy(2),Range(0.142,0.15));
    file << "cut: " << cut << "GeV" << endl;
    file << "status: " << fit->status() << endl;
    file << "covQual: " << fit->covQual() << endl;
    file << "edm: " << fit->edm() <<  endl;
    file << "Yield: " <<  S.getVal() << " " << S.getError() << endl;
    file << "Yield2: " << S2.getVal() << " " << S2.getError() << endl; 
    file << "sigma: " << sigma.getVal() <<  " " << sigma.getError() << endl;
    file << "mean: " << mean.getVal() << " " << mean.getError() << endl;
    file << "sigma2: " << sigma2.getVal() << " " << sigma2.getError() << endl;
    file << "mean2: " << mean2.getVal() << " " << mean2.getError() << endl;
    file << endl; 
  }
  
  RooPlot* xFrame = x.frame(Bins(40));
  xFrame->SetTitle("D* #rightarrow D^{0}(K#pi)#pi");
  data->plotOn(xFrame);
  if(isMC == 0)
  {
    sum.plotOn(xFrame);
    sum.plotOn(xFrame,RooFit::Components(bkg),RooFit::LineStyle(kDashed));
  }
  else
  {
    sumMC.plotOn(xFrame, Range(0.139,0.159));
    //      sumMC.plotOn(xFrame,RooFit::Components(gauss2),RooFit::LineStyle(kDashed),Range(0.139,0.159));
  }
  data->plotOn(xFrame);
  file << xFrame->chiSquare() << endl;
  TCanvas c;
  TPaveText* ptext = 0;
  TPaveText* ptex = 0;
  
  if(usePixel == 0)
  {
    ptext = new TPaveText(0.47,0.30,0.9,0.40,"TRNDC");
    ptex = new TPaveText(0.47,0.15,0.9,0.25,"NDC"); 
  }
  else
  {
    if(isMC == 0)
    {
      ptext = new TPaveText(0.47,0.8,0.9,0.9,"TRNDC");
      ptex = new TPaveText(0.47,0.8,0.9,0.9,"NDC"); 
    }
    else
    {
      ptext = new TPaveText(0.47,0.78,0.9,0.88,"TRNDC");
      ptex = new TPaveText(0.47,0.63,0.9,0.73,"NDC"); 
    }
  }
  
  ptext->SetFillColor(0);
  ptext->SetTextSize(0.04);
  ptext->SetTextAlign(13);
  ptext->AddText("CMS Preliminary");
  //ptext->AddText("#sqrt{s} = 13 TeV, 40.0 pb^{-1}");
  ptext->AddText("#sqrt{s} = 13 TeV, Spring15 MinBias MC");
  xFrame->SetYTitle("Events / 0.5 MeV/c^{ 2}");
  xFrame->GetYaxis()->SetTitleOffset(1.3);
  xFrame->GetYaxis()->SetLabelSize(0.03);
  xFrame->GetXaxis()->SetLabelSize(0.03);
  
  ptex->SetFillColor(0);
  ptex->SetTextSize(0.033);
  ptex->SetTextAlign(13);
  char theyield[50];
  char themean[50];
  char thesigma[50];
  if(isMC == 0)
  {
    sprintf(theyield,"Yield = %u #pm %u",(unsigned)S.getVal(),(unsigned)S.getError());
    sprintf(themean,"Mean = (%.3f #pm %.3f) MeV/c^{2}",mean.getVal()*1000.0,mean.getError()*1000.0);
    sprintf(thesigma,"Sigma = (%.3f #pm %.3f) MeV/c^{2}",sigma.getVal()*1000.0,sigma.getError()*1000.0);
  }
  else
  {
    sprintf(theyield,"Yield = %u #pm %u",(unsigned)(S.getVal()+S2.getVal()),(unsigned)(sqrt(pow(S.getError(),2)+pow(S2.getError(),2))));
    sprintf(themean,"Mean = (%.3f #pm %.3f) MeV/c^{2}",mean.getVal()*1000.0,mean.getError()*1000.0);
    sprintf(thesigma,"Sigma = (%.3f #pm %.3f) MeV/c^{2}",sigma.getVal()*1000.0,sigma.getError()*1000.0);
  } 
  ptex->AddText(theyield);
  ptex->AddText(themean);
  ptex->AddText(thesigma);
  
  xFrame->Draw();
  ptext->Draw("same");
  ptex->Draw("same");
  
  c.SaveAs("dM.png");
  c.SaveAs("dM.pdf");

  S.Print();

  // compute integrals
  double sfactor;
  if(isMC)
    sfactor = 0.0;
  else
  {
    double sigbkg, bkg1, bkg2;
    double base = 0.145421;
    
    x.setRange("signal",base-0.0013,base+0.0013);
    x.setRange("background1",base-0.0051,base-0.0025);
    x.setRange("background2",base+0.0025,base+0.0051);    

    RooAbsReal* iSB = B.createIntegral(x,Range("signal"));
    RooAbsReal* iB1 = B.createIntegral(x,Range("background1"));
    RooAbsReal* iB2 = B.createIntegral(x,Range("background2"));

    sigbkg = iSB->getVal();
    bkg1 = iB1->getVal();
    bkg2 = iB2->getVal();
 
    sfactor = -1.0 * (sigbkg / (bkg1 + bkg2));
  }

  return sfactor;
}

