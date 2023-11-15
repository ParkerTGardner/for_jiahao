#ifndef FOURIER_FIT
#define FOURIER_FIT

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TTree.h"
#include "TF1.h"
#include "TH2.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TKey.h"
#include "TClass.h"
#include "TROOT.h"
#include "TGraph.h" 
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "math.h"
#include "TGraphErrors.h" 
#include "TLatex.h"
#include "TLine.h"
#include "TMinuit.h"

TF1* fourier_fit(TH2D* hSig, TH2D* hBack, double DelEta_lo = 2.0, double V1_init = -0.1, double V2_init = 0.0225, double V3_init = -0.01) {

    TH1D *hSig1D = (TH1D*) hSig->ProjectionY("h1D", hSig->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
    TH1D *hBack1D = (TH1D*) hBack->ProjectionY("h1D", hBack->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
    hSig1D->SetDefaultSumw2(kTRUE);
    hBack1D->SetDefaultSumw2(kTRUE);
    
    hSig1D->Divide(hBack1D);
    hSig1D->Scale(hBack1D->GetMaximum());
    TH1D* h1D = hSig1D;
    h1D->SetDefaultSumw2(kTRUE);
    
    std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";
    TF1* func = new TF1("func", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
    
    
    func->SetParameter(0, h1D->GetMaximum());
    func->SetParameter(1, V1_init);
    
    func->SetParameter(2, V2_init);
    func->SetParameter(3, V3_init);
    
    h1D->Fit(func, "q 0");
    h1D->Fit(func, "q 0");
    h1D->Fit(func, "m q 0");
    h1D->Fit(func, "m q 0");
    h1D->Fit(func, "m q E 0");
    h1D->Fit(func, "m E q 0");
    return func;
}


















































#endif