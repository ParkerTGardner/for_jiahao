#include <iostream>
#include <fstream>
#include "math.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TKey.h"
#include "TClass.h"
#include "TROOT.h"
#include "../include/sherpa_constants.h"


void Fourier_analysis() {
    // Open the file
    TFile* file = new TFile("newroot/dijet/AB_MC_2.root");
    // Loop over all objects in the file
    TIter next(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        TH2D* h = (TH2D*)key->ReadObj();
         
        // // Check if the histogram name matches the pattern
        std::string name = h->GetName();
        if (name.find("hSignalS") != 0) continue;

        // // Extract the parameters from the histogram name
        int wtrk=0, wppt=0, wpPU=0;
        sscanf(name.c_str(), "hSignalS_trk_%d_ppt_%d_PU_%d", &wtrk, &wppt, &wpPU);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //  main code


        TH2D* hSig = (TH2D*)file->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
        TH2D* hBack = (TH2D*)file->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
        
        hSig->SetDefaultSumw2(kTRUE);
        hBack->SetDefaultSumw2(kTRUE);


        // here I use the 1.6 deltaeta cut, equals to no cut for dijet, single jet should be 2.0

        double DelEta_lo = 1.6;
        TH1D *hSig1D = (TH1D*) hSig->ProjectionY("h1D", hSig->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
        TH1D *hBack1D = (TH1D*) hBack->ProjectionY("h1D", hBack->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
        hSig1D->SetDefaultSumw2(kTRUE);
        hBack1D->SetDefaultSumw2(kTRUE);
        hSig1D->Divide(hBack1D);
        hSig1D->Scale(hBack->GetMaximum());
        TH1D* h1D = hSig1D;
        h1D->SetDefaultSumw2(kTRUE);
        
        std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)))";
        TF1* func = new TF1("func", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
        

        func->SetParameter(0, h1D->GetMaximum());
        func->SetParameter(1, 0);
        func->SetParameter(2, 0);
        func->SetParameter(3, 0);
        h1D->Fit(func, "q 0");
        h1D->Fit(func, "q 0");
        h1D->Fit(func, "m q 0");
        h1D->Fit(func, "m q 0");
        h1D->Fit(func, "m q E 0");
        h1D->Fit(func, "m E q 0");


        // Get the parameters
        double v0 = func->GetParameter(0);
        double v1 = func->GetParameter(1);
        double v2 = func->GetParameter(2);
        double v3 = func->GetParameter(3);

        h1D->Scale(M_PI/v0);
        func->SetParameter(0,M_PI);

        // Create the functions
        TF1* func1 = new TF1("v1", "[0]*TMath::Cos(x)+0.5", h1D->GetXaxis()->GetXmin(), h1D->GetXaxis()->GetXmax());
        func1->SetParameter(0, v1);
        TF1* func2 = new TF1("v2", "[0]*TMath::Cos(2*x)+0.5", h1D->GetXaxis()->GetXmin(), h1D->GetXaxis()->GetXmax());
        func2->SetParameter(0, v2);
        TF1* func3 = new TF1("v3", "[0]*TMath::Cos(3*x)+0.5", h1D->GetXaxis()->GetXmin(), h1D->GetXaxis()->GetXmax());
        func3->SetParameter(0, v3);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////      
        

        TCanvas* c1 = new TCanvas("c1", "Fourier_fit", 1600, 1200);
        TString graphtitle = Form("AB_MC %d_to_%d_pt_%.2f_to_%.2f.png", trackbinbounds[wtrk-1], trackbinboundsUpper[wtrk-1], ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
        
    
        h1D->SetStats(false);
        h1D->SetTitle(graphtitle);
        h1D->Draw();
        
        h1D->GetXaxis()->SetTitle("#Delta#phi^{*}");
        func->SetLineColor(kBlack);
        func->Draw("same");
        
        h1D->GetYaxis()->SetRangeUser(0.465, 0.56);
        // Draw func1, func2, and func3
        func1->SetLineColor(kRed);
        func1->SetLineWidth(2);
        func1->Draw("same");
        func2->SetLineColor(kGreen);
        func2->SetLineWidth(2);
        func2->Draw("same");
        func3->SetLineColor(kBlue);
        func2->SetLineWidth(2);
        func3->Draw("same");

       


        // Create a legend
        TLegend* leg = new TLegend(0.15,0.75,0.35,0.9);
        leg->AddEntry(h1D,"Original h1D","l");
        leg->AddEntry(func,"Fit function","l");
        leg->AddEntry(func1,"V1*cos(x)","l");
        leg->AddEntry(func2,"V2*cos(2x)","l");
        leg->AddEntry(func3,"V3*cos(3x)","l");
        leg->Draw();

        // Create the folder if it doesn't exist
        gSystem->mkdir("output_buffer/Fourier_fit_curves/AB_MC", kTRUE);

        TString filename = "output_buffer/Fourier_fit_curves/AB_MC/"+ graphtitle;
        c1->SaveAs(filename);
        // Clean up
        delete c1;
        delete func;
        delete func1;
        delete func2;
        delete func3;
        delete h1D;
        delete hSig;
        delete hBack;
    }
    delete file;
}

int main() {
    Fourier_analysis();
    return 0;
}

