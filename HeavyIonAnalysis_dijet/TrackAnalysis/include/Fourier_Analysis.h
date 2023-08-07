//
// Created by Jiahao on 2023/7/26.
//Fourier_Analysis
#ifndef FOURIER_ANALYSIS
#define FOURIER_ANALYSIS


#include <iostream>
#include <fstream>
#include "math.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "tuple"
 


void FitHarmonic( TH2D* hSig, TH2D* hBack, double v[4][][], int wtrk, int wppt, int wpPU) {
    double DelEta*_lo = 2.0;
   
    hSig->Divide(hBack);
    // Create a 1D projection for DeltaEta > 2
    TH1D* h1D = hSig->ProjectionY("h1D", hSig->GetXaxis()->FindBin(DelEta*_lo), -1);
    
    // Define a 3 term harmonic function
    TF1* func = new TF1("func", "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))", h1D->GetXaxis()->GetXmin(), h1D->GetXaxis()->GetXmax());
    func->SetParameter(0, h1D->GetMaximum());
    // Fit the function to the histogram
    h1D->Fit(func, "Q");
    v[0][wtrk-1][wppt*wpPU-1] = func->GetParameter(0);
    v[1][wtrk-1][wppt*wpPU-1] = func->GetParameter(1);
    v[2][wtrk-1][wppt*wpPU-1] = func->GetParameter(2);
    v[3][wtrk-1][wppt*wpPU-1] = func->GetParameter(3); 

    TF1* func1 = new TF1("v1", "v[1][wtrk-1][wppt*wpPU-1] *cos(x) ", h1D->GetXaxis()->GetXmin(), h1D->GetXaxis()->GetXmax());
    TF1* func2 = new TF1("v2", "v[2][wtrk-1][wppt*wpPU-1] *cos(x) ", h1D->GetXaxis()->GetXmin(), h1D->GetXaxis()->GetXmax()); 
    TF1* func3 = new TF1("v3", "v[3][wtrk-1][wppt*wpPU-1] *cos(x) ", h1D->GetXaxis()->GetXmin(), h1D->GetXaxis()->GetXmax()); 
    // Draw the histogram and the fit
    TCanvas* c1 = new TCanvas("c1", "Fourier_fit", 800, 600);
    h1D->Draw();
    func->Draw("same");
    func1->Draw("same");
    func2->Draw("same");
    func3->Draw("same");
    //Save the canvas
    c1->SaveAs(Form("Fourier_fit_for_pt:%2f_to_%2f.png",wtrk,wppt*wpPU));
    // Clean up
    delete c1;
    delete func;
    delete func1;
    delete func2;
    delete func3; 
    delete h1D;
    delete hSig;
    delete hBack;
    // delete file;
    // Open the file and get the histograms
    // TFile* file = new TFile("job_all.root");
    // TH2D* hSig = (TH2D*)file->Get("hSigS_70_to_1000_and_3_to_30_w_PU_1");
    // TH2D* hBack = (TH2D*)file->Get("hBckS_70_to_1000_and_3_to_30_w_PU_1");

    // Divide the signal histogram by the background histogram
    // TH2D* hSig2D = (TH2D*)hSig->Clone("hSig2D");
    // TH2D* hBack2D = (TH2D*)hBack->Clone("hBack2D");

    // // Clear the bins for DeltaEta <= 2
    // for (int i = 1; i <=hSig2D->GetXaxis()->FindBin(DelEta*_lo) ; i++) {
    //     for (int j = 1; j <= hSig2D->GetNbinsY(); j++) {
    //         hSig2D->SetBinContent(i, j, 0);
    //         hBack2D->SetBinContent(i, j, 0);
    //     }
    // }
    // hSig2D->Divide(hBack2D);
}
#endif