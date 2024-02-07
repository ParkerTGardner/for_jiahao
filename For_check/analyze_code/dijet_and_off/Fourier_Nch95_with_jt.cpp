
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
#include "../include/sherpa_constants_t.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "math.h"
#include "TGraphErrors.h" 
#include "TLine.h"
void processRootFile() {
   
    // TFile* file = new TFile("newroot/dijet/AB_data_3.root");
    // TFile* file2 = new TFile("newroot/dijet/AB_data_7-20.root");
    // TFile* file3 = new TFile("newroot/dijet/AB_MC_2.root");
    // TFile* file4 = new TFile("newroot/dijet/AB_MC_7-20.root");

    TFile* file = new TFile("newroot/dijet/AA_data.root");
    TFile* file2 = new TFile("newroot/dijet/AA_data_7-20.root");
    TFile* file3 = new TFile("newroot/dijet/AA_MC_2.root");
    TFile* file4 = new TFile("newroot/dijet/AA_MC_7-20.root");

    TGraphErrors* g1 = new TGraphErrors(5);
    TGraphErrors* g2 = new TGraphErrors(5);
    TGraphErrors* g3 = new TGraphErrors(5);

    TGraphErrors* g12 = new TGraphErrors(5);
    TGraphErrors* g22 = new TGraphErrors(5);
    TGraphErrors* g32 = new TGraphErrors(5);



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    // Nch>95 evolve with jt, AA and AB,  MC and data


    TString v_for_pt = Form("Singlejet Fourier_components_for_Nch>95 for 2<#Delta#eta*<6");
    // Variables to hold the v1, v2, v3 values
    TCanvas* c1 = new TCanvas("c1", v_for_pt, 2400, 1800);

    double DelEta_lo = 2.0;
    double DelEta_mid = 6.0;
    double DelEta_hi  =  6.0;

    
    int point = 0;
    int point_2 = 0;
    int point_3 = 0;
    int wtrk = trackbin;
    for (int wppt = 1; wppt < 4; wppt++) {
        
        TH2D* hSig = (TH2D*)file->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
        TH2D* hBack = (TH2D*)file->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
        TH1D* hBin = (TH1D*)file->Get(Form("hBinDist_gen_%d",wtrk));
        hSig->SetDefaultSumw2(kTRUE);
        hBack->SetDefaultSumw2(kTRUE);
        

        TH1D *hSig1D = (TH1D*) hSig->ProjectionY("h1D", hSig->GetXaxis()->FindBin(DelEta_lo), hSig->GetXaxis()->FindBin(DelEta_mid))->Clone();
        TH1D *hBack1D = (TH1D*) hBack->ProjectionY("h1D", hBack->GetXaxis()->FindBin(DelEta_lo), hBack->GetXaxis()->FindBin(DelEta_mid))->Clone();

        hSig1D->SetDefaultSumw2(kTRUE);
        hBack1D->SetDefaultSumw2(kTRUE);
        hSig1D->Divide(hBack1D);
        hSig1D->Scale(hBack1D->GetMaximum());
        TH1D* h1D = hSig1D;
        h1D->SetDefaultSumw2(kTRUE);
        
        std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)))";
        TF1* func = new TF1("func", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
        func->SetParameter(0, h1D->GetMaximum());
        func->SetParameter(1, 0.1);
        func->SetParameter(2, 0.1);
        func->SetParameter(3, 0.1);
        func->SetParameter(4, 0.1);
        h1D->Fit(func, "q 0");
        h1D->Fit(func, "q 0");
        h1D->Fit(func, "m q 0");
        h1D->Fit(func, "m q 0");
        h1D->Fit(func, "m q E 0");
        h1D->Fit(func, "m E q 0");
        



        g1->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(1));
        g1->SetPointError(point, 0, func->GetParError(1));
        g2->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(2));
        g2->SetPointError(point, 0, func->GetParError(2));
        g3->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(3));
        g3->SetPointError(point, 0, func->GetParError(3));
        point++;


    }

    

    for (int wppt = 4; wppt < 6; wppt++) {
        
        
        TH2D* hSig = (TH2D*)file2->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
        TH2D* hBack = (TH2D*)file2->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
        TH1D* hBin = (TH1D*)file2->Get(Form("hBinDist_gen_%d",wtrk));
        hSig->SetDefaultSumw2(kTRUE);
        hBack->SetDefaultSumw2(kTRUE);

        

        TH1D *hSig1D = (TH1D*) hSig->ProjectionY("h1D", hSig->GetXaxis()->FindBin(DelEta_lo), hSig->GetXaxis()->FindBin(DelEta_mid))->Clone();
        TH1D *hBack1D = (TH1D*) hBack->ProjectionY("h1D", hBack->GetXaxis()->FindBin(DelEta_lo), hBack->GetXaxis()->FindBin(DelEta_mid))->Clone();
        hSig1D->SetDefaultSumw2(kTRUE);
        hBack1D->SetDefaultSumw2(kTRUE);
        hSig1D->Divide(hBack1D);
        hSig1D->Scale(hBack1D->GetMaximum());
        TH1D* h1D = hSig1D;
        h1D->SetDefaultSumw2(kTRUE);
        
        std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)))";
        TF1* func = new TF1("func", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
        

        func->SetParameter(0, h1D->GetMaximum());
        func->SetParameter(1, 0.1);
        
        func->SetParameter(2, 0.1);
        func->SetParameter(3, 0.1);
        func->SetParameter(4, 0.1);
        h1D->Fit(func, "q 0");
        h1D->Fit(func, "q 0");
        h1D->Fit(func, "m q 0");
        h1D->Fit(func, "m q 0");
        h1D->Fit(func, "m q E 0");
        h1D->Fit(func, "m E q 0");
        



        g1->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(1));
        g1->SetPointError(point, 0, func->GetParError(1));
        g1->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(1));
        g2->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(2));
        g2->SetPointError(point, 0, func->GetParError(2));
        g3->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(3));
        g3->SetPointError(point, 0, func->GetParError(3));
        point++;


    }


    for (int wppt = 1; wppt < 4; wppt++) {
        
        TH2D* hSig = (TH2D*)file3->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
        TH2D* hBack = (TH2D*)file3->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
        TH1D* hBin = (TH1D*)file3->Get(Form("hBinDist_gen_%d",wtrk));
        hSig->SetDefaultSumw2(kTRUE);
        hBack->SetDefaultSumw2(kTRUE);
        

        TH1D *hSig1D = (TH1D*) hSig->ProjectionY("h1D", hSig->GetXaxis()->FindBin(DelEta_lo), hSig->GetXaxis()->FindBin(DelEta_mid))->Clone();
        TH1D *hBack1D = (TH1D*) hBack->ProjectionY("h1D", hBack->GetXaxis()->FindBin(DelEta_lo), hBack->GetXaxis()->FindBin(DelEta_mid))->Clone();

        hSig1D->SetDefaultSumw2(kTRUE);
        hBack1D->SetDefaultSumw2(kTRUE);
        hSig1D->Divide(hBack1D);
        hSig1D->Scale(hBack1D->GetMaximum());
        TH1D* h1D = hSig1D;
        h1D->SetDefaultSumw2(kTRUE);
        
        std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)))";
        TF1* func = new TF1("func", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
        func->SetParameter(0, h1D->GetMaximum());
        func->SetParameter(1, 0.1);
        func->SetParameter(2, 0.1);
        func->SetParameter(3, 0.1);
        func->SetParameter(4, 0.1);
        h1D->Fit(func, "q 0");
        h1D->Fit(func, "q 0");
        h1D->Fit(func, "m q 0");
        h1D->Fit(func, "m q 0");
        h1D->Fit(func, "m q E 0");
        h1D->Fit(func, "m E q 0");
        



        g12->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(1));
        g12->SetPointError(point, 0, func->GetParError(1));
        g22->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(2));
        g22->SetPointError(point, 0, func->GetParError(2));
        g32->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(3));
        g32->SetPointError(point, 0, func->GetParError(3));
        point_2++;


    }

    

    for (int wppt = 4; wppt < 6; wppt++) {
        
        
        TH2D* hSig = (TH2D*)file4->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
        TH2D* hBack = (TH2D*)file4->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
        TH1D* hBin = (TH1D*)file4->Get(Form("hBinDist_gen_%d",wtrk));
        hSig->SetDefaultSumw2(kTRUE);
        hBack->SetDefaultSumw2(kTRUE);

        

        TH1D *hSig1D = (TH1D*) hSig->ProjectionY("h1D", hSig->GetXaxis()->FindBin(DelEta_lo), hSig->GetXaxis()->FindBin(DelEta_mid))->Clone();
        TH1D *hBack1D = (TH1D*) hBack->ProjectionY("h1D", hBack->GetXaxis()->FindBin(DelEta_lo), hBack->GetXaxis()->FindBin(DelEta_mid))->Clone();
        hSig1D->SetDefaultSumw2(kTRUE);
        hBack1D->SetDefaultSumw2(kTRUE);
        hSig1D->Divide(hBack1D);
        hSig1D->Scale(hBack1D->GetMaximum());
        TH1D* h1D = hSig1D;
        h1D->SetDefaultSumw2(kTRUE);
        
        std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)))";
        TF1* func = new TF1("func", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
        

        func->SetParameter(0, h1D->GetMaximum());
        func->SetParameter(1, 0.1);
        
        func->SetParameter(2, 0.1);
        func->SetParameter(3, 0.1);
        func->SetParameter(4, 0.1);
        h1D->Fit(func, "q 0");
        h1D->Fit(func, "q 0");
        h1D->Fit(func, "m q 0");
        h1D->Fit(func, "m q 0");
        h1D->Fit(func, "m q E 0");
        h1D->Fit(func, "m E q 0");
        



        g12->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(1));
        g12->SetPointError(point, 0, func->GetParError(1));
        g12->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(1));
        g22->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(2));
        g22->SetPointError(point, 0, func->GetParError(2));
        g32->SetPoint(point, ptbinbounds_lo[wppt-1], func->GetParameter(3));
        g32->SetPointError(point, 0, func->GetParError(3));
        point_2++;


    }



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






            
    TLegend* legend = new TLegend(0.7,0.75,0.9,0.9); // Adjust these parameters to move the legend to the desired position




    legend->AddEntry(g1, "V1_data", "p");
    legend->AddEntry(g2, "V2_data", "p");
    legend->AddEntry(g3, "V3_data", "p");
    legend->AddEntry(g12, "V1_MC", "l");
    legend->AddEntry(g22, "V2_MC", "l");
    legend->AddEntry(g32, "V3_MC", "l");

   



    double ymin = -0.1;
    double ymax = 0.1;


    TH1F* frame = c1->DrawFrame(-0.1, ymin, 1.1, ymax);
    frame->GetXaxis()->SetTitle("jT (GeV/c)");
    frame->SetTitle(v_for_pt);
    
    g1->SetMarkerColor(kRed);
    g1->SetMarkerSize(5);
    g1->SetMarkerStyle(21); 
    g1->Draw("PE");

    g2->SetMarkerColor(kBlue);
    g2->SetMarkerSize(5);
    g2->SetMarkerStyle(34); 
    g2->Draw("PE same");

    g3->SetMarkerColor(kGreen);
    g3->SetMarkerSize(5);
    g3->SetMarkerStyle(47); 
    g3->Draw("PE same");



    g12->SetLineColor(kRed);
    g12->SetLineWidth(5);
    g12->SetLineStyle(9); 
    g12->Draw("CE same");

    g22->SetLineColor(kBlue);
    g22->SetLineWidth(5);
    g22->SetLineStyle(9); 
    g22->Draw("CE same");

    g32->SetLineColor(kGreen);
    g32->SetLineWidth(5);
    g32->SetLineStyle(9); 
    g32->Draw("CE same");


    double x_min = -0.1;
    // double x_max = g1->GetXaxis()->GetXmax();
    double x_max = 1.1;

    TLine* line1 = new TLine(x_min, 0, x_max, 0);
    line1->SetLineColor(kBlack);
    // line1->SetLineStyle(9);
    line1->SetLineWidth(3);
    line1->Draw("same");

    TLine* line2 = new TLine(x_min, -0.0025, x_max, -0.0025);
    line2->SetLineColor(kBlue);
    // line2->SetLineStyle(9);
    line2->SetLineWidth(3);
    line2->Draw("same");

    TLine* line3 = new TLine(x_min, -0.01, x_max, -0.01);
    line3->SetLineColor(kGreen);
    // line3->SetLineStyle(9);
    line3->SetLineWidth(3);
    line3->Draw("same");

    TLine* line4 = new TLine(x_min, -0.0225, x_max, -0.0225);
    line4->SetLineColor(kViolet);
    // line4->SetLineStyle(9);
    line4->SetLineWidth(3);
    line4->Draw("same");

    TLine* line5 = new TLine(x_min, -0.04, x_max, -0.04);
    line5->SetLineColor(kCyan+3);
    // line5->SetLineStyle(9);
    line5->SetLineWidth(3);
    line5->Draw("same");


    TLine* line6 = new TLine(x_min, 0.0025, x_max, 0.0025);
    line6->SetLineColor(kBlue);
    // line2->SetLineStyle(9);
    line6->SetLineWidth(3);
    line6->Draw("same");

    TLine* line7 = new TLine(x_min, 0.01, x_max, 0.01);
    line7->SetLineColor(kGreen);
    // lne3->SetLineStyle(9);
    line7->SetLineWidth(3);
    line7->Draw("same");

    TLine* line8 = new TLine(x_min, 0.0225, x_max, 0.0225);
    line8->SetLineColor(kViolet);
    // line4->SetLineStyle(9);
    line8->SetLineWidth(3);
    line8->Draw("same");

    TLine* line9 = new TLine(x_min, 0.04, x_max, 0.04);
    line9->SetLineColor(kCyan+3);
    // line5->SetLineStyle(9);
    line9->SetLineWidth(3);
    line9->Draw("same");
    
            
          


    legend->Draw();
    c1->Update();

    gSystem->mkdir("output_buffer/fouriers/dijet", kTRUE);

    TString filename = Form("output_buffer/fouriers/dijet/Singlejet_Fourier_components_for_Nch>95_26.png");
    // std::cout << "Saving canvas to: " << filename << std::endl;
    c1->SaveAs(filename);
    delete c1;
    // delete g1;
    // delete g2;
    // delete g3;
    
    
    
    
}

int main() {
    processRootFile();
    return 0;
}



