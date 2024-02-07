
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
#include "../include/sherpa_constants.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "math.h"
#include "TGraphErrors.h" 
#include "TLine.h"
void processRootFile() {
   
    TFile* file = new TFile("newroot/dijet/AB_data_3.root");
    TFile* file2 = new TFile("newroot/dijet/AB_MC_2.root");
    TFile* file3 = new TFile("newroot/dijet/AB_data_3.root");
    TFile* file4 = new TFile("newroot/dijet/AA_data.root");
    TFile* file5 = new TFile("newroot/dijet/AA_MC_2.root");
    TFile* file6 = new TFile("newroot/dijet/AA_data.root");

    // TFile* file3 = new TFile("newroot/MultAB_A/trkbin_8/AB_A_trig_2.root");

    
   
        for (int wppt = 1; wppt <ptbin+1 ; wppt++) {
            
            

            // Create a graph for each array
            TGraphErrors* g1   = new TGraphErrors(trackbin-1);
            TGraphErrors* g1_2 = new TGraphErrors(trackbin-1);
            TGraphErrors* g1_3 = new TGraphErrors(trackbin-1);
            TGraphErrors* g1_4 = new TGraphErrors(trackbin-1);
            TGraphErrors* g1_5 = new TGraphErrors(trackbin-1);
            TGraphErrors* g1_6 = new TGraphErrors(trackbin-1);

            TGraphErrors* g2   = new TGraphErrors(trackbin-1);
            TGraphErrors* g2_2 = new TGraphErrors(trackbin-1);
            TGraphErrors* g2_3 = new TGraphErrors(trackbin-1);
            TGraphErrors* g2_4 = new TGraphErrors(trackbin-1);
            TGraphErrors* g2_5 = new TGraphErrors(trackbin-1);
            TGraphErrors* g2_6 = new TGraphErrors(trackbin-1);
            
            TGraphErrors* g3   = new TGraphErrors(trackbin-1);
            TGraphErrors* g3_2 = new TGraphErrors(trackbin-1);
            TGraphErrors* g3_3 = new TGraphErrors(trackbin-1);
            TGraphErrors* g3_4 = new TGraphErrors(trackbin-1);
            TGraphErrors* g3_5 = new TGraphErrors(trackbin-1);
            TGraphErrors* g3_6 = new TGraphErrors(trackbin-1);


            int point = 0;
            int point_2 = 0;
            int point_3 = 0;
            int point_4 = 0;
            int point_5 = 0;
            int point_6 = 0;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



            //This code is for AA, AB, comparison
            //It is AB_data, AB_MC, AB_data_bck; AA_data, AB_MC, AA_data_bck



            for (int wtrk = 2; wtrk < trackbin+1; wtrk++) {

                
                TH2D* hSig = (TH2D*)file->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH2D* hBack = (TH2D*)file->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH1D* hBin = (TH1D*)file->Get(Form("hBinDist_gen_%d",wtrk));
                
                hSig->SetDefaultSumw2(kTRUE);
                hBack->SetDefaultSumw2(kTRUE);

                double DelEta_lo_di = 1.6;
                double DelEta_lo_sing = 2.0;
                TH1D *hSig1D = (TH1D*) hSig->ProjectionY("h1D", hSig->GetXaxis()->FindBin(DelEta_lo_di), -1)->Clone();

                TH1D *hBack1D = (TH1D*) hBack->ProjectionY("h1D", hBack->GetXaxis()->FindBin(DelEta_lo_di), -1)->Clone();
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
                h1D->Fit(func, "q 0");
                h1D->Fit(func, "q 0");
                h1D->Fit(func, "m q 0");
                h1D->Fit(func, "m q 0");
                h1D->Fit(func, "m q E 0");
                h1D->Fit(func, "m E q 0");
                


                g1->SetPoint(point, hBin->GetMean(1), func->GetParameter(1));
                g1->SetPointError(point, 0, func->GetParError(1));
        
                                
                g2->SetPoint(point, hBin->GetMean(1), func->GetParameter(2));
                g2->SetPointError(point, 0, func->GetParError(2));
                // std::cout << "err for AB_data v2: " << point << ": (" << hBin->GetMean(1) << ", " << func->GetParError(2) << ")" << std::endl;
                g3->SetPoint(point, hBin->GetMean(1), func->GetParameter(3));
                g3->SetPointError(point, 0, func->GetParError(3));
                point++;




                TH2D* hSig2 = (TH2D*)file2->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH2D* hBack2 = (TH2D*)file2->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH1D* hBin_2 = (TH1D*)file2->Get(Form("hBinDist_gen_%d",wtrk));
                
                hSig2->SetDefaultSumw2(kTRUE);
                hBack2->SetDefaultSumw2(kTRUE);

                TH1D *hSig1D_2 = (TH1D*) hSig2->ProjectionY("h1D", hSig2->GetXaxis()->FindBin(DelEta_lo_di), -1)->Clone();
                TH1D *hBack1D_2 = (TH1D*) hBack2->ProjectionY("h1D", hBack2->GetXaxis()->FindBin(DelEta_lo_di), -1)->Clone();
                hSig1D_2->SetDefaultSumw2(kTRUE);
                hBack1D_2->SetDefaultSumw2(kTRUE);
                hSig1D_2->Divide(hBack1D_2);
                hSig1D_2->Scale(hBack1D_2->GetMaximum());
                TH1D* h1D_2 = hSig1D_2;
                h1D_2->SetDefaultSumw2(kTRUE);
                
                
                TF1* func_2 = new TF1("func_2", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                func_2->SetParameter(0, h1D_2->GetMaximum());
                func_2->SetParameter(1, 0.1);
               
                func_2->SetParameter(2, 0.1);
                func_2->SetParameter(3, 0.1);
                h1D_2->Fit(func_2, "q 0");
                h1D_2->Fit(func_2, "q 0");
                h1D_2->Fit(func_2, "m q 0");
                h1D_2->Fit(func_2, "m q 0");
                h1D_2->Fit(func_2, "m q E 0");
                h1D_2->Fit(func_2, "m E q 0");
              
                g1_2->SetPoint(point_2, hBin_2->GetMean(1), func_2->GetParameter(1));
                g1_2->SetPointError(point_2, 0, func_2->GetParError(1));




                g2_2->SetPoint(point_2, hBin_2->GetMean(1), func_2->GetParameter(2));
                g2_2->SetPointError(point_2, 0, func_2->GetParError(2));
                g3_2->SetPoint(point_2, hBin_2->GetMean(1), func_2->GetParameter(3));
                g3_2->SetPointError(point_2, 0, func_2->GetParError(3));
                point_2++;

               

               

                // TH2D* hSig3 = (TH2D*)file3->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH2D* hBack3 = (TH2D*)file3->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH1D* hBin_3 = (TH1D*)file3->Get(Form("hBinDist_gen_%d",wtrk));
                hBack3->SetDefaultSumw2(kTRUE);

                TH1D *hBack1D_3 = (TH1D*) hBack3->ProjectionY("h1D", hBack3->GetXaxis()->FindBin(DelEta_lo_di), -1)->Clone();
            
                hBack1D_3->SetDefaultSumw2(kTRUE);
                hBack1D_3->Scale(hBack1D_3->GetMaximum());
                TH1D* h1D_3 = hBack1D_3;
                h1D_3->SetDefaultSumw2(kTRUE);
                
                
                TF1* func_3 = new TF1("func_3", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                func_3->SetParameter(0, h1D_3->GetMaximum());
                func_3->SetParameter(1, 0.1);
               
                func_3->SetParameter(2, 0.1);
                func_3->SetParameter(3, 0.1);
                h1D_3->Fit(func_3, "q 0");
                h1D_3->Fit(func_3, "q 0");
                h1D_3->Fit(func_3, "m q 0");
                h1D_3->Fit(func_3, "m q 0");
                h1D_3->Fit(func_3, "m q E 0");
                h1D_3->Fit(func_3, "m E q 0");
         



                g1_3->SetPoint(point_3, hBin_3->GetMean(1), func_3->GetParameter(1));
                g1_3->SetPointError(point_3, 0, func_3->GetParError(1));
                g2_3->SetPoint(point_3, hBin_3->GetMean(1), func_3->GetParameter(2));
                g2_3->SetPointError(point_3, 0, func_3->GetParError(2));
                g3_3->SetPoint(point_3, hBin_3->GetMean(1), func_3->GetParameter(3));
                g3_3->SetPointError(point_3, 0, func_3->GetParError(3));
                point_3++;


                

                

                TH2D* hSig4 = (TH2D*)file4->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH2D* hBack4 = (TH2D*)file4->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH1D* hBin_4 = (TH1D*)file4->Get(Form("hBinDist_gen_%d",wtrk));

                
                hSig4->SetDefaultSumw2(kTRUE);
                hBack4->SetDefaultSumw2(kTRUE);

                
                
                TH1D *hSig1D_4 = (TH1D*) hSig4->ProjectionY("h1D", hSig4->GetXaxis()->FindBin(DelEta_lo_sing), -1)->Clone();
                TH1D *hBack1D_4 = (TH1D*) hBack4->ProjectionY("h1D", hBack4->GetXaxis()->FindBin(DelEta_lo_sing), -1)->Clone();

                
            
                hSig1D_4->SetDefaultSumw2(kTRUE);
                hBack1D_4->SetDefaultSumw2(kTRUE);
                hSig1D_4->Divide(hBack1D_4);
                hSig1D_4->Scale(hBack1D_4->GetMaximum());
                TH1D* h1D_4 = hSig1D_4;
                h1D_4->SetDefaultSumw2(kTRUE);

                
                
                
                TF1* func_4 = new TF1("func_4", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                func_4->SetParameter(0, h1D_4->GetMaximum());
                func_4->SetParameter(1, 0.0);
               
                func_4->SetParameter(2, 0.0);
                func_4->SetParameter(3, 0.0);
                h1D_4->Fit(func_4, "q 0");
                h1D_4->Fit(func_4, "q 0");
                h1D_4->Fit(func_4, "m q 0");
                h1D_4->Fit(func_4, "m q 0");
                h1D_4->Fit(func_4, "m q E 0");
                h1D_4->Fit(func_4, "m E q 0");

         



                g1_4->SetPoint(point_4, hBin_4->GetMean(1), func_4->GetParameter(1));
                g1_4->SetPointError(point_4, 0, func_4->GetParError(1));


                g2_4->SetPoint(point_4, hBin_4->GetMean(1), func_4->GetParameter(2));
                g2_4->SetPointError(point_4, 0, func_4->GetParError(2));
                g3_4->SetPoint(point_4, hBin_4->GetMean(1), func_4->GetParameter(3));
                g3_4->SetPointError(point_4, 0, func_4->GetParError(3));
                point_4++;


                TH2D* hSig5 = (TH2D*)file5->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH2D* hBack5 = (TH2D*)file5->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH1D* hBin_5 = (TH1D*)file5->Get(Form("hBinDist_gen_%d",wtrk));
                hSig5->SetDefaultSumw2(kTRUE);
                hBack5->SetDefaultSumw2(kTRUE);

                TH1D *hSig1D_5 = (TH1D*) hSig5->ProjectionY("h1D", hSig5->GetXaxis()->FindBin(DelEta_lo_sing), -1)->Clone();
                TH1D *hBack1D_5 = (TH1D*) hBack5->ProjectionY("h1D", hBack5->GetXaxis()->FindBin(DelEta_lo_sing), -1)->Clone();
            
                hSig1D_5->SetDefaultSumw2(kTRUE);
                hBack1D_5->SetDefaultSumw2(kTRUE);
                hSig1D_5->Divide(hBack1D_5);
                hSig1D_5->Scale(hBack1D_5->GetMaximum());
                TH1D* h1D_5 = hSig1D_5;
                h1D_5->SetDefaultSumw2(kTRUE);
                
                
                TF1* func_5 = new TF1("func_5", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                func_5->SetParameter(0, h1D_5->GetMaximum());
                func_5->SetParameter(1, 0.1);
               
                func_5->SetParameter(2, 0.1);
                func_5->SetParameter(3, 0.1);
                h1D_5->Fit(func_5, "q 0");
                h1D_5->Fit(func_5, "q 0");
                h1D_5->Fit(func_5, "m q 0");
                h1D_5->Fit(func_5, "m q 0");
                h1D_5->Fit(func_5, "m q E 0");
                h1D_5->Fit(func_5, "m E q 0");
         



                g1_5->SetPoint(point_5, hBin_5->GetMean(1), func_5->GetParameter(1));
                g1_5->SetPointError(point_5, 0, func_5->GetParError(1));
                g2_5->SetPoint(point_5, hBin_5->GetMean(1), func_5->GetParameter(2));
                g2_5->SetPointError(point_5, 0, func_5->GetParError(2));
                g3_5->SetPoint(point_5, hBin_5->GetMean(1), func_5->GetParameter(3));
                g3_5->SetPointError(point_5, 0, func_5->GetParError(3));
                point_5++;






                // TH2D* hSig6 = (TH2D*)file6->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH2D* hBack6 = (TH2D*)file6->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH1D* hBin_6 = (TH1D*)file6->Get(Form("hBinDist_gen_%d",wtrk));
                
                hBack6->SetDefaultSumw2(kTRUE);

            
                TH1D *hBack1D_6 = (TH1D*) hBack6->ProjectionY("h1D", hBack6->GetXaxis()->FindBin(DelEta_lo_sing), -1)->Clone();
            
                hBack1D_6->SetDefaultSumw2(kTRUE);
                hBack1D_6->Scale(hBack1D_6->GetMaximum());
                TH1D* h1D_6 = hBack1D_6;
                h1D_6->SetDefaultSumw2(kTRUE);
                
                
                TF1* func_6 = new TF1("func_6", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                func_6->SetParameter(0, h1D_6->GetMaximum());
                func_6->SetParameter(1, 0.1);
               
                func_6->SetParameter(2, 0.1);
                func_6->SetParameter(3, 0.1);                // std::cout<< "fit3 before"<< std::endl;
                h1D_6->Fit(func_6, "q 0");
                h1D_6->Fit(func_6, "q 0");
                h1D_6->Fit(func_6, "m q 0");
                h1D_6->Fit(func_6, "m q 0");
                h1D_6->Fit(func_6, "m q E 0");
                h1D_6->Fit(func_6, "m E q 0");
         



                g1_6->SetPoint(point_6, hBin_6->GetMean(1), func_6->GetParameter(1));
                g1_6->SetPointError(point_6, 0, func_6->GetParError(1));
                g2_6->SetPoint(point_6, hBin_6->GetMean(1), func_6->GetParameter(2));
                g2_6->SetPointError(point_6, 0, func_6->GetParError(2));
                g3_6->SetPoint(point_6, hBin_5->GetMean(1), func_6->GetParameter(3));
                g3_6->SetPointError(point_6, 0, func_6->GetParError(3));
                point_6++;

            }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            
            TString v1_for_pt = Form("V1_for_pt:%.2f_to_%.2f",ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
            // Variables to hold the v1, v2, v3 values
            TCanvas* c1 = new TCanvas("c1", v1_for_pt, 2400, 1800);
            
            TLegend* legend = new TLegend(0.6,0.65,0.9,0.9); 


            

            //max error of each curve





            int n1 = g1->GetN();  // Number of points in the graph
            double *xErrors = g1->GetEX();  // Array of x errors
            double *yErrors = g1->GetEY();  // Array of y errors
            
            int n1_2 = g1_2->GetN();  // Number of points in the graph
            double *xErrors1_2 = g1_2->GetEX();  // Array of x errors
            double *yErrors1_2 = g1_2->GetEY();  // Array of y errors

            int n1_3 = g1_3->GetN();  // Number of points in the graph
            double *xErrors1_3 = g1_3->GetEX();  // Array of x errors
            double *yErrors1_3 = g1_3->GetEY();  // Array of y errors

            int n1_4 = g1_4->GetN();  // Number of points in the graph
            double *xErrors1_4 = g1_4->GetEX();  // Array of x errors
            double *yErrors1_4 = g1_4->GetEY();  // Array of y errors

            int n1_5 = g1_5->GetN();  // Number of points in the graph
            double *xErrors1_5 = g1_5->GetEX();  // Array of x errors
            double *yErrors1_5 = g1_5->GetEY();  // Array of y errors

            int n1_6 = g1_6->GetN();  // Number of points in the graph
            double *xErrors1_6 = g1_6->GetEX();  // Array of x errors
            double *yErrors1_6 = g1_6->GetEY();  // Array of y errors

            double ymaxError1_1 = *std::max_element(yErrors, yErrors + n1);
            double ymaxError1_2 = *std::max_element(yErrors1_2, yErrors1_2 + n1_2);
            double ymaxError1_3 = *std::max_element(yErrors1_3, yErrors1_3 + n1_3);
            double ymaxError1_4 = *std::max_element(yErrors1_4, yErrors1_4 + n1_4);
            double ymaxError1_5 = *std::max_element(yErrors1_5, yErrors1_5 + n1_5);
            double ymaxError1_6 = *std::max_element(yErrors1_6, yErrors1_6 + n1_6);

            legend->AddEntry(g1, Form("V1_AB_data,err:%.3e",ymaxError1_1), "p");
            legend->AddEntry(g1_2, Form("V1_AB_MC,err:%.3e", ymaxError1_2), "l");
            legend->AddEntry(g1_3, Form("V1_AB_bck,err:%.3e", ymaxError1_3), "l");
            legend->AddEntry(g1_4, Form("V1_AA_data,err:%.3e", ymaxError1_4), "p");
            legend->AddEntry(g1_5, Form("V1_AA_MC,err:%.3e", ymaxError1_5), "l");
            legend->AddEntry(g1_6, Form("V1_AA_bck,err:%.3e", ymaxError1_6), "l");


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
           

            double ymin = -0.3;
            double ymax = 0.3;

            g1->GetYaxis()->SetRangeUser(ymin, ymax);
            g1_2->GetYaxis()->SetRangeUser(ymin, ymax);
            g1_3->GetYaxis()->SetRangeUser(ymin, ymax);
            g1_4->GetYaxis()->SetRangeUser(ymin, ymax);
            g1_5->GetYaxis()->SetRangeUser(ymin, ymax);
            g1_6->GetYaxis()->SetRangeUser(ymin, ymax);

            
            
            

            g1->SetTitle(v1_for_pt);
            g1->GetXaxis()->SetTitle("Multiplicity");
            g1->SetMarkerColor(kBlue);
            g1->SetMarkerSize(5);
            g1->SetMarkerStyle(21); 
            g1->Draw("APE");

            g1_2->SetLineColor(kBlue);
            g1_2->SetLineWidth(5);
            g1_2->SetLineStyle(9);
            g1_2->Draw("CPE same");

            g1_3->SetLineWidth(3);
            g1_3->SetLineColor(kGreen);
            g1_3->Draw("CPE same");

            g1_4->SetMarkerColor(kRed);
            g1_4->SetMarkerSize(5);
            g1_4->SetMarkerStyle(34); 
            g1_4->Draw("PE same");

            g1_5->SetLineColor(kRed);
            g1_5->SetLineWidth(5);
            g1_5->SetLineStyle(9);
            g1_5->Draw("CPE same");

            g1_6->SetLineColor(kViolet);
            g1_6->Draw("CPE same");

            legend->Draw();
            c1->Update();

            gSystem->mkdir("output_buffer/fouriers/dijet_16", kTRUE);
            TString filename = Form("output_buffer/fouriers/dijet_16/v1_for_pt_%.2f_to_%.2f.png",  ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
            c1->SaveAs(filename);
            delete c1;




            TString v2_for_pt = Form("V2_for_pt:%.2f_to_%.2f",ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
            // Variables to hold the v1, v2, v3 values
            TCanvas* c2 = new TCanvas("c2", v2_for_pt, 2400, 1800);
            
            TLegend* legend2 = new TLegend(0.6,0.65,0.9,0.9); // Adjust these parameters to move the legend to the desired position


            int n2 = g2->GetN();  // Number of points in the graph
            double *xErrors2 = g2->GetEX();  // Array of x errors
            double *yErrors2 = g2->GetEY();  // Array of y errors
            
            int n2_2 = g2_2->GetN();  // Number of points in the graph
            double *xErrors2_2 = g2_2->GetEX();  // Array of x errors
            double *yErrors2_2 = g2_2->GetEY();  // Array of y errors

            int n2_3 = g2_3->GetN();  // Number of points in the graph
            double *xErrors2_3 = g2_3->GetEX();  // Array of x errors
            double *yErrors2_3 = g2_3->GetEY();  // Array of y errors

            int n2_4 = g2_4->GetN();  // Number of points in the graph
            double *xErrors2_4 = g2_4->GetEX();  // Array of x errors
            double *yErrors2_4 = g2_4->GetEY();  // Array of y errors

            int n2_5 = g2_5->GetN();  // Number of points in the graph
            double *xErrors2_5 = g2_5->GetEX();  // Array of x errors
            double *yErrors2_5 = g2_5->GetEY();  // Array of y errors

            int n2_6 = g2_6->GetN();  // Number of points in the graph
            double *xErrors2_6 = g2_6->GetEX();  // Array of x errors
            double *yErrors2_6 = g2_6->GetEY();  // Array of y errors

            double ymaxError2_1 = *std::max_element(yErrors2, yErrors2 + n2);
            double ymaxError2_2 = *std::max_element(yErrors2_2, yErrors2_2 + n2_2);
            double ymaxError2_3 = *std::max_element(yErrors2_3, yErrors2_3 + n2_3);
            double ymaxError2_4 = *std::max_element(yErrors2_4, yErrors2_4 + n2_4);
            double ymaxError2_5 = *std::max_element(yErrors2_5, yErrors2_5 + n2_5);
            double ymaxError2_6 = *std::max_element(yErrors2_6, yErrors2_6 + n2_6);
            
            legend2->AddEntry(g2, Form("V2_AB_data,err:%.3e",ymaxError2_1), "p");
            legend2->AddEntry(g2_2, Form("V2_AB_MC,err:%.3e", ymaxError2_2), "l");
            legend2->AddEntry(g2_3, Form("V2_AB_bck,err:%.3e", ymaxError2_3), "l");
            legend2->AddEntry(g2_4, Form("V2_AA_data,err:%.3e", ymaxError2_4), "p");
            legend2->AddEntry(g2_5, Form("V2_AA_MC,err:%.3e", ymaxError2_5), "l");
            legend2->AddEntry(g2_6, Form("V2_AA_bck,err:%.3e", ymaxError2_6), "l");
            

            double ymin2 = -0.05;
            double ymax2 = 0.15;

            g2->GetYaxis()->SetRangeUser(ymin2, ymax2);
            g2_2->GetYaxis()->SetRangeUser(ymin2, ymax2);
            g2_3->GetYaxis()->SetRangeUser(ymin2, ymax2);
            g2_4->GetYaxis()->SetRangeUser(ymin2, ymax2);
            g2_5->GetYaxis()->SetRangeUser(ymin2, ymax2);
            g2_6->GetYaxis()->SetRangeUser(ymin2, ymax2);
            

            g2->SetTitle(v2_for_pt);
            g2->GetXaxis()->SetTitle("Multiplicity");
            g2->SetMarkerColor(kBlue);
            g2->SetMarkerSize(5);
            g2->SetMarkerStyle(21); 
            g2->Draw("APE");

            g2_2->SetLineColor(kBlue);
            g2_2->SetLineWidth(5);
            g2_2->SetLineStyle(9);
            g2_2->Draw("CPE same");

            g2_3->SetLineWidth(3);
            g2_3->SetLineColor(kGreen);
            g2_3->Draw("CPE same");

            g2_4->SetMarkerColor(kRed);
            g2_4->SetMarkerSize(5);
            g2_4->SetMarkerStyle(34); 
            g2_4->Draw("PE same");

            g2_5->SetLineColor(kRed);
            g2_5->SetLineWidth(5);
            g2_5->SetLineStyle(9);
            g2_5->Draw("CPE same");

            g2_6->SetLineColor(kViolet);
            g2_6->Draw("CPE same");

            double x_min = g2->GetXaxis()->GetXmin();
            // double x_max = g1->GetXaxis()->GetXmax();
            double x_max = g2->GetX()[g2->GetN()-1];

            TLine* line1 = new TLine(x_min, 0, x_max, 0);
            line1->SetLineColor(kBlack);
            line1->SetLineStyle(9);
            line1->Draw();

            TLine* line2 = new TLine(x_min, 0.0025, x_max, 0.0025);
            line2->SetLineColor(kBlue);
            line2->SetLineStyle(9);
            line2->Draw();

            TLine* line3 = new TLine(x_min, 0.01, x_max, 0.01);
            line3->SetLineColor(kGreen);
            line3->SetLineStyle(9);
            line3->Draw();

            TLine* line4 = new TLine(x_min, 0.0225, x_max, 0.0225);
            line4->SetLineColor(kViolet);
            line4->SetLineStyle(9);
            line4->Draw();

            TLine* line5 = new TLine(x_min, 0.04, x_max, 0.04);
            line5->SetLineColor(kCyan+3);
            line5->SetLineStyle(9);
            line5->Draw();

            legend2->Draw();
            c2->Update();

            // gSystem->mkdir("output_buffer/fouriers/trkbin_4/Fourier_fit_each_Mult_20_20_CM_Pair", kTRUE);
            TString filename2 = Form("output_buffer/fouriers/dijet_16/v2_for_pt_%.2f_to_%.2f.png",  ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
            c2->SaveAs(filename2);
            delete c2;






            TString v3_for_pt = Form("V3_for_pt:%.2f_to_%.2f",ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
            // Variables to hold the v1, v2, v3 values
            TCanvas* c3 = new TCanvas("c3", v3_for_pt, 2400, 1800);
            
            TLegend* legend3 = new TLegend(0.6,0.65,0.9,0.9); // Adjust these parameters to move the legend to the desired position

            
            int n3 = g3->GetN();  // Number of points in the graph
            double *xErrors3 = g3->GetEX();  // Array of x errors
            double *yErrors3 = g3->GetEY();  // Array of y errors
            
            int n3_2 = g3->GetN();  // Number of points in the graph
            double *xErrors3_2 = g3_2->GetEX();  // Array of x errors
            double *yErrors3_2 = g3_2->GetEY();  // Array of y errors

            int n3_3 = g3_3->GetN();  // Number of points in the graph
            double *xErrors3_3 = g3_3->GetEX();  // Array of x errors
            double *yErrors3_3 = g3_3->GetEY();  // Array of y errors

            int n3_4 = g3_4->GetN();  // Number of points in the graph
            double *xErrors3_4 = g3_4->GetEX();  // Array of x errors
            double *yErrors3_4 = g3_4->GetEY();  // Array of y errors

            int n3_5 = g3_5->GetN();  // Number of points in the graph
            double *xErrors3_5 = g3_5->GetEX();  // Array of x errors
            double *yErrors3_5 = g3_5->GetEY();  // Array of y errors

            int n3_6 = g3_6->GetN();  // Number of points in the graph
            double *xErrors3_6 = g3_6->GetEX();  // Array of x errors
            double *yErrors3_6 = g3_6->GetEY();  // Array of y errors

            double ymaxError3_1 = *std::max_element(yErrors3, yErrors3 + n3);
            double ymaxError3_2 = *std::max_element(yErrors3_2, yErrors3_2 + n3_2);
            double ymaxError3_3 = *std::max_element(yErrors3_3, yErrors3_3 + n3_3);
            double ymaxError3_4 = *std::max_element(yErrors3_4, yErrors3_4 + n3_4);
            double ymaxError3_5 = *std::max_element(yErrors3_5, yErrors3_5 + n3_5);
            double ymaxError3_6 = *std::max_element(yErrors3_6, yErrors3_6 + n3_6);
            
            legend3->AddEntry(g3, Form("V3_AB_data,err:%.3e",ymaxError3_1), "p");
            legend3->AddEntry(g3_2, Form("V3_AB_MC,err:%.3e", ymaxError3_2), "l");
            legend3->AddEntry(g3_3, Form("V3_AB_bck,err:%.3e", ymaxError3_3), "l");
            legend3->AddEntry(g3_4, Form("V3_AA_data,err:%.3e", ymaxError3_4), "p");
            legend3->AddEntry(g3_5, Form("V3_AA_MC,err:%.3e", ymaxError3_5), "l");
            legend3->AddEntry(g3_6, Form("V3_AA_bck,err:%.3e", ymaxError3_6), "l");
            
            

            double ymin3 = -0.05;
            double ymax3 = 0.05;

            g3->GetYaxis()->SetRangeUser(ymin3, ymax3);
            g3_2->GetYaxis()->SetRangeUser(ymin3, ymax3);
            g3_3->GetYaxis()->SetRangeUser(ymin3, ymax3);
            g3_4->GetYaxis()->SetRangeUser(ymin3, ymax3);
            g3_5->GetYaxis()->SetRangeUser(ymin3, ymax3);
            g3_6->GetYaxis()->SetRangeUser(ymin3, ymax3);
            

            g3->SetTitle(v3_for_pt);
            g3->GetXaxis()->SetTitle("Multiplicity");
            g3->SetMarkerColor(kBlue);
            g3->SetMarkerSize(5);
            g3->SetMarkerStyle(21); 
            g3->Draw("APE");

            g3_2->SetLineColor(kBlue);
            g3_2->SetLineWidth(5);
            g3_2->SetLineStyle(9);
            g3_2->Draw("CPE same");

            g3_3->SetLineWidth(3);
            g3_3->SetLineColor(kGreen);
            g3_3->Draw("CPE same");

            g3_4->SetMarkerColor(kRed);
            g3_4->SetMarkerSize(5);
            g3_4->SetMarkerStyle(34); 
            g3_4->Draw("PE same");

            g3_5->SetLineColor(kRed);
            g3_5->SetLineWidth(5);
            g3_5->SetLineStyle(9);
            g3_5->Draw("CPE same");

            g3_6->SetLineColor(kViolet);
            g3_6->Draw("CPE same");

            legend3->Draw();
            c3->Update();

            // gSystem->mkdir("output_buffer/fouriers/trkbin_4/Fourier_fit_each_Mult_20_20_CM_Pair", kTRUE);
            TString filename3 = Form("output_buffer/fouriers/dijet_16/v3_for_pt_%.2f_to_%.2f.png",  ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
            c3->SaveAs(filename3);

            delete c3;


            
            

        }





            
    // }

        
        
    // }// outFile->Close();
    
    
}

int main() {
    processRootFile();
    return 0;
}



