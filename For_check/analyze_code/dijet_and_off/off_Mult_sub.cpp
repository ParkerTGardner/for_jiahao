
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
#include "../include/sherpa_constants_2_2.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "math.h"
#include "TGraphErrors.h" 
#include "TLatex.h"
#include "TLine.h"
#include "TMinuit.h"
void processRootFile() {









///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    //This file I need to mddify a little bit. r1 and r3 are invaild.
    // I want to overlay data and MC of r2 and adjust Markersize and textsize
    // but the main body should be this

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    
    // open files and adjust pad
    
    
    
    TFile* file = new TFile("newroot/off/off_data_02.root");

    TString graph_title = "Off Data";

    TCanvas* c1 = new TCanvas("c1", graph_title, 6400, 4800);

        // Create a pad that spans the entire canvas
    TPad* title_pad = new TPad("title_pad", "", 0, 0, 1, 1);
    title_pad->SetFillStyle(4000);  // Make the pad transparent

    // Draw the pad on the canvas
    title_pad->Draw();
    title_pad->cd();

    // Create dummy histograms for the x and y titles
    TH1F* h1 = new TH1F("h1", "", 1, 0, 1);
    h1->SetStats(kFALSE);  // Hide the statistics box
    h1->GetXaxis()->SetTitle("j_{T}^{a}-j_{T}^{b}(GeV/c)");
    h1->GetYaxis()->SetTitle("r_{n}(j_{T}^{a},j_{T}^{b})");
    // Hide the x-axis values
    h1->GetXaxis()->SetLabelSize(0);

    // Hide the y-axis values
    h1->GetYaxis()->SetLabelSize(0);

        // Hide the axis lines
    h1->GetXaxis()->SetAxisColor(0);
    h1->GetYaxis()->SetAxisColor(0);

    h1->GetXaxis()->SetTitleOffset(1.3); // Adjust the factor as needed
    h1->GetYaxis()->SetTitleOffset(1.2); // Adjust the factor as needed
    h1->GetYaxis()->CenterTitle(true);
    h1->GetXaxis()->CenterTitle(true);
    h1->SetTitle(graph_title);
    h1->SetTitleOffset(1.3);
    h1->SetFillColor(0);
    h1->GetYaxis()->SetTitleSize(0.04);
    h1->Draw();


   

    // Create a 2D array of TPad objects
    TPad*** pads = new TPad**[ptbin-1];
    for (int i = 0; i < ptbin-1; i++) {
        pads[i] = new TPad*[trackbin];
    }

    double margin = 0.05;

    // Calculate the width and height of the pads
    double pad_width = (1.0 - 2*margin) / (ptbin-1);
    double pad_height = (1.0 - 2*margin) / trackbin;

    for (int wppt2 = 0; wppt2 < ptbin-1; wppt2++) {

        for (int wtrk = 0; wtrk < trackbin; wtrk++) {

            // Calculate the position of the pad
            double y1 = margin + wtrk * pad_height;
            double x1 = margin + wppt2 * pad_width;
            double y2 = y1 + pad_height;
            double x2 = x1 + pad_width;

            // Create the pad
            
            pads[wppt2][wtrk] = new TPad(Form("pad_%d_%d", wppt2, wtrk), "", x1, y1, x2, y2);
            
            if ( wppt2!=0 && wtrk!=0 ) pads[wppt2][wtrk]->SetMargin(0, 0, 0, 0);
            
            if (wppt2==0 && wtrk!=0 ) pads[wppt2][wtrk]->SetMargin(0.05, 0, 0, 0);

            if (wtrk ==0 && wppt2!=0 ) pads[wppt2][wtrk]->SetMargin(0, 0, 0.05, 0);

            if (wtrk ==0 && wppt2 ==0 ) pads[wppt2][wtrk]->SetMargin(0.05, 0, 0.05, 0);
            // Draw the pad on the canvas
            pads[wppt2][wtrk]->Draw();


        
        }
            
    }

    c1->Update();



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    // main body

    // wppt2 corresponds jTb, and wppt for jta

    

    
    for (int wtrk = 1 ; wtrk < trackbin+1; wtrk++ ){

        for (int wppt2 = ptbin; wppt2 > 1; wppt2--) {


            TGraphErrors* g1 = new TGraphErrors(wppt2);
            TGraphErrors* g2 = new TGraphErrors(wppt2);
            TGraphErrors* g3 = new TGraphErrors(wppt2);
            TGraphErrors* g4 = new TGraphErrors(wppt2);

            int point = 0;


            for (int wppt = wppt2 ; wppt > 0; wppt--) {
            
                TH2D* hSig = (TH2D*)file->Get(Form("hSigS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                TH2D* hBack = (TH2D*)file->Get(Form("hBckS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                TH1D* hJt_high = (TH1D*)file->Get(Form("hJt_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                TH1D* hJt_low = (TH1D*)file->Get(Form("hJt_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                
                double xvalue = hJt_high->GetMean() - hJt_low->GetMean();
                double err_high = hJt_high->GetMeanError();
                double err_low = hJt_low->GetMeanError();
                double err_xvalue = sqrt(pow(err_high, 2) + pow(err_low, 2));



                hSig->SetDefaultSumw2(kTRUE);
                hBack->SetDefaultSumw2(kTRUE);

                double DelEta_lo = 2.0;
                // double DelEta_mid = 4.7;
                double DelEta_hi  =  6.0;
                TH1D *hSig1D = (TH1D*) hSig->ProjectionY("h1D", hSig->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
                TH1D *hBack1D = (TH1D*) hBack->ProjectionY("h1D", hBack->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
                hSig1D->SetDefaultSumw2(kTRUE);
                hBack1D->SetDefaultSumw2(kTRUE);
                hSig1D->Divide(hBack1D);
                hSig1D->Scale(hBack1D->GetMaximum());
                TH1D* h1D = hSig1D;
                h1D->SetDefaultSumw2(kTRUE);
                
                std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)))";
                TF1* func = new TF1("func", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                

                
                func->SetParameter(0, h1D->GetMaximum());
                func->SetParameter(1, 0.05);
                
                func->SetParameter(2, 0.05);
                func->SetParameter(3, 0.01);
                func->SetParameter(4, 0.01);
                h1D->Fit(func, "q 0");
                h1D->Fit(func, "q 0");
                h1D->Fit(func, "m q 0");
                h1D->Fit(func, "m q 0");
                h1D->Fit(func, "m q E 0");
                h1D->Fit(func, "m E q 0");



                //////////////////////////////////////////////////////////////////////////////////




                TH2D* hSig2 = (TH2D*)file->Get(Form("hSigS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                TH2D* hBack2 = (TH2D*)file->Get(Form("hBckS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                
                hSig2->SetDefaultSumw2(kTRUE);
                hBack2->SetDefaultSumw2(kTRUE);

                

            
                TH1D *hSig1D_2 = (TH1D*) hSig2->ProjectionY("h1D", hSig2->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
                TH1D *hBack1D_2 = (TH1D*) hBack2->ProjectionY("h1D", hBack2->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
                hSig1D_2->SetDefaultSumw2(kTRUE);
                hBack1D_2->SetDefaultSumw2(kTRUE);
                hSig1D_2->Divide(hBack1D_2);
                hSig1D_2->Scale(hBack1D_2->GetMaximum());
                TH1D* h1D_2 = hSig1D_2;
                h1D_2->SetDefaultSumw2(kTRUE);
                
                // std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)))";
                TF1* func_2 = new TF1("func_2", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                


                func_2->SetParameter(0, h1D_2->GetMaximum());
                func_2->SetParameter(1, 0.05);
                func_2->SetParameter(2, 0.05);
                func_2->SetParameter(3, 0.01);
                func_2->SetParameter(4, 0.01);
                
                
                h1D_2->Fit(func_2, "q 0");
                h1D_2->Fit(func_2, "q 0");
                h1D_2->Fit(func_2, "m q 0");
                h1D_2->Fit(func_2, "m q 0");
                h1D_2->Fit(func_2, "m q E 0");
                h1D_2->Fit(func_2, "m E q 0");
                



                /////////////////////////////////////////////////



                TH2D* hSig3 = (TH2D*)file->Get(Form("hSigS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                TH2D* hBack3 = (TH2D*)file->Get(Form("hBckS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                
                hSig3->SetDefaultSumw2(kTRUE);
                hBack3->SetDefaultSumw2(kTRUE);


                
                TH1D *hSig1D_3 = (TH1D*) hSig3->ProjectionY("h1D", hSig3->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
                TH1D *hBack1D_3 = (TH1D*) hBack3->ProjectionY("h1D", hBack3->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
                hSig1D_3->SetDefaultSumw2(kTRUE);
                hBack1D_3->SetDefaultSumw2(kTRUE);
                hSig1D_3->Divide(hBack1D_3);
                hSig1D_3->Scale(hBack1D_3->GetMaximum());
                TH1D* h1D_3 = hSig1D_3;
                h1D_3->SetDefaultSumw2(kTRUE);
                TF1* func_3 = new TF1("func_3", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                

                func_3->SetParameter(0, h1D_3->GetMaximum());
                func_3->SetParameter(1, 0.05);
                func_3->SetParameter(2, 0.05);
                func_3->SetParameter(3, 0.01);
                func_3->SetParameter(4, 0.01);
                
                
                h1D_3->Fit(func_3, "q 0");
                h1D_3->Fit(func_3, "q 0");
                h1D_3->Fit(func_3, "m q 0");
                h1D_3->Fit(func_3, "m q 0");
                h1D_3->Fit(func_3, "m q E 0");
                h1D_3->Fit(func_3, "m E q 0");
            
                

                double r1 = fabs(func->GetParameter(1)/sqrt(func_2->GetParameter(1)*func_3->GetParameter(1)));
                double r2 = fabs(func->GetParameter(2)/sqrt(func_2->GetParameter(2)*func_3->GetParameter(2)));
                double r3 = fabs(func->GetParameter(3)/sqrt(func_2->GetParameter(3)*func_3->GetParameter(3)));
                double err1 = r1*sqrt(pow((func->GetParError(1)/func->GetParameter(1)),2) + 0.25*pow((func_2->GetParError(1)/func->GetParameter(1)),2) + 0.25*pow((func_3->GetParError(1)/func_3->GetParameter(1)),2));
                double err2 = r2*sqrt(pow((func->GetParError(2)/func->GetParameter(2)),2) + 0.25*pow((func_2->GetParError(2)/func->GetParameter(2)),2) + 0.25*pow((func_3->GetParError(2)/func_3->GetParameter(2)),2));
                double err3 = r3*sqrt(pow((func->GetParError(3)/func->GetParameter(3)),2) + 0.25*pow((func_2->GetParError(3)/func->GetParameter(3)),2) + 0.25*pow((func_3->GetParError(3)/func_3->GetParameter(3)),2));


                g1->SetPoint(point, xvalue, r1);
                g1->SetPointError(point, err_xvalue, err1);
                                
                g2->SetPoint(point, xvalue, r2);
                g2->SetPointError(point, err_xvalue, err2);

                g3->SetPoint(point, xvalue, r3);
                g3->SetPointError(point, err_xvalue, err3);

                
                
                point++;


            
            
            }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            

            
            
            double ymin = 0.0;
            double ymax = 1.1;
            
            g1->GetYaxis()->SetRangeUser(ymin, ymax);
            g2->GetYaxis()->SetRangeUser(ymin, ymax);
            g3->GetYaxis()->SetRangeUser(ymin, ymax);
            g4->GetYaxis()->SetRangeUser(ymin, ymax);
          

            
            g1->GetXaxis()->SetTitleSize(0.04); // Increase the size as needed
            g1->GetYaxis()->SetTitleSize(0.04); // Increase the size as needed
            // 
            // g1->GetYaxis()->SetTitle("r_{n}(j_{T}^{a}-j_{T}^{b})");
            g1->GetXaxis()->SetTitleOffset(0.95); // Adjust the factor as needed
            g1->GetYaxis()->SetTitleOffset(0.95); // Adjust the factor as needed
            g1->GetYaxis()->CenterTitle(true);
            g1->GetXaxis()->CenterTitle(true);

            // if(wtrk==1)


            pads[wppt2-2][wtrk-1]->cd();

            g1->SetMarkerColor(kRed);
            g1->SetLineColor(kRed);
            g1->SetMarkerStyle(20);
            // g1->SetLineStyle(2);
            // g1_2->SetLineWidth(3);

            g1->SetTitle(""); 
            g1->Draw("ACP");

            g2->SetMarkerStyle(21);
            g2->SetMarkerColor(kBlue);
            g2->SetLineColor(kBlue);
            // g2->SetLineStyle(2);
            // g2_2->SetLineWidth(3); 
            g2->Draw("CP same");

            if(wtrk < trackbin) {
                g3->SetMarkerStyle(22);
                g3->SetMarkerColor(kGreen);
                g3->SetLineColor(kGreen);
                g3->Draw("CP same");
            }
            

            if (wtrk==1){
                TLatex latex;
                latex.SetTextSize(0.1);
                latex.DrawLatexNDC(0.35, 0.15, Form("j_{T}^{a}: %.2f to %.2f GeV/c", ptbinbounds_lo[wppt2-1], ptbinbounds_hi[wppt2-1]));  // Adjust the coordinates and text as needed

            }

            if (wppt2==2){
                TLatex latex;
                latex.SetTextSize(0.1);
                latex.DrawLatexNDC(0.1, 0.4, Form("N_{ch}: %d to %d", trackbinbounds[wtrk-1], trackbinboundsUpper[wtrk-1]));  // Adjust the coordinates and text as needed

            }

            double x_min = g1->GetXaxis()->GetXmin();
            // double x_max = g1->GetXaxis()->GetXmax();
            double x_max = g1->GetX()[g1->GetN()-1];

            TLine* line1 = new TLine(x_min, 1, x_max, 1);
            line1->SetLineColor(kBlack);
            line1->SetLineStyle(9);
            line1->Draw();

            TLine* line2 = new TLine(x_min, 0.6, x_max, 0.6);
            line2->SetLineColor(kBlue);
            line2->SetLineStyle(9);
            line2->Draw();

            TLine* line3 = new TLine(x_min, 0.4, x_max, 0.4);
            line3->SetLineColor(kGreen);
            line3->SetLineStyle(9);
            line3->Draw();

            // TLine* line4 = new TLine(x_min, 0.9, x_max, 0.9);
            // line4->SetLineColor(kBlue);
            // line4->SetLineStyle(9);
            // line4->Draw();

            


            // legend->Draw();
            c1->Update();

            
            

            
        }

                       
    }

    


    

    

    






    




    c1->Update();

    gSystem->mkdir(Form("output_buffer/off/"), kTRUE);

    TString filename = Form("output_buffer/off/off_data_02.png");
    // std::cout << "Saving canvas to: " << filename << std::endl;
    c1->SaveAs(filename);
    delete c1;
    // }
    
    
}

int main() {
    processRootFile();
    return 0;
}



