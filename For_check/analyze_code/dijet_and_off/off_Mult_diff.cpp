
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
   
    // TFile* file = new TFile("newroot/off/off_data_02.root");
    // TFile* file2 = new TFile("newroot/off/off_MC_02.root");
    


        // Determine the number of wtrk and wppt2 values

    TString graph_title = "Differential j_{T}";

    TCanvas* c1 = new TCanvas("c1", graph_title, 2000, 1600);

        // Create a pad that spans the entire canvas
    TPad* title_pad = new TPad("title_pad", "", 0.0, 0.0, 1.0, 1.0);
    title_pad->SetFillStyle(4000);  // Make the pad transparent

    // Draw the pad on the canvas
    title_pad->Draw();
    title_pad->cd();

    // Create dummy histograms for the x and y titles
    TH1F* h1 = new TH1F("h1", "", 1, 0, 1);
    h1->SetTitle(graph_title);
    h1->SetStats(kFALSE);  // Hide the statistics box
    h1->GetXaxis()->SetTitle("j_{T}^{trig}(GeV/c)");
    h1->GetYaxis()->SetTitle("V_{n}(j_{T}^{trig}, j_{T}^{ref}) / #sqrt{V_{n}(j_{T}^{ref}, j_{T}^{ref})}");
    // Hide the x-axis values
    h1->GetXaxis()->SetLabelSize(0);

    // Hide the y-axis values
    h1->GetYaxis()->SetLabelSize(0);

        // Hide the axis lines
    h1->GetXaxis()->SetAxisColor(0);
    h1->GetYaxis()->SetAxisColor(0);

    // h1->GetXaxis()->SetTitleOffset(1.02); // Adjust the factor as needed
    // h1->GetYaxis()->SetTitleOffset(1.04); // Adjust the factor as needed
    h1->GetYaxis()->CenterTitle(true);
    h1->GetXaxis()->CenterTitle(true);
    
    // title_pad->SetTopMargin(0.15); // Increase this value if the title is cut off
    gStyle->SetTitleFontSize(0.04); // Adjust the size as needed
    gStyle->SetTitleOffset(1.1, "t"); // Adjust the offset as needed
    h1->SetFillColor(0);
    h1->GetYaxis()->SetTitleSize(0.027);
    h1->GetXaxis()->SetTitleSize(0.027);
    
    // h1->SetTitleOffset(0.95);
    h1->Draw();


   

    // // Create a 2D array of TPad objects
    TPad** pads = new TPad*[trackbin];

    double margin = 0.09;

    // Calculate the width and height of the pads
    double pad_width = (1.0 - 2*margin) / 3;
    double pad_height = (1.0 - 2*margin) / 3;

    for (int wtrk = 0; wtrk < trackbin; wtrk++) {
        // Calculate the position of the pad
        double x1 = margin + (wtrk % 3) * pad_width;
        double y1 = margin + (wtrk / 3) * pad_height;
        double x2 = x1 + pad_width;
        double y2 = y1 + pad_height;

        // Create the pad
        pads[wtrk] = new TPad(Form("pad_%d",wtrk), "", x1, y1, x2, y2);
        if (wtrk==1||wtrk==2){
            pads[wtrk]->SetMargin(0, 0, 0.1, 0);  // Set all margins to 0
        }
        else if (wtrk==3||wtrk==6) pads[wtrk]->SetMargin(0.1, 0, 0, 0);
        else if (wtrk==0) pads[wtrk]->SetMargin(0.1, 0, 0.1, 0);
        else pads[wtrk]->SetMargin(0, 0, 0, 0);
        pads[wtrk]->Draw();
    }

    c1->Update();


    

    std::string folders[] = {"data","MC"};

    for (const auto& folder : folders) {
        
        // TString filename = Form("newroot/off/off_diff_%s_02.root", folder.c_str());
        TString filename = Form("newroot/off/off_diff_%s_4.root", folder.c_str());
        TFile* file = new TFile(filename);
    
        for (int wtrk = 1 ; wtrk < trackbin+1; wtrk++ ){

            int wppt2 = ptbin;
            TGraphErrors* g2 = new TGraphErrors(ptbin-1);
            

            int point = 0;
            // int wtrk = 1;
            
            

            for (int wppt = 1 ; wppt < wppt2; wppt++) {
            
                TH2D* hSig = (TH2D*)file->Get(Form("hSigS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                TH2D* hBack = (TH2D*)file->Get(Form("hBckS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                // TH1D* hJt_high = (TH1D*)file->Get(Form("hJt_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                TH1D* hJt_low = (TH1D*)file->Get(Form("hJt_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                
                double xvalue = hJt_low->GetMean();
                double err_xvalue = hJt_low->GetMeanError();


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
                
                std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";
                TF1* func = new TF1("func", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                
                
                func->SetParameter(0, h1D->GetMaximum());
                func->SetParameter(1, -0.1);
                
                func->SetParameter(2, 0.0225);
                func->SetParameter(3, -0.01);
                
                h1D->Fit(func, "q 0");
                h1D->Fit(func, "q 0");
                h1D->Fit(func, "m q 0");
                h1D->Fit(func, "m q 0");
                h1D->Fit(func, "m q E 0");
                h1D->Fit(func, "m E q 0");

                // // std::cout<< Form("hSigS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])) + Form(" v2 = %.3f",func->GetParameter(2))<<std::endl;
                // std::cout<< std::string(Form("hSigS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1]))) + std::string(Form(" V2 = %.3f",func->GetParameter(2))) << std::endl;
                // hSigS_0_to_20_pt_3_to_7_pt_3_to_30
                // hSigS_0_to_20_pt_3_to_7_pt_3_to_30

                //////////////////////////////////////////////////////////////////////////////////




                TH2D* hSig2 = (TH2D*)file->Get(Form("hSigS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                TH2D* hBack2 = (TH2D*)file->Get(Form("hBckS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                
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
                
                TF1* func_2 = new TF1("func_2", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                

                func_2->SetParameter(0, h1D_2->GetMaximum());
                func_2->SetParameter(1, -0.1);
                
                func_2->SetParameter(2, 0.0225);
                func_2->SetParameter(3, -0.01);
                
                
                h1D_2->Fit(func_2, "q 0");
                h1D_2->Fit(func_2, "q 0");
                h1D_2->Fit(func_2, "m q 0");
                h1D_2->Fit(func_2, "m q 0");
                h1D_2->Fit(func_2, "m q E 0");
                h1D_2->Fit(func_2, "m E q 0");

                
                
            
                double r2 = func->GetParameter(2)/sqrt(func_2->GetParameter(2));
                double err2 = r2*sqrt(pow((sqrt(2)*func->GetParError(2)/func->GetParameter(2)),2) + 0.25*pow((sqrt(2)*func_2->GetParError(2)/func_2->GetParameter(2)),2) );

                                
                g2->SetPoint(point, xvalue, r2);
                g2->SetPointError(point, err_xvalue, err2);
                
                point++;

            
            }



            // std::cout<<"111"<<std::endl;
            
            double ymin = 0.0;
            double ymax = 1.0;
            
            g2->GetYaxis()->SetRangeUser(ymin, ymax);

        
            if(folder=="data"){
            
                g2->GetXaxis()->SetTitleSize(0.04); // Increase the size as needed
                g2->GetYaxis()->SetTitleSize(0.04); // Increase the size as needed
                g2->GetXaxis()->SetLabelSize(0.05); // Set the x-axis label size
                g2->GetYaxis()->SetLabelSize(0.05); // Set the y-axis label size
                // g2->GetXaxis()->SetTitleOffset(0.95); // Adjust the factor as needed
                // g2->GetYaxis()->SetTitleOffset(0.95); // Adjust the factor as needed
                g2->GetYaxis()->CenterTitle(true);
                g2->GetXaxis()->CenterTitle(true);

                // if(wtrk==1)


                pads[wtrk-1]->cd();

                g2->SetMarkerColor(kBlue);
                g2->SetLineColor(kBlue);
                g2->SetMarkerStyle(20);
                g2->SetLineWidth(2);
                g2->SetMarkerSize(3);
                g2->SetTitle(""); 
                g2->Draw("ACP");
                // g1_2->SetLineWidth(3);
                if (wtrk==trackbin)
                {   TLegend* legend = new TLegend (0.75,0.85,0.95,0.95);
                    legend->AddEntry(g2, "data", "p");
                    legend->SetBorderSize(0);
                    legend->Draw();
                }
                

            }
            else {

                pads[wtrk-1]->cd();
    


                g2->SetMarkerColor(kBlue);
                g2->SetLineColor(kBlue);
                g2->SetMarkerStyle(4);
                g2->SetLineStyle(9);
                g2->SetLineWidth(2);
                g2->SetMarkerSize(3);

                if (wtrk==trackbin)
                {   TLegend* legend_MC = new TLegend (0.75,0.75,0.95,0.85);
                    legend_MC->AddEntry(g2, "PYTHIA", "p");
                    legend_MC->SetBorderSize(0);
                    legend_MC->Draw("same");
                }
                g2->SetTitle(""); 
                g2->Draw("CP same");
            }

            

            
        
            // if (wtrk==1){
            //     TLatex latex;
            //     latex.SetTextSize(0.15);
            //     latex.DrawLatexNDC(0.4, 0.2, Form("j_{T}^{a}: %.2f to %.2f GeV/c", ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]));  // Adjust the coordinates and text as needed

            // }
            

            if (!(wtrk%3-1)){
                TLatex latex;
                latex.SetTextSize(0.07);
                latex.DrawLatexNDC(0.15, 0.8, Form("N_{ch}: %d to %d", trackbinbounds[wtrk-1], trackbinboundsUpper[wtrk-1]));  // Adjust the coordinates and text as needed

            }
            else{
                TLatex latex;
                latex.SetTextSize(0.07);
                latex.DrawLatexNDC(0.05, 0.8, Form("N_{ch}: %d to %d", trackbinbounds[wtrk-1], trackbinboundsUpper[wtrk-1]));  // Adjust the coordinates and text as needed
            }

            double x_min = g2->GetXaxis()->GetXmin();
            // double x_max = g1->GetXaxis()->GetXmax();
            double x_max = g2->GetX()[g2->GetN()-1];

            TLine* line1 = new TLine(x_min, 0.20, x_max, 0.20);
            line1->SetLineColor(kGreen-1);
            line1->SetLineStyle(9);
            line1->SetLineWidth(3);
            line1->Draw("same");


            TLine* line2 = new TLine(x_min, 0.40, x_max, 0.40);
            line2->SetLineColor(kGreen-1);
            line2->SetLineStyle(9);
            line2->SetLineWidth(3);
            line2->Draw("same");

            TLine* line3 = new TLine(x_min, 0.60, x_max, 0.60);
            line3->SetLineColor(kGreen-1);
            line3->SetLineStyle(9);
            line3->SetLineWidth(3);
            line3->Draw("same");

            // TLine* line4 = new TLine(x_min, 0.20, x_max, 0.20);
            // line4->SetLineColor(kRed);
            // line4->SetLineStyle(9);
            // line4->SetLineWidth(2);
            // line4->Draw("same");

            


            // legend->Draw();
            c1->Update();


                

                 

                
                

                
            // }

                        
        }

    }

    


    

    

    






    




    c1->Update();

    gSystem->mkdir(Form("output_buffer/off/"), kTRUE);

    TString filename = Form("output_buffer/off/off_diff_4.png");
    // std::cout << "Saving canvas to: " << filename << std::endl;
    c1->SaveAs(filename);
    delete c1;
    // }
    
    
}

int main() {
    processRootFile();
    return 0;
}



