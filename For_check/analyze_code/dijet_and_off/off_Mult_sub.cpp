
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
#include "../include/sherpa_constants_2_3.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "math.h"
#include "TGraphErrors.h" 
#include "TLatex.h"
#include "TLine.h"
#include "TMinuit.h"
#include "../include/fourier_fit.h"
void processRootFile() {
   
    
        // Determine the number of wtrk and wppt2 values

    TString graph_title = "Factorization Ratio ";

    TCanvas* c1 = new TCanvas("c1", graph_title, 3300, 1600);

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
    h1->GetYaxis()->SetTitle("r_{2}(j_{T}^{a},j_{T}^{b})");
    // Hide the x-axis values
    h1->GetXaxis()->SetLabelSize(0);

    // Hide the y-axis values
    h1->GetYaxis()->SetLabelSize(0);

        // Hide the axis lines
    h1->GetXaxis()->SetAxisColor(0);
    h1->GetYaxis()->SetAxisColor(0);

    h1->GetXaxis()->SetTitleOffset(1.3); // Adjust the factor as needed
    h1->GetYaxis()->SetTitleOffset(1); // Adjust the factor as needed
    h1->GetYaxis()->CenterTitle(true);
    h1->GetXaxis()->CenterTitle(true);
    h1->SetTitle(graph_title);
    gStyle->SetTitleFontSize(0.06); // Adjust the size as needed
    gStyle->SetTitleOffset(1.2, "t"); // Adjust the offset as needed
    h1->SetFillColor(0);
    h1->GetYaxis()->SetTitleSize(0.04);
    // h1->SetTitleSize(0.05);
    // h1->SetTitleOffset(0.95);
    h1->Draw();







    // Create a 2D array of TPad objects
int trackbin_plot = 3;
int wtrk_values[] = {3, 6, 8};
TPad*** pads = new TPad**[trackbin_plot];
for (int i = 0; i < trackbin_plot; i++) {
    pads[i] = new TPad*[ptbin-2];
}

double margin = 0.07;

// Calculate the width and height of the pads
double pad_width = (1.0 - 2*margin) / (ptbin-2);
double pad_height = (1.0 - 2*margin) / trackbin_plot;

for (int wtrk = 0; wtrk < trackbin_plot; wtrk++) {
    for (int wppt2 = 0; wppt2 < ptbin-2; wppt2++) {
    
        // Calculate the position of the pad
        double x1 = margin + wppt2 * pad_width;
        double y1 = margin + wtrk * pad_height;
        double x2 = x1 + pad_width;
        double y2 = y1 + pad_height;

        // Create the pad
        
        pads[wtrk][wppt2] = new TPad(Form("pad_%d_%d", wtrk, wppt2), "", x1, y1, x2, y2);
        
        if ( wtrk!=0 && wppt2!=0 ) pads[wtrk][wppt2]->SetMargin(0, 0, 0, 0);
        if (wtrk==0 && wppt2!=0 ) pads[wtrk][wppt2]->SetMargin(0, 0, 0.1, 0);
        if (wppt2 ==0 && wtrk!=0 ) pads[wtrk][wppt2]->SetMargin(0.1, 0, 0, 0);
        if (wppt2 ==0 && wtrk ==0 ) pads[wtrk][wppt2]->SetMargin(0.1, 0, 0.1, 0);
         
        // Draw the pad on the canvas
        pads[wtrk][wppt2]->Draw();
    }
}


    

    // Update the canvas
    c1->Update();



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // This is the main code



    // The differential jT expression is V(a, b)/sqrt(V(a,a))*sqrt(V(b,b))
    // So I read 2D from (a,b), (a,a) and (b,b)
    // Here wppt:= a and wppt2:=b 
    // The x := the mean of jT(b)-jT(a)
    // By convention, we usually consider a>b, so the x lable is written as a-b

    std::string folders[] = {"data","MC"};

    for (const auto& folder : folders) {
        
        TString filename = Form("newroot/off/off_%s_02.root", folder.c_str());
        TFile* file = new TFile(filename);

    // trackbin loop
        for (int i = 0 ; i < trackbin_plot; i++ ){

            int wtrk = wtrk_values[i];
    //  for each trackbin, for each b here is a b-a curve
            for (int wppt2 = ptbin; wppt2 > 2; wppt2--) {

                TGraphErrors* g2 = new TGraphErrors(wppt2);
                

                int point = 0;


    // loop over a
                for (int wppt = wppt2 ; wppt > 0; wppt--) {

                    // x value
                    TH1D* hJt_high = (TH1D*)file->Get(Form("hJt_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                    TH1D* hJt_low = (TH1D*)file->Get(Form("hJt_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                    
                    double xvalue = hJt_high->GetMean() - hJt_low->GetMean();
                    double err_high = hJt_high->GetMeanError();
                    double err_low = hJt_low->GetMeanError();
                    double err_xvalue = sqrt(pow(err_high, 2) + pow(err_low, 2));


                    //V(a,b)
                
                    TH2D* hSig = (TH2D*)file->Get(Form("hSigS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                    TH2D* hBack = (TH2D*)file->Get(Form("hBckS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                    

                    hSig->SetDefaultSumw2(kTRUE);
                    hBack->SetDefaultSumw2(kTRUE);

                    TF1* func = fourier_fit(hSig,hBack); 



                    //////////////////////////////////////////////////////////////////////////////////

                    // V(a,a)


                    TH2D* hSig2 = (TH2D*)file->Get(Form("hSigS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1]), (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                    TH2D* hBack2 = (TH2D*)file->Get(Form("hBckS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1]), (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                    
                    hSig2->SetDefaultSumw2(kTRUE);
                    hBack2->SetDefaultSumw2(kTRUE);

                    TF1* func_2 = fourier_fit(hSig2,hBack2); 

                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    //V(b,b)

                    TH2D* hSig3 = (TH2D*)file->Get(Form("hSigS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                    TH2D* hBack3 = (TH2D*)file->Get(Form("hBckS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1]), (int)(10*ptbinbounds_lo[wppt2-1]), (int)(10*ptbinbounds_hi[wppt2-1])));
                    
                    hSig3->SetDefaultSumw2(kTRUE);
                    hBack3->SetDefaultSumw2(kTRUE);

                    TF1* func_3 = fourier_fit(hSig3,hBack3);
                    
                    
                
                    double r2 = func->GetParameter(2)/sqrt(func_2->GetParameter(2)*func_3->GetParameter(2));
                    double err2 = r2*sqrt(pow((sqrt(2)*func->GetParError(2)/func->GetParameter(2)),2) + 0.25*pow((sqrt(2)*func_2->GetParError(2)/func_2->GetParameter(2)),2) + 0.25*pow((sqrt(2)*func_3->GetParError(2)/func_3->GetParameter(2)),2));

                    if (wppt!=wppt2){
                        g2->SetPoint(point, xvalue, r2);
                        g2->SetPointError(point, err_xvalue, err2);
                    }
                    else {
                        g2->SetPoint(point, xvalue, r2);
                        g2->SetPointError(point, err_xvalue, 0);
                    }              
                    
                    
                    point++;

                
                }

                // end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                
                double ymin = 0.6;
                double ymax = 1.4;
                
                g2->GetYaxis()->SetRangeUser(ymin, ymax);

            
                if(folder=="data"){
                
                    g2->GetXaxis()->SetTitleSize(0.04); // Increase the size as needed
                    g2->GetYaxis()->SetTitleSize(0.04); // Increase the size as needed
                    g2->GetXaxis()->SetLabelSize(0.08); // Set the x-axis label size
                    g2->GetYaxis()->SetLabelSize(0.08); // Set the y-axis label size
                    g2->GetXaxis()->SetTitleOffset(1.1); // Adjust the factor as needed
                    g2->GetYaxis()->SetTitleOffset(1.1); // Adjust the factor as needed
                    g2->GetYaxis()->CenterTitle(true);
                    g2->GetXaxis()->CenterTitle(true);

                    


                    pads[i][wppt2-3]->cd();
                    
                    g2->SetMarkerColor(kBlue);
                    g2->SetLineColor(kBlue);
                    g2->SetMarkerStyle(20);
                    g2->SetLineWidth(2);
                    g2->SetMarkerSize(4);
                    // g1_2->SetLineWidth(3);

                    g2->SetTitle(""); 
                    g2->Draw("ACP");
                    if (wtrk==8 && wppt2 == ptbin)
                    {   TLegend* legend = new TLegend (0.7,0.90,0.9,0.98);
                        legend->AddEntry(g2, "data", "p");
                        legend->SetBorderSize(0);
                        legend->Draw("same");
                    }
                    

                }
                else {
                    std::cout<<i<<std::endl;
std::cout<<wppt2-3<<std::endl;

                    pads[i][wppt2-3]->cd();

                    g2->SetMarkerColor(kBlue);
                    g2->SetLineColor(kBlue);
                    g2->SetMarkerStyle(4);
                    g2->SetLineStyle(9);
                    g2->SetLineWidth(2);
                    g2->SetMarkerSize(4);
                    // g1_2->SetLineWidth(3);

                    g2->SetTitle(""); 
                    g2->Draw("CP same");
                    if (wtrk==8 && wppt2==ptbin)
                    {   TLegend* legend_MC = new TLegend (0.7,0.80,0.9,0.90);
                        legend_MC->AddEntry(g2, "PYTHIA", "p");
                        legend_MC->SetBorderSize(0);
                        legend_MC->Draw("same");
                    }
                }

                

                
                if (i==0){
                    if(wppt2==3){
                        TLatex latex;
                        latex.SetTextSize(0.08);
                        latex.DrawLatexNDC(0.2, 0.4, Form("N_{ch}: %d to %d", trackbinbounds[wtrk-1], trackbinboundsUpper[wtrk-1]));  // Adjust the coordinates and text as needed

                        TLatex latex2;
                        latex2.SetTextSize(0.08);
                        latex2.DrawLatexNDC(0.4, 0.2, Form("j_{T}^{a}: %.2f to %.2f GeV/c", ptbinbounds_lo[wppt2-1], ptbinbounds_hi[wppt2-1]));  // Adjust the coordinates and text as needed
                    }
                    else{
                        TLatex latex;
                        latex.SetTextSize(0.08);
                        latex.DrawLatexNDC(0.3, 0.2, Form("j_{T}^{a}: %.2f to %.2f GeV/c", ptbinbounds_lo[wppt2-1], ptbinbounds_hi[wppt2-1]));  // Adjust the coordinates and text as needed
                    }
                    

                } else if (wppt2==3){
                    TLatex latex;
                    latex.SetTextSize(0.08);
                    latex.DrawLatexNDC(0.2, 0.4, Form("N_{ch}: %d to %d", trackbinbounds[wtrk-1], trackbinboundsUpper[wtrk-1]));  // Adjust the coordinates and text as needed
                    
                }

                double x_min = g2->GetXaxis()->GetXmin();
                // double x_max = g1->GetXaxis()->GetXmax();
                double x_max = g2->GetX()[g2->GetN()-1];

                TLine* line1 = new TLine(x_min, 1, x_max, 1);
                line1->SetLineColor(kBlack);
                line1->SetLineStyle(9);
                line1->SetLineWidth(5);
                line1->Draw("same");


                TLine* line2 = new TLine(x_min, 0.8, x_max, 0.8);
                line2->SetLineColor(kGreen-1);
                line2->SetLineStyle(9);
                line2->SetLineWidth(4);
                line2->Draw("same");

                TLine* line3 = new TLine(x_min, 1.2, x_max, 1.2);
                line3->SetLineColor(kGreen-1);
                line3->SetLineStyle(9);
                line3->SetLineWidth(4);
                line3->Draw("same");

                // TLine* line4 = new TLine(x_min, 0.8, x_max, 0.8);
                // line4->SetLineColor(kGreen-1);
                // line4->SetLineStyle(9);
                // line4->SetLineWidth(3);
                // line4->Draw("same");

                


                // legend->Draw();
                c1->Update();

                
                

                
            }

                        
        }

    }

    


    

    

    






    




    c1->Update();

    gSystem->mkdir(Form("output_buffer/off/"), kTRUE);

    TString filename = Form("output_buffer/off/off_data_4_00.png");
    // std::cout << "Saving canvas to: " << filename << std::endl;
    c1->SaveAs(filename);
    delete c1;
    // }
    
    
}

int main() {
    processRootFile();
    return 0;
}



