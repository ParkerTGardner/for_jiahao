
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

    TString graph_title = "Nomalized j_{T} Spectrum";

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
    // h1->GetYaxis()->SetTitle("V_{n}(j_{T}^{trig}, j_{T}^{ref}) / #sqrt{V_{n}(j_{T}^{ref}, j_{T}^{ref})}");
    // Hide the x-axis values
    h1->GetXaxis()->SetLabelSize(0);

    // Hide the y-axis values
    h1->GetYaxis()->SetLabelSize(0);

        // Hide the axis lines
    h1->GetXaxis()->SetAxisColor(0);
    h1->GetYaxis()->SetAxisColor(0);

    h1->GetXaxis()->SetTitleOffset(1.15); // Adjust the factor as needed
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

    double margin = 0.07;

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

                
            TH1D* hJt_0_3 = (TH1D*)file->Get(Form("hJt_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],0, 3));
            TH1D* hJt_3_30 = (TH1D*)file->Get(Form("hJt_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], 3, 30));

            TH1D* hJt_0_30 = (TH1D*)hJt_3_30->Clone("hCombined");
            hJt_0_30->Add(hJt_0_3);
            hJt_0_30->SetStats(false);
            hJt_0_30->Scale(1.0/hJt_0_30->GetMaximum());
            //////////////////////////////////////////////////////////////////////////////////


        
            if(folder=="data"){
            
                hJt_0_30->GetXaxis()->SetTitleSize(0.04); // Increase the size as needed
                hJt_0_30->GetYaxis()->SetTitleSize(0.04); // Increase the size as needed
                hJt_0_30->GetXaxis()->SetLabelSize(0.08); // Set the x-axis label size
                hJt_0_30->GetYaxis()->SetLabelSize(0.08); // Set the y-axis label size
                // g2->GetXaxis()->SetTitleOffset(0.95); // Adjust the factor as needed
                // g2->GetYaxis()->SetTitleOffset(0.95); // Adjust the factor as needed
                hJt_0_30->GetYaxis()->CenterTitle(true);
                hJt_0_30->GetXaxis()->CenterTitle(true);



                pads[wtrk-1]->cd();

                hJt_0_30->SetMarkerColor(kBlue);
                hJt_0_30->SetLineColor(kBlue);
                hJt_0_30->SetMarkerStyle(20);
                hJt_0_30->SetLineWidth(2);
                hJt_0_30->SetMarkerSize(2);
                hJt_0_30->SetTitle(""); 
                hJt_0_30->Draw("CP");




                // g1_2->SetLineWidth(3);
                if (wtrk==trackbin)
                {   TLegend* legend = new TLegend (0.75,0.85,0.95,0.95);
                    legend->AddEntry(hJt_0_30, "data", "p");
                    legend->SetBorderSize(0);
                    legend->Draw("same");
                }
                

            }
            else {

                pads[wtrk-1]->cd();

                hJt_0_30->SetMarkerColor(kBlue);
                hJt_0_30->SetLineColor(kBlue);
                hJt_0_30->SetMarkerStyle(4);
                hJt_0_30->SetLineStyle(9);
                hJt_0_30->SetLineWidth(2);
                hJt_0_30->SetMarkerSize(2);
                hJt_0_30->SetTitle(""); 
                hJt_0_30->Draw("ACP same");


                if (wtrk==trackbin)
                {   TLegend* legend_MC = new TLegend (0.75,0.75,0.95,0.85);
                    legend_MC->AddEntry(hJt_0_30, "PYTHIA", "p");
                    legend_MC->SetBorderSize(0);
                    legend_MC->Draw("same");
                }
                
            }

            

            
        
            // if (wtrk==1){
            //     TLatex latex;
            //     latex.SetTextSize(0.15);
            //     latex.DrawLatexNDC(0.4, 0.2, Form("j_{T}^{a}: %.2f to %.2f GeV/c", ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]));  // Adjust the coordinates and text as needed

            // }
            

            if (!(wtrk%3-1)){
                TLatex latex;
                latex.SetTextSize(0.07);
                latex.DrawLatexNDC(0.5, 0.8, Form("N_{ch}: %d to %d", trackbinbounds[wtrk-1], trackbinboundsUpper[wtrk-1]));  // Adjust the coordinates and text as needed

            }
            else{
                TLatex latex;
                latex.SetTextSize(0.07);
                latex.DrawLatexNDC(0.4, 0.8, Form("N_{ch}: %d to %d", trackbinbounds[wtrk-1], trackbinboundsUpper[wtrk-1]));  // Adjust the coordinates and text as needed
            }

            // double x_min = g2->GetXaxis()->GetXmin();
            // // double x_max = g1->GetXaxis()->GetXmax();
            // double x_max = g2->GetX()[g2->GetN()-1];

            // TLine* line1 = new TLine(x_min, 0.20, x_max, 0.20);
            // line1->SetLineColor(kGreen-1);
            // line1->SetLineStyle(9);
            // line1->SetLineWidth(3);
            // line1->Draw("same");


            // TLine* line2 = new TLine(x_min, 0.40, x_max, 0.40);
            // line2->SetLineColor(kGreen-1);
            // line2->SetLineStyle(9);
            // line2->SetLineWidth(3);
            // line2->Draw("same");

            // TLine* line3 = new TLine(x_min, 0.60, x_max, 0.60);
            // line3->SetLineColor(kGreen-1);
            // line3->SetLineStyle(9);
            // line3->SetLineWidth(3);
            // line3->Draw("same");

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

    TString filename = Form("output_buffer/off/off_diff_jt_4.png");
    // std::cout << "Saving canvas to: " << filename << std::endl;
    c1->SaveAs(filename);
    delete c1;
    // }
    
    
}

int main() {
    processRootFile();
    return 0;
}


