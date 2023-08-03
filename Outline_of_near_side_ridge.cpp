
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
#include "sherpa_constants.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "math.h"
#include "TGraphErrors.h" 
#include "TGraphSmooth.h"
void processRootFile() {

    double DelEta_lo = 2.0;
    double phi_lo    = -M_PI/16;
    double phi_hi    = M_PI/16;
    
    for (int wppt = 1; wppt < ptbin+1; wppt++) {
      
        for (int wtrk = 1; wtrk < trackbin+1; wtrk++) {
            TString Outline = Form("Outline_for_track_%d_to_%d_pt_%.2f_to_%.2f.png", trackbinbounds[wtrk-1], trackbinboundsUpper[wtrk-1], ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
            // Variables to hold the v1, v2, v3 values
            TCanvas* c1 = new TCanvas("c1", Outline, 800, 600);
            TFile* file = new TFile("jobs_more_all.root");
            // // // Get the histogram
            TH2D* hSig = (TH2D*)file->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
            TH2D* hBack = (TH2D*)file->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
            hSig->SetDefaultSumw2(kTRUE);
            hBack->SetDefaultSumw2(kTRUE);

        
            TH1D *hSig1D = (TH1D*) hSig->ProjectionX("h1D", phi_lo, phi_hi)->Clone();
            TH1D *hBack1D = (TH1D*) hBack->ProjectionX("h1D", phi_lo, phi_hi)->Clone();
            hSig1D->SetDefaultSumw2(kTRUE);
            hBack1D->SetDefaultSumw2(kTRUE);
            hSig1D->Divide(hBack1D);
            // hSig1D->Scale(hBack1D->GetMaximum());
            // // hBack1D->Scale(hBack1D->GetMaximum());
            TH1D* h1D = hSig1D;
            h1D->SetDefaultSumw2(kTRUE);
            std::cout<< h1D->GetMinimum()<< std::endl;
            // h1D->Draw("hist CP");
                // std::cout << "draw" << std::endl;
            
        
            TFile* file2 = new TFile("jobs_more_all_2.root");
            
            TH2D* hSig2 = (TH2D*)file2->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
            TH2D* hBack2 = (TH2D*)file2->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
            hSig2->SetDefaultSumw2(kTRUE);
            hBack2->SetDefaultSumw2(kTRUE);

            // // double DelEta_lo = 2.0;
            TH1D *hSig1D_2 = (TH1D*) hSig2->ProjectionX("h1D", phi_lo, phi_hi)->Clone();
            TH1D *hBack1D_2 = (TH1D*) hBack2->ProjectionX("h1D", phi_lo, phi_hi)->Clone();
            hSig1D_2->SetDefaultSumw2(kTRUE);
            hBack1D_2->SetDefaultSumw2(kTRUE);
            hSig1D_2->Divide(hBack1D_2);
            hSig1D_2->Scale(hBack1D_2->GetMaximum());
            // // hBack1D->Scale(hBack1D_2->GetMaximum());
            TH1D* h1D_2 = hSig1D_2;
            h1D_2->SetDefaultSumw2(kTRUE);
            
            
            

            
            TFile* file3 = new TFile("jobs_more_all_3.root");
            TH2D* hSig3 = (TH2D*)file3->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
            TH2D* hBack3 = (TH2D*)file3->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
            hSig3->SetDefaultSumw2(kTRUE);
            hBack3->SetDefaultSumw2(kTRUE);

            // double DelEta_lo = 2.0;
            TH1D *hSig1D_3 = (TH1D*) hSig3->ProjectionX("h1D", phi_lo, phi_lo)->Clone();
            TH1D *hBack1D_3 = (TH1D*) hBack3->ProjectionX("h1D", phi_lo, phi_hi)->Clone();
            
            hSig1D_3->SetDefaultSumw2(kTRUE);
            hBack1D_3->SetDefaultSumw2(kTRUE);
            hSig1D_3->Divide(hBack1D_3);
            hSig1D_3->Scale(hBack1D_3->GetMaximum());
            // hBack1D->Scale(hBack1D_2->GetMaximum());
            TH1D* h1D_3 = hSig1D_3;
            h1D_3->SetDefaultSumw2(kTRUE);

            
            h1D->Divide(h1D_3);
            h1D_2->Divide(h1D_3);
            
            
            // file3->Close();




            // std::cout<< "legend" << std::endl;
            
           
            
            // std::cout<< h1D->GetMinimum() << std::endl;

            double ymin = std::min({h1D->GetMinimum(), h1D_2->GetMinimum()});
            double ymax = std::max({h1D->GetMaximum(), h1D_2->GetMaximum()});
            
            h1D->Scale(1/h1D->GetMaximum());
            h1D_2->Scale(1/h1D_2->GetMaximum());
            h1D_3->Scale(1/h1D_3->GetMaximum());
            TGraph* graph1 = new TGraph(h1D);
            TGraph* graph2 = new TGraph(h1D_2);
            TGraph* graph3 = new TGraph(h1D_3);
            TGraphSmooth* gs1 = new TGraphSmooth("normal");
            TGraphSmooth* gs2 = new TGraphSmooth("normal");
            TGraphSmooth* gs3 = new TGraphSmooth("normal");
            TGraph* g1 = gs1->SmoothKern(graph1, "normal", 0.1);
            TGraph* g2 = gs2->SmoothKern(graph2, "normal", 0.1);
            TGraph* g3 = gs3->SmoothKern(graph3, "normal", 0.1);

            // // Add a small buffer to the y range

            // std::cout<< "legend0.5" << std::endl;

            // double buffer = 0.3 * (ymax - ymin);
            // ymin -= buffer;
            // ymax += buffer;

            // // Set the y range for each graph
            // std::cout<< "legend1" << std::endl;

            // h1D->GetYaxis()->SetRangeUser(0, 1.2);
            // h1D_2->GetYaxis()->SetRangeUser(0, 1.2);
            // h1D_3->GetYaxis()->SetRangeUser(0, 1.2);


            // h1D->SetTitle(Outline);
            // h1D->GetXaxis()->SetTitle("#Delta#eta");

            // h1D->SetLineColor(kRed);
            // h1D->SetLineStyle(1);
            // h1D->Draw("hist");


            // h1D_2->SetLineColor(kBlue);
            // h1D_2->SetLineStyle(1);
            // h1D_2->Draw("hist same");


            // h1D_3->SetLineStyle(1);
            // h1D_3->SetLineColor(kBlack);
            // h1D_3->Draw("hist same");


            TLegend* legend = new TLegend(0.75,0.75,0.9,0.9); // Adjust these parameters to move the legend to the desired position
            legend->AddEntry(g1, "jt", "l");
            legend->AddEntry(g2, "jt_nor", "l");
            legend->AddEntry(g3, "unweighted", "l");


            g1->GetYaxis()->SetRangeUser(0, 1.2);
            g2->GetYaxis()->SetRangeUser(0, 1.2);
            g3->GetYaxis()->SetRangeUser(0, 1.2);

            // std::cout<< "legend2" << std::endl;

            g1->SetTitle(Outline);
            g1->GetXaxis()->SetTitle("#Delta#eta");

            g1->SetLineColor(kRed);
            g1->SetLineStyle(1);
            g1->Draw("AC");


            g2->SetLineColor(kBlue);
            g2->SetLineStyle(1);
            g2->Draw("C same");


            g3->SetLineStyle(1);
            g3->SetLineColor(kBlack);
            g3->Draw("C same");


            legend->Draw();
            c1->Update();

            gSystem->mkdir("Outline_of_nearside", kTRUE);

            TString filename = "Outline_of_nearside/"+Outline;
            std::cout << "Saving canvas to: " << filename << std::endl;
            c1->SaveAs(filename);
            delete c1;

            file->Close();
            file2->Close();
            file3->Close();
        }
    }

    
}

int main() {
    processRootFile();
    return 0;
}



