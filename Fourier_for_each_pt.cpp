
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
#include "include/sherpa_constants.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "math.h"

void processRootFile() {
   

    
    TFile* file = new TFile("job_morebins_all_2.root");
    std::cout << "Opened file: job_morebins_all_2.root" << std::endl;
    TIter next(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
    // Loop over ptbin and trackbin
        TClass *cl = gROOT->GetClass(key->GetClassName());
        // if (!cl->InheritsFrom("TH2D")) continue;
        if (cl->InheritsFrom("TH2D")) {
        std::cout << "Histogram: " << key->GetName() << std::endl;}
        TH2D* h = (TH2D*)key->ReadObj();
         
        // Check if the histogram name matches the pattern
        std::string name = h->GetName();
        std::cout << "Histogram Name: " << name << std::endl;
        if (name.find("hjtw_SignalS") != 0) continue;
        std::cout << "hjtw_SignalS" << std::endl;
        for (int wppt = 1; wppt < ptbin+1; wppt++) {
            TString v_for_pt = Form("v_for_pt:%.2f_to_%.2f",ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
            // Variables to hold the v1, v2, v3 values
            TCanvas* c1 = new TCanvas("c1", v_for_pt, 800, 600);

            // Create a graph for each array
            TGraph* g1 = new TGraph(trackbin);
            TGraph* g2 = new TGraph(trackbin);
            TGraph* g3 = new TGraph(trackbin);
            int point = 0;
        
            for (int wtrk = 1; wtrk < trackbin+1; wtrk++) {
                // Get the histogram
                TH2D* hSig = (TH2D*)file->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH2D* hBack = (TH2D*)file->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                hSig->Divide(hBack);
                // Create a 1D projection for DeltaEta > 2
                double DelEta_lo = 2.0;
                TH1D* h1D = hSig->ProjectionY("h1D", hSig->GetXaxis()->FindBin(DelEta_lo), -1);
                
                // Define a 3 term harmonic function
                double normalizationFactor = h1D->GetMaximum();
                TF1* func = new TF1("func", "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))", h1D->GetXaxis()->GetXmin(), h1D->GetXaxis()->GetXmax());
                func->SetParameter(0, normalizationFactor);

                // Fit the function to the histogram
                h1D->Fit(func, "NQ");
                // std::cout << "fit" << std::endl;
                // Get the parameters
                g1->SetPoint(point, trackbinbounds[wtrk-1], func->GetParameter(1));
                g2->SetPoint(point, trackbinbounds[wtrk-1], func->GetParameter(2));
                g3->SetPoint(point, trackbinbounds[wtrk-1], func->GetParameter(3));
                point++;
                // std::cout << "draw" << std::endl;
            
            }

            TLegend* legend = new TLegend(0.15,0.75,0.3,0.9); // Adjust these parameters to move the legend to the desired position
            // Add entries to the legend
            legend->AddEntry(g1, "v1", "l");
            legend->AddEntry(g2, "v2", "l");
            legend->AddEntry(g3, "v3", "l");
            // Draw the graphs on the same canvas

            double ymin = std::min({g1->GetHistogram()->GetMinimum(), g2->GetHistogram()->GetMinimum(), g3->GetHistogram()->GetMinimum()});
            double ymax = std::max({g1->GetHistogram()->GetMaximum(), g2->GetHistogram()->GetMaximum(), g3->GetHistogram()->GetMaximum()});

            // Add a small buffer to the y range
            double buffer = 0.1 * (ymax - ymin);
            ymin -= buffer;
            ymax += buffer;

            // Set the y range for each graph
            g1->GetYaxis()->SetRangeUser(ymin, ymax);
            g2->GetYaxis()->SetRangeUser(ymin, ymax);
            g3->GetYaxis()->SetRangeUser(ymin, ymax);
            g1->SetTitle(v_for_pt);
            g1->GetXaxis()->SetTitle("Multiplicity");
            g1->SetLineColor(kBlue);
            g1->Draw("AL");
            g2->SetLineColor(kGreen);
            g2->Draw("L same");
            g3->SetLineColor(kMagenta);
            g3->Draw("L same");
            legend->Draw();
            c1->Update();

            gSystem->mkdir("Fourier_fit_each", kTRUE);

            TString filename = Form("Fourier_fit_each/Fourier_components_for_pt_%.2f_to_%.2f.png",  ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
            std::cout << "Saving canvas to: " << filename << std::endl;
            c1->SaveAs(filename);
            delete c1;
            delete g1;
            delete g2;
            delete g3;
        }
        
        
    }
    // outFile->Close();
    
    file->Close();
}

int main() {
    processRootFile();
    return 0;
}