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
#include "include/sherpa_constants.h"

void Fourier_Analysis() {
    // Open the file
    TFile* file = new TFile("job_morebins_all.root");
    std::cout << "Opened file: job_morebins_all.root" << std::endl;
    // Loop over all objects in the file
    TIter next(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        // std::cout << "Tkey" << std::endl;
        TClass *cl = gROOT->GetClass(key->GetClassName());
        // if (!cl->InheritsFrom("TH2D")) continue;
        TH2D* h = (TH2D*)key->ReadObj();
         
        // Check if the histogram name matches the pattern
        std::string name = h->GetName();
        std::cout << "Histogram Name: " << name << std::endl;
        if (name.find("hjtw_SignalS") != 0) continue;
        std::cout << "hjtw_SignalS" << std::endl;
        // Extract the parameters from the histogram name
        int wtrk, wppt, wpPU;
        sscanf(name.c_str(), "hjtw_SignalS_trk_%d_ppt_%d_PU_%d", &wtrk, &wppt, &wpPU);
        std::cout << "Extracted parameters: " << wtrk << ", " << wppt << ", "  << wpPU << std::endl;
        // Perform the analysis for this histogram
        // ...
        TH2D* hSig = (TH2D*)file->Get(Form("hjtw_SignalS_trk_%d_ppt_%d_PU_%d" , wtrk, wppt, wpPU));
        TH2D* hBack = (TH2D*)file->Get(Form("hjtw_BckrndS_trk_%d_ppt_%d_PU_%d" , wtrk, wppt, wpPU));
        std::cout << Form("hjtw_SignalS_trk_%d_ppt_%d_PU_%d" , wtrk, wppt, wpPU) << std::endl;
        // Divide the signal histogram by the background histogram
        hSig->Divide(hBack);

        // Create a 1D projection for DeltaEta > 2
        double DelEta_lo = 2.0;
        TH1D* h1D = hSig->ProjectionY("h1D", hSig->GetXaxis()->FindBin(DelEta_lo), -1);
        
        // Define a 3 term harmonic function
        // double normalizationFactor = h1D->GetMaximum();
        TF1* func = new TF1("func", "[0]/(TMath::Pi()*2)*([1]+2*([2]*TMath::Cos(x)+[3]*TMath::Cos(2*x)+[4]*TMath::Cos(3*x)))", h1D->GetXaxis()->GetXmin(), h1D->GetXaxis()->GetXmax());
        // func->FixParameter(0, normalizationFactor);
        func->FixParameter(1, 1);
        // Fit the function to the histogram
        h1D->Fit(func, "NQ");

        // Get the parameters
        double normalizationFactor = func->GetParameter(0);
        double v1 = func->GetParameter(2);
        double v2 = func->GetParameter(3);
        double v3 = func->GetParameter(4);

        

        // Create the functions
        TF1* func1 = new TF1("v1", "[0]*TMath::Cos(x)", h1D->GetXaxis()->GetXmin(), h1D->GetXaxis()->GetXmax());
        // std:: cout<< v1 << std::endl;
        func1->SetParameter(0, v1);
        TF1* func2 = new TF1("v2", "[0]*TMath::Cos(2*x)", h1D->GetXaxis()->GetXmin(), h1D->GetXaxis()->GetXmax());
        // std:: cout<< v2 << std::endl;
        func2->SetParameter(0, v2);
        TF1* func3 = new TF1("v3", "[0]*TMath::Cos(3*x)", h1D->GetXaxis()->GetXmin(), h1D->GetXaxis()->GetXmax());
        func3->SetParameter(0, v3);
        // std:: cout<< v3 << std::endl;
       
        

        // Normalize func
        // double scaleParam = func->GetParameter(0);
        func->SetParameter(0, M_PI);
        // std::cout << scaleParam / normalizationFactor <<std::endl;
        // std::cout << func->GetParameter(1) <<std::endl;
        
        func->SetParameter(1, 0);
        h1D->Scale( M_PI / normalizationFactor);
        // Reduce the size of the error bars
        for (int i = 1; i <= h1D->GetNbinsX(); i++) {
            double binError = h1D->GetBinError(i);
            double binContent = h1D->GetBinContent(i);
            h1D->SetBinContent(i, binContent - 0.5);
            h1D->SetBinError(i, binError / normalizationFactor *M_PI);
            }       
        
        // Draw the histogram and the fit
        TCanvas* c1 = new TCanvas("c1", "Fourier_fit", 800, 600);
        h1D->SetMarkerStyle(20);
        h1D->SetMarkerColor(kRed);
        h1D->SetLineColor(kRed);
        h1D->Draw();
        func->SetLineColor(kBlack);
        func->Draw("same");
        
        // Adjust the y-axis range
        double yMin = TMath::Min(h1D->GetMinimum(), TMath::Min(TMath::Min(func1->GetMinimum(), func2->GetMinimum()), func3->GetMinimum()));
        double yMax = TMath::Max(h1D->GetMaximum(), TMath::Max(TMath::Max(func1->GetMaximum(), func2->GetMaximum()), func3->GetMaximum()));
        h1D->GetYaxis()->SetRangeUser(yMin, yMax);
        // Draw func1, func2, and func3
        func1->SetLineColor(kRed);
        func1->SetLineStyle(2);
        func1->SetLineWidth(2);
        func1->Draw("same");
        func2->SetLineColor(kGreen);
        func2->SetLineWidth(2);
        func2->SetLineStyle(2);
        func2->Draw("same");
        func3->SetLineColor(kBlue);
        func3->SetLineStyle(2);
        func2->SetLineWidth(2);
        func3->Draw("same");

        // Create a legend
        TLegend* leg = new TLegend(0.15,0.75,0.3,0.9);
        leg->AddEntry(h1D,"Original h1D","l");
        leg->AddEntry(func,"Fit function","l");
        leg->AddEntry(func1,"v1*cos(x)","l");
        leg->AddEntry(func2,"v2*cos(2x)","l");
        leg->AddEntry(func3,"v3*cos(3x)","l");
        leg->Draw();

        // Create the folder if it doesn't exist
        gSystem->mkdir("Fourier_fit", kTRUE);

        // Save the canvas
        // c1->SaveAs(Form("Fourier_fit/Fourier_fit_for_pt_%d_to_%d.png", ptbin_lo, ptbin_hi));
        TString filename = Form("Fourier_fit/Fourier_fit_for_track_%d_to%d_pt_%.2f_to_%.2f.png", trackbinbounds[wtrk-1], trackbinboundsUpper[wtrk-1], ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
        std::cout << "Saving canvas to: " << filename << std::endl;
        c1->SaveAs(filename);
        // Clean up
        delete c1;
        delete func;
        // delete normalized_func;
        delete func1;
        delete func2;
        delete func3;
        delete h1D;
        delete hSig;
        delete hBack;
    }
    std::cout << "Closing file." << std::endl;
    // Close the file
    delete file;
}

int main() {
    Fourier_Analysis();
    return 0;
}