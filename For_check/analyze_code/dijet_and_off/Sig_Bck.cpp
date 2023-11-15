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
#include "../include/sherpa_constants.h"
#include "TView.h"
// TH1::SetDefaultSumw2(kTRUE);
// TH2::SetDefaultSumw2(kTRUE);

void Fourier_analysis() {
    // Open the file
    TFile* file = new TFile("newroot/dijet/AB_MC_2.root");
    // std::cout << "Opened file: dijobs.root" << std::endl;
    // Loop over all objects in the file
    TIter next(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        TH2D* h = (TH2D*)key->ReadObj();
         
        std::string name = h->GetName();
        if (name.find("hSignalS") != 0) continue;
        int wtrk=0, wppt=0, wpPU=0;
        sscanf(name.c_str(), "hSignalS_trk_%d_ppt_%d_PU_%d", &wtrk, &wppt, &wpPU);
        std::cout << "Extracted parameters: " << wtrk << ", " << wppt << ", "  << wpPU << std::endl;

        TH2D* hSig = (TH2D*)file->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
        TH2D* hBack = (TH2D*)file->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
        // hSig->GetXaxis()->SetRangeUser(-6, -2);
        // hBack->GetXaxis()->SetRangeUser(-6, -2);
        hSig->SetDefaultSumw2(kTRUE);
        hBack->SetDefaultSumw2(kTRUE);
        hSig->Divide(hBack);
        // hSig->Scale(hBack->GetMaximum());
        TH2D* h1D = hSig;
        
       
        

        TCanvas* c1 = new TCanvas("c1", "signal", 800, 600);
        TString graphtitle = Form("signal_AB_MC for_track_%d_to_%d_pt_%.2f_to_%.2f.png", trackbinbounds[wtrk-1], trackbinboundsUpper[wtrk-1], ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
        h1D->Draw("surf7");
        h1D->SetTitle(graphtitle);
        h1D->GetXaxis()->SetTitle("#Delta#eta");
        c1->Update();
        TView *view = c1->GetView();
        if (view) {
            // Int_t irep;
            view->RotateView(60, 70);
            c1->SetView(view);
            c1->Modified();
            c1->Update();
        } else {
            printf("Invalid TView\n");
        }
        gSystem->mkdir("output_buffer/signals/AB_MC", kTRUE);

        TString filename = "output_buffer/signals/AB_MC/"+ graphtitle;
        c1->SaveAs(filename);
        // Clean up
        delete c1;
    }
    // std::cout << "Closing file." << std::endl;
    // Close the file
    delete file;
}

int main() {
    Fourier_analysis();
    return 0;
}
