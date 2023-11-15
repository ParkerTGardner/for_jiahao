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
#include "../../include/sherpa_constants_cumulant.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "math.h"
#include "TGraphErrors.h" 
#include "TLine.h"

void processRootFile() {
    

    TString folder1 = "v2 = 0%";
    TString folder2 = "v2 = 5%";
    TString folder3 = "v2 = 10%";
    TString folder4 = "v2 = 15%";
    TString folder5 = "v2 = 20%";
    TString folder6 = "data";
    

    std::string folders[] = {"0","5","10", "15","20","80"};
    std::vector<TGraphErrors*> averages;
    for (const auto& folder : folders) {

        TString bootstrap_file = Form("c24/0.5-3.0/4sub/%s/all.root", folder.c_str());
        TFile* file = new TFile(bootstrap_file);
        
        TString add_all = Form("newroot/new_cumulant/0.5-3.0/4sub/%s/all_add.root", folder.c_str());
        TFile* file2 = new TFile(add_all);

        TGraphErrors* hc24 = new TGraphErrors(trackbin);
        hc24->SetStats(false);
        
        for (int bin = 1; bin <= trackbin; bin++) {

            TH1D* h_gen = (TH1D*)file2->Get(Form("hBinDist_gen_%d", bin));  // Get the corresponding hBinDist_gen histogram
            
            double x_mean = h_gen->GetMean();
            double x_mean_error = h_gen->GetRMS() / sqrt(h_gen->GetEntries());
            
            TH1D* hc24_bin = (TH1D*)file->Get(Form("hc24_%d",bin-1));
            double c24 = hc24_bin->GetMean();
            double c24_err = hc24_bin->GetStdDev();
            std::cout<<c24<<std::endl;
            hc24->SetPoint(bin-1, x_mean, c24);
            hc24->SetPointError(bin-1, x_mean_error, c24_err);
              
        }
        averages.push_back(hc24);
    }

    TCanvas* c = new TCanvas("c", "hc24", 1600, 1200);
    // TLegend* legend = new TLegend(0.12,0.75,0.4,0.9);
    TLegend* legend = new TLegend(0.6,0.75,0.88,0.9);
    legend->AddEntry(averages[0], folder1, "p");
    legend->AddEntry(averages[1], folder2, "p");
    legend->AddEntry(averages[2], folder3, "p");
    legend->AddEntry(averages[3], folder4, "p");
    legend->AddEntry(averages[4], folder5, "p");
    legend->AddEntry(averages[5], folder6, "p");
    // legend->AddEntry(averages[6], folder7+Form("       x = %.1f, y = %.5f #pm %.5f", x_9[6], y_9[6], err_9[6]), "p");
    
    double ymin = -0.005;
    double ymax = 0.003;

    // double ymin = -0.03;
    // double ymax = 0.02;
    
    averages[0]->GetYaxis()->SetRangeUser(ymin, ymax);
    averages[1]->GetYaxis()->SetRangeUser(ymin, ymax);
    averages[2]->GetYaxis()->SetRangeUser(ymin, ymax);
    averages[3]->GetYaxis()->SetRangeUser(ymin, ymax);
    averages[4]->GetYaxis()->SetRangeUser(ymin, ymax);
    averages[5]->GetYaxis()->SetRangeUser(ymin, ymax);

    


    averages[0]->SetLineColor(kRed);
    averages[0]->SetLineWidth(2);
    averages[0]->SetMarkerSize(2);
    averages[0]->SetMarkerStyle(20);
    averages[0]->SetMarkerColor(kRed);
    averages[0]->GetXaxis()->SetTitle("Multiplicity");

    // averages[0]->SetTitle("hc24 0.86, 2.25, 5.0 ");                //2sub
    // averages[0]->SetTitle("hc24 0.86, 1.8, 2.6, 5.0 ");            //3sub
    averages[0]->SetTitle("hc24 0.86, 1.45, 1.85, 2.3, 5.0 ");     //4sub
    
    averages[0]->Draw("APEC ");


    averages[1]->SetLineColor(kBlue);
    averages[1]->SetLineWidth(2);
    averages[1]->SetMarkerSize(2);
    averages[1]->SetMarkerStyle(22);
    averages[1]->SetMarkerColor(kBlue);
    averages[1]->Draw("PEC same");

    averages[2]->SetLineColor(kGreen);
    averages[2]->SetLineWidth(2);
    averages[2]->SetMarkerSize(2);
    averages[2]->SetMarkerStyle(45);
    averages[2]->SetMarkerColor(kGreen);
    averages[2]->Draw("PEC same");

    averages[3]->SetLineColor(kViolet);
    averages[3]->SetLineWidth(2);
    averages[3]->SetMarkerSize(2);
    averages[3]->SetMarkerStyle(43);
    averages[3]->SetMarkerColor(kViolet);
    averages[3]->Draw("PEC same");

    

    averages[4]->SetLineColor(kOrange-3);
    averages[4]->SetLineWidth(2);
    averages[4]->SetMarkerSize(2);
    averages[4]->SetMarkerStyle(45);
    averages[4]->SetMarkerColor(kOrange-3);
    averages[4]->Draw("PEC same");

    averages[5]->SetLineColor(kCyan+3);
    averages[5]->SetLineWidth(2);
    averages[5]->SetMarkerStyle(33);
    averages[5]->SetMarkerColor(kCyan+3);
    averages[5]->SetMarkerSize(3);
    averages[5]->Draw("PEC same");









    double x_min = averages[5]->GetXaxis()->GetXmin();
    // double x_max = g1->GetXaxis()->GetXmax();
    double x_max = averages[5]->GetX()[averages[5]->GetN()-1];

    TLine* line1 = new TLine(x_min, 0, x_max, 0);
    line1->SetLineColor(kBlack);
    line1->SetLineStyle(9);
    line1->Draw();

    TLine* line2 = new TLine(x_min, -0.00000625, x_max, -0.00000625);
    line2->SetLineColor(kBlue);
    line2->SetLineStyle(9);
    line2->Draw();

    TLine* line3 = new TLine(x_min, -0.0001, x_max, -0.0001);
    line3->SetLineColor(kGreen);
    line3->SetLineStyle(9);
    line3->Draw();

    TLine* line4 = new TLine(x_min, -0.00050625, x_max, -0.00050625);
    line4->SetLineColor(kViolet);
    line4->SetLineStyle(9);
    line4->Draw();

    TLine* line5 = new TLine(x_min, -0.0016, x_max, -0.0016);
    line5->SetLineColor(kOrange-3);
    line5->SetLineStyle(9);
    line5->Draw();

    
    
    legend->Draw();
    c->Update();


    gSystem->mkdir("out_new_cumulant/0.5-3.0", kTRUE);

    TString filename = Form("out_new_cumulant/0.5-3.0/hc24_4sub.png");

    c->SaveAs(filename);
    // delete file;
    delete c;
    averages.clear();
  
}

int main() {
    processRootFile();
    return 0;
}
