
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
void processRootFile() {
   

    
   
        for (int wppt = 1; wppt < ptbin+1; wppt++) {
            
            TString v_for_pt = Form("v_for_pt:%.2f_to_%.2f",ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
            // Variables to hold the v1, v2, v3 values
            TCanvas* c1 = new TCanvas("c1", v_for_pt, 800, 600);

            // Create a graph for each array
            TGraphErrors* g1 = new TGraphErrors(trackbin);
            TGraphErrors* g2 = new TGraphErrors(trackbin);
            TGraphErrors* g3 = new TGraphErrors(trackbin);
            TGraphErrors* g1_2 = new TGraphErrors(trackbin);
            TGraphErrors* g2_2 = new TGraphErrors(trackbin);
            TGraphErrors* g3_2 = new TGraphErrors(trackbin);
            TGraphErrors* g1_3 = new TGraphErrors(trackbin);
            TGraphErrors* g2_3 = new TGraphErrors(trackbin);
            TGraphErrors* g3_3 = new TGraphErrors(trackbin);
            int point = 0;
            int point_2 = 0;
            int point_3 = 0;
            for (int wtrk = 1; wtrk < trackbin+1; wtrk++) {
                TFile* file = new TFile("jobs_more_all.root");
                // // // Get the histogram
                TH2D* hSig = (TH2D*)file->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH2D* hBack = (TH2D*)file->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                hSig->SetDefaultSumw2(kTRUE);
                hBack->SetDefaultSumw2(kTRUE);

                double DelEta_lo = 2.0;
                TH1D *hSig1D = (TH1D*) hSig->ProjectionY("h1D", hSig->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
                TH1D *hBack1D = (TH1D*) hBack->ProjectionY("h1D", hBack->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
                hSig1D->SetDefaultSumw2(kTRUE);
                hBack1D->SetDefaultSumw2(kTRUE);
                hSig1D->Divide(hBack1D);
                hSig1D->Scale(hBack1D->GetMaximum());
                // // hBack1D->Scale(hBack1D->GetMaximum());
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
                // std::cout<< "fit1"<< std::endl;
                // std::cout<< func->GetParameter(1)<< std::endl;
                // //    std::cout<< func->GetParameter(1)<< std::endl;
                // //     std::cout<< "fit1"<< std::endl;


                // // Fit the function to the histogram
                // // h1D->Fit(func, "NQ");
                // // std::cout << "fit" << std::endl;
                // // Get the parameters
                g1->SetPoint(point, trackbinbounds[wtrk-1], func->GetParameter(1));
                g1->SetPointError(point, 0, func->GetParError(1));
                
                g2->SetPoint(point, trackbinbounds[wtrk-1], func->GetParameter(2));
                g2->SetPointError(point, 0, func->GetParError(2));
                g3->SetPoint(point, trackbinbounds[wtrk-1], func->GetParameter(3));
                g3->SetPointError(point, 0, func->GetParError(3));
                point++;
                    // std::cout << "draw" << std::endl;
                file->Close();
            
                TFile* file2 = new TFile("jobs_more_all_2.root");
                // // std::cout << "Opened file: job_morebins_all_2.root" << std::endl;
                
                
            
                // // for (int wtrk2 = 1; wtrk2 < trackbin+1; wtrk2++) {
                //     // Get the histogram
                TH2D* hSig2 = (TH2D*)file2->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH2D* hBack2 = (TH2D*)file2->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                hSig2->SetDefaultSumw2(kTRUE);
                hBack2->SetDefaultSumw2(kTRUE);

                // // double DelEta_lo = 2.0;
                TH1D *hSig1D_2 = (TH1D*) hSig2->ProjectionY("h1D", hSig2->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
                TH1D *hBack1D_2 = (TH1D*) hBack2->ProjectionY("h1D", hBack2->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
                hSig1D_2->SetDefaultSumw2(kTRUE);
                hBack1D_2->SetDefaultSumw2(kTRUE);
                hSig1D_2->Divide(hBack1D_2);
                hSig1D_2->Scale(hBack1D_2->GetMaximum());
                // // hBack1D->Scale(hBack1D_2->GetMaximum());
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
                std::cout<< "fit2"<< std::endl;
                std::cout<< func_2->GetParameter(1)<< std::endl;
                // // std::cout<< "fit2"<< std::endl;
                // // func_2->GetParameter(1);
                // // Get the parameters
                g1_2->SetPoint(point_2, trackbinbounds[wtrk-1], func_2->GetParameter(1));
                g1_2->SetPointError(point_2, 0, func_2->GetParError(1));
                g2_2->SetPoint(point_2, trackbinbounds[wtrk-1], func_2->GetParameter(2));
                g2_2->SetPointError(point_2, 0, func_2->GetParError(2));
                g3_2->SetPoint(point_2, trackbinbounds[wtrk-1], func_2->GetParameter(3));
                g3_2->SetPointError(point_2, 0, func_2->GetParError(3));
                point_2++;
                file2->Close();





                
                TFile* file3 = new TFile("jobs_more_all_3.root");
                // std::cout << "Opened file: job_morebins_all_2.root" << std::endl;
                
                
            
                // for (int wtrk2 = 1; wtrk2 < trackbin+1; wtrk2++) {
                    // Get the histogram
                TH2D* hSig3 = (TH2D*)file3->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                TH2D* hBack3 = (TH2D*)file3->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1", trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1], (int)(10*ptbinbounds_lo[wppt-1]), (int)(10*ptbinbounds_hi[wppt-1])));
                hSig3->SetDefaultSumw2(kTRUE);
                hBack3->SetDefaultSumw2(kTRUE);

                // double DelEta_lo = 2.0;
                TH1D *hSig1D_3 = (TH1D*) hSig3->ProjectionY("h1D", hSig3->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
                TH1D *hBack1D_3 = (TH1D*) hBack3->ProjectionY("h1D", hBack3->GetXaxis()->FindBin(DelEta_lo), -1)->Clone();
                if(!hSig1D_3){
                    std::cout<< "null3"<< std::endl;
                }
                hSig1D_3->SetDefaultSumw2(kTRUE);
                hBack1D_3->SetDefaultSumw2(kTRUE);
                hSig1D_3->Divide(hBack1D_3);
                hSig1D_3->Scale(hBack1D_3->GetMaximum());
                // hBack1D->Scale(hBack1D_2->GetMaximum());
                TH1D* h1D_3 = hSig1D_3;
                h1D_3->SetDefaultSumw2(kTRUE);
                
                
                TF1* func_3 = new TF1("func_3", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                func_3->SetParameter(0, h1D_3->GetMaximum());
                func_3->SetParameter(1, 0.1);
               
                func_3->SetParameter(2, 0.1);
                func_3->SetParameter(3, 0.1);
                std::cout<< "fit3 before"<< std::endl;
                h1D_3->Fit(func_3, "q 0");
                h1D_3->Fit(func_3, "q 0");
                h1D_3->Fit(func_3, "m q 0");
                h1D_3->Fit(func_3, "m q 0");
                h1D_3->Fit(func_3, "m q E 0");
                h1D_3->Fit(func_3, "m E q 0");
                std::cout<< "fit3"<< std::endl;
                std::cout<< func_3->GetParameter(1)<< std::endl;
                // func_2->GetParameter(1);
                // Get the parameters
                g1_3->SetPoint(point_3, trackbinbounds[wtrk-1], func_3->GetParameter(1));
                g1_3->SetPointError(point_3, 0, func_3->GetParError(1));
                g2_3->SetPoint(point_3, trackbinbounds[wtrk-1], func_3->GetParameter(2));
                g2_3->SetPointError(point_3, 0, func_3->GetParError(2));
                g3_3->SetPoint(point_3, trackbinbounds[wtrk-1], func_3->GetParameter(3));
                g3_3->SetPointError(point_3, 0, func_3->GetParError(3));
                point_3++;
                file3->Close();


                // std::cout << "draw" << std::endl;
            }
            // }

            // TLegend* legend = new TLegend(0.15,0.75,0.3,0.9); // Adjust these parameters to move the legend to the desired position
            // Add entries to the legend
            





            
            TLegend* legend = new TLegend(0.7,0.65,0.85,0.9); // Adjust these parameters to move the legend to the desired position
                // Add entries to the legend
            legend->AddEntry(g1_2, "v1_jtnor", "l");
            legend->AddEntry(g2_2, "v2_jtnor", "l");
            legend->AddEntry(g3_2, "v3_jtnor", "l");

            legend->AddEntry(g1, "v1_jt", "l");
            legend->AddEntry(g2, "v2_jt", "l");
            legend->AddEntry(g3, "v3_jt", "l");

            legend->AddEntry(g1_3, "v1", "l");
            legend->AddEntry(g2_3, "v2", "l");
            legend->AddEntry(g3_3, "v3", "l");
            // Draw the graphs on the same canvasg1_2->GetHistogram()->GetMinimum(), g2_2->GetHistogram()->GetMinimum(), g3_2->GetHistogram()->GetMinimum()
            // g1_2->GetHistogram()->GetMinimum(), g2_2->GetHistogram()->GetMinimum(), g3_2->GetHistogram()->GetMinimum()
            double ymin = std::min({g1->GetYaxis()->GetXmin(), g2->GetYaxis()->GetXmin(), g3->GetYaxis()->GetXmin(), g1_2->GetYaxis()->GetXmin(), g2_2->GetYaxis()->GetXmin(), g3_2->GetYaxis()->GetXmin(), g1_3->GetYaxis()->GetXmin(), g2_3->GetYaxis()->GetXmin(), g3_3->GetYaxis()->GetXmin()});
            double ymax = std::max({g1->GetYaxis()->GetXmax(), g2->GetYaxis()->GetXmax(), g3->GetYaxis()->GetXmax(), g1_2->GetYaxis()->GetXmax(), g2_2->GetYaxis()->GetXmax(), g3_2->GetYaxis()->GetXmax(), g1_3->GetYaxis()->GetXmin(), g2_3->GetYaxis()->GetXmin(), g3_3->GetYaxis()->GetXmin()});
            // // Add a small buffer to the y range
            double buffer = 0.3 * (ymax - ymin);
            ymin -= buffer;
            ymax += buffer;

            // // Set the y range for each graph
            g1_2->GetYaxis()->SetRangeUser(ymin, ymax);
            g2_2->GetYaxis()->SetRangeUser(ymin, ymax);
            g3_2->GetYaxis()->SetRangeUser(ymin, ymax);
            g1->GetYaxis()->SetRangeUser(ymin, ymax);
            g2->GetYaxis()->SetRangeUser(ymin, ymax);
            g3->GetYaxis()->SetRangeUser(ymin, ymax);
            g1_3->GetYaxis()->SetRangeUser(ymin, ymax);
            g2_3->GetYaxis()->SetRangeUser(ymin, ymax);
            g3_3->GetYaxis()->SetRangeUser(ymin, ymax);
            g1_2->SetTitle(v_for_pt);
            g1_2->GetXaxis()->SetTitle("Multiplicity");
            g1_2->SetLineColor(kRed);
            g1_2->SetLineStyle(2);
            // g1_2->SetLineWidth(3); 
            g1_2->Draw("ACP");
            g2_2->SetLineColor(kBlue);
            g2_2->SetLineStyle(2);
            // g2_2->SetLineWidth(3); 
            g2_2->Draw("CP same");
            g3_2->SetLineColor(kGreen);
            g3_2->SetLineStyle(2);
            // g3_2->SetLineWidth(3);
            g3_2->Draw("CP same");
            g1->SetLineColor(kRed);
            g1->Draw("CP same");
            g2->SetLineColor(kBlue);
            g2->Draw("CP same");
            g3->SetLineColor(kGreen);
            g3->Draw("CP same");
            g1_3->SetLineColor(kRed);
            g1_3->SetLineStyle(8);
            g1_3->Draw("CP same");
            g2_3->SetLineColor(kBlue);
            g2_3->SetLineStyle(8);
            g2_3->Draw("CP same");
            g3_3->SetLineColor(kGreen);
            g3_3->SetLineStyle(8);
            g3_3->Draw("CP same");


            legend->Draw();
            c1->Update();

            gSystem->mkdir("Fourier_fit_each", kTRUE);

            TString filename = Form("Fourier_fit_each/Fourier_components_for_pt_%.2f_to_%.2f.png",  ptbinbounds_lo[wppt-1], ptbinbounds_hi[wppt-1]);
            // std::cout << "Saving canvas to: " << filename << std::endl;
            c1->SaveAs(filename);
            delete c1;
            delete g1;
            delete g2;
            delete g3;
            delete g1_2;
            delete g2_2;
            delete g3_2;
            delete g1_3;
            delete g2_3;
            delete g3_3;
            

        }





            
    // }

        
        
    // }// outFile->Close();
    
    
}

int main() {
    processRootFile();
    return 0;
}



