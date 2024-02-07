// #ifndef BOOTSTRAP
// #define BOOTSTRAP
// # define N 1000
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
#include "../include/sherpa_constants_cumulant.h"
#include "math.h"
#include "TGraphErrors.h" 
#include "TLine.h"
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <fstream>
#include <filesystem>
#include <thread>
#include <chrono>
#include <ctime>
// The index is for parallel jobs.
// Randomly hadd n files from new_cumulant, generate temp_%02d.root
// For MC, n = 19; For Data, n =80
#define INJ 80

void haddRandomFiles(int n, int index) {
    namespace fs = std::filesystem;
    std::string dir = Form("newroot/new_cumulant/0.0-3.0/3sub/%d", INJ);
    std::vector<std::string> files;

    // Get all .root files in the directory
    for (const auto &entry : fs::directory_iterator(dir)) {
        if (entry.path().extension() == ".root") {
            files.push_back(entry.path().string());
        }
    }
    // Resample:
    // seed from the current time and the thread ID
    unsigned int time_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::seed_seq ss{time_seed, static_cast<unsigned int>(std::hash<std::thread::id>{}(std::this_thread::get_id()))};
    std::mt19937 g(ss);

    std::uniform_int_distribution<> dis(0, files.size() - 1);

    std::vector<std::string> resampled_files;
    for(int i = 0; i < n; ++i) {
        int random_index = dis(g);
        resampled_files.push_back(files[random_index]);
    }

    //hadd the resamples to calculate c24
    
    std::string command = Form("hadd temp_%02d.root ",index);
    for (int i = 0; i < n && i < resampled_files.size(); ++i) {
        command += resampled_files[i] + " ";
    }
    system(command.c_str());

}


// read temp.root, calculate c24 for one resample, and fill the hc24
//This is for 2sub
void c24_2sub_single(TH1D* hc24[trackbin], int index){

    TString tempfile = Form("temp_%02d.root",index);
    TFile* file = new TFile(tempfile);   

    TH1D* hjet_avg_numerator_two = (TH1D*)file->Get("hjet_avg_numerator_two");
    TH1D* hjet_avg_denominat_two = (TH1D*)file->Get("hjet_avg_denominat_two");
    hjet_avg_numerator_two->SetDefaultSumw2(kTRUE);
    hjet_avg_denominat_two->SetDefaultSumw2(kTRUE);


    TH1D* hjet_avg_numerator_four = (TH1D*)file->Get("hjet_avg_numerator_four");
    TH1D* hjet_avg_denominat_four = (TH1D*)file->Get("hjet_avg_denominat_four");
    hjet_avg_numerator_four->SetDefaultSumw2(kTRUE);
    hjet_avg_denominat_four->SetDefaultSumw2(kTRUE);



    for (int bin = 1; bin <= trackbin; bin++) {

        double numerator_two = hjet_avg_numerator_two->GetBinContent(bin);
        double denominat_two = hjet_avg_denominat_two->GetBinContent(bin);
        // double x_mean = h_gen->GetMean();

        double c22 = numerator_two/denominat_two;

        double numerator_four = hjet_avg_numerator_four->GetBinContent(bin);
        double denominat_four = hjet_avg_denominat_four->GetBinContent(bin);

        double c4 = numerator_four/denominat_four;
        
        double c24 = c4 - 2 * c22 * c22;

    
        hc24[bin-1]->Fill(c24,1);
        
    }

    // Always remember to close the files, or errors occur
    file->Close();
    std::filesystem::remove(Form("temp_%02d.root",index));
}

//Three 

void c24_3sub_single(TH1D* hc24[trackbin], int index){

    TString tempfile = Form("temp_%02d.root",index);
    TFile* file = new TFile(tempfile);   

    TH1D* hjet_avg_numerator_two_AT = (TH1D*)file->Get("hjet_avg_numerator_two_AT");
    TH1D* hjet_avg_denominat_two_AT = (TH1D*)file->Get("hjet_avg_denominat_two_AT");
    hjet_avg_numerator_two_AT->SetDefaultSumw2(kTRUE);
    hjet_avg_denominat_two_AT->SetDefaultSumw2(kTRUE);

    TH1D* hjet_avg_numerator_two_AP = (TH1D*)file->Get("hjet_avg_numerator_two_AP");
    TH1D* hjet_avg_denominat_two_AP = (TH1D*)file->Get("hjet_avg_denominat_two_AP");
    hjet_avg_numerator_two_AP->SetDefaultSumw2(kTRUE);
    hjet_avg_denominat_two_AP->SetDefaultSumw2(kTRUE);


    TH1D* hjet_avg_numerator_four = (TH1D*)file->Get("hjet_avg_numerator_four");
    TH1D* hjet_avg_denominat_four = (TH1D*)file->Get("hjet_avg_denominat_four");
    hjet_avg_numerator_four->SetDefaultSumw2(kTRUE);
    hjet_avg_denominat_four->SetDefaultSumw2(kTRUE);
    
    for (int bin = 1; bin <= trackbin; bin++) {

        double numerator_two_AT = hjet_avg_numerator_two_AT->GetBinContent(bin);
        double denominat_two_AT = hjet_avg_denominat_two_AT->GetBinContent(bin);
        double c22_AT = numerator_two_AT/denominat_two_AT;

        double numerator_two_AP = hjet_avg_numerator_two_AP->GetBinContent(bin);
        double denominat_two_AP = hjet_avg_denominat_two_AP->GetBinContent(bin);
        double c22_AP = numerator_two_AP/denominat_two_AP;
        
        double numerator_four = hjet_avg_numerator_four->GetBinContent(bin);
        double denominat_four = hjet_avg_denominat_four->GetBinContent(bin);
        double c4 = numerator_four/denominat_four;
        
        double c24 = c4 - 2 * c22_AP * c22_AT ;
        hc24[bin-1]->Fill(c24,1);
        
    }

    // Always remember close the files, or errors occur
    file->Close();
    std::filesystem::remove(Form("temp_%02d.root",index));

}





// 4sub
void c24_4sub_single(TH1D* hc24[trackbin], int index){
    TString tempfile = Form("temp_%02d.root",index);
    TFile* file = new TFile(tempfile);  
    TH1D* hjet_avg_numerator_two_AP = (TH1D*)file->Get("hjet_avg_numerator_two_AP");
    TH1D* hjet_avg_denominat_two_AP = (TH1D*)file->Get("hjet_avg_denominat_two_AP");
    hjet_avg_numerator_two_AP->SetDefaultSumw2(kTRUE);
    hjet_avg_denominat_two_AP->SetDefaultSumw2(kTRUE);

    TH1D* hjet_avg_numerator_two_TQ = (TH1D*)file->Get("hjet_avg_numerator_two_TQ");
    TH1D* hjet_avg_denominat_two_TQ = (TH1D*)file->Get("hjet_avg_denominat_two_TQ");
    hjet_avg_numerator_two_TQ->SetDefaultSumw2(kTRUE);
    hjet_avg_denominat_two_TQ->SetDefaultSumw2(kTRUE);

    TH1D* hjet_avg_numerator_two_AQ = (TH1D*)file->Get("hjet_avg_numerator_two_AQ");
    TH1D* hjet_avg_denominat_two_AQ = (TH1D*)file->Get("hjet_avg_denominat_two_AQ");
    hjet_avg_numerator_two_AQ->SetDefaultSumw2(kTRUE);
    hjet_avg_denominat_two_AQ->SetDefaultSumw2(kTRUE);

    TH1D* hjet_avg_numerator_two_TP = (TH1D*)file->Get("hjet_avg_numerator_two_TP");
    TH1D* hjet_avg_denominat_two_TP = (TH1D*)file->Get("hjet_avg_denominat_two_TP");
    hjet_avg_numerator_two_TP->SetDefaultSumw2(kTRUE);
    hjet_avg_denominat_two_TP->SetDefaultSumw2(kTRUE);


    TH1D* hjet_avg_numerator_four = (TH1D*)file->Get("hjet_avg_numerator_four");
    TH1D* hjet_avg_denominat_four = (TH1D*)file->Get("hjet_avg_denominat_four");
    hjet_avg_numerator_four->SetDefaultSumw2(kTRUE);
    hjet_avg_denominat_four->SetDefaultSumw2(kTRUE);

    
    for (int bin = 1; bin <= trackbin; bin++) {

        double numerator_two_AP = hjet_avg_numerator_two_AP->GetBinContent(bin);
        double denominat_two_AP = hjet_avg_denominat_two_AP->GetBinContent(bin);
        double c22_AP = numerator_two_AP/denominat_two_AP;

        double numerator_two_TQ = hjet_avg_numerator_two_TQ->GetBinContent(bin);
        double denominat_two_TQ = hjet_avg_denominat_two_TQ->GetBinContent(bin);
        double c22_TQ = numerator_two_TQ/denominat_two_TQ;

        double numerator_two_AQ = hjet_avg_numerator_two_AQ->GetBinContent(bin);
        double denominat_two_AQ = hjet_avg_denominat_two_AQ->GetBinContent(bin);
        double c22_AQ = numerator_two_AQ/denominat_two_AQ;

        double numerator_two_TP = hjet_avg_numerator_two_TP->GetBinContent(bin);
        double denominat_two_TP = hjet_avg_denominat_two_TP->GetBinContent(bin);
        double c22_TP = numerator_two_TP/denominat_two_TP;
        
        double numerator_four = hjet_avg_numerator_four->GetBinContent(bin);
        double denominat_four = hjet_avg_denominat_four->GetBinContent(bin);
        double c4 = numerator_four/denominat_four;

        double c24 = c4 - c22_AP * c22_TQ - c22_AQ * c22_TP ;

        hc24[bin-1]->Fill(c24,1);
        
    }
    file->Close();
    std::filesystem::remove(Form("temp_%02d.root",index));
}
// N times bootstrap
void make_bootstrap(int index){

    gSystem->mkdir(Form("c24/0.0-3.0/3sub/%d", INJ), kTRUE);
    TFile *file = new TFile(Form("c24/0.0-3.0/3sub/%d/%02d.root", INJ,index), "RECREATE");
    TH1D* hc24[trackbin];


    for(int i=0; i<trackbin; i++){
        hc24[i] = new TH1D(Form("hc24_%d",i), Form("hc24_%d",i),10000,-0.01,0.01);  // 100 bins from 0 to 100
        hc24[i]-> SetDefaultSumw2(kTRUE);
    }
    // 300 per job, I usually take 10 jobs
    for(int i=0; i<1000; i++){
        haddRandomFiles(80, index);
        c24_3sub_single(hc24, index);
    }
    file->cd();
    for(int i=0; i<trackbin; i++){
        hc24[i] -> Write();  
    }

    file->Close();
}

int main(int argc, char* argv[]) {
    int index = std::atoi(argv[1]);
    make_bootstrap(index);
    return 0;
}






















































































// #endif