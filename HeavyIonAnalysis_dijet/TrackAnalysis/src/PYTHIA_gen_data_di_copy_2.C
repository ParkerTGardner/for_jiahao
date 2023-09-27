#define MyClass_cxx
// touch



#include "include/TrimPythia_data.h"
//below allows for rotation into jet frame
#include "include/coordinateTools.h"
//defines constants
#include "include/sherpa_constants.h"

#include <iostream>
#include <iomanip>

#include <vector>
#include "math.h"
#include <numeric>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <complex>
#include <cmath>

#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TComplex.h"


using TMath::ATan;
using TMath::ACos;
using TMath::Exp;



//Bool folds
//{{{
    bool F_eventpass(std::vector< float > *jetPt, int jetnumber, float jetPtCut){

        if(jetPt->size() < 1)           return false;
        if(jetnumber > 200)              return false;
        if((*jetPt)[0] < jetPtCut) return false;
        return true;
    }

    std::vector<int> F_JetSort(std::vector< float > * vphi, std::vector< float > * vphi2) {
        std::vector<int> result(vphi->size());
        float mR = *std::min_element(vphi->begin(), vphi->end());
        float mG = *std::min_element(vphi2->begin(), vphi2->end());
        float xR = *std::max_element(vphi->begin(), vphi->end());
        float xG = *std::max_element(vphi2->begin(), vphi2->end());
        if(fabs(xR)> 3.07 || fabs(mR) < .07 || fabs(xG) > 3.07 || fabs(mG) < .07 ){
            std::iota(result.begin(), result.end(), 0);
            std::sort(result.begin(), result.end(),
                    [&](int A, int B) -> bool {
                    return fabs((*vphi)[A]) > fabs((*vphi)[B]);
                    });
        }else{
            std::iota(result.begin(), result.end(), 0);
            std::sort(result.begin(), result.end(),
                    [&](int A, int B) -> bool {
                    return (*vphi)[A] > (*vphi)[B];
                    });
        }
        return result;
    }

    std::vector<int> F_JetSortMult(std::vector< int > * vmult) {
        std::vector<int> result(vmult->size());
        std::iota(result.begin(), result.end(), 0);
        std::sort(result.begin(), result.end(),
                [&](int A, int B) -> bool {
                return fabs((*vmult)[A]) > fabs((*vmult)[B]);
                });
        return result;
    }

    void  keepMatch(std::vector<int> keep, std::vector<int>& indices){
        for(int i = keep.size() - 1; i >= 0; i--){
            if(keep[i] < 2){
                keep.erase(keep.begin() + i); 
                indices.erase(indices.begin() + i);
            }
        }
    }

    bool F_jetpass(std::vector< float > * jetEta, std::vector< float > * jetPt, int     ijet, float   jetPtCut){
        if(fabs( (*jetEta)[ijet])   >jetEtaCut)   return false;
        if((*jetPt)[ijet]           <jetPtCut)    return false;
        return true;
    }
  //}}}

// THIS IS WHERE THE MAIN CODE BEGINS




void MyClass::Loop(int job, std::string fList){
// cout << "TEST 1" << endl;
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    bool PU_bool = 0;
    int PUbin;
    if(PU_bool==1){
        PUbin = 3;
    }else{
        PUbin = 1;
    }
    bool smear_bool = 0;



          double track_eta_lim = 5.0;
          const int i40 = 40;
          const int i470 = 470;
          const int i3500 = 3500;
          const int i50 = 50;
          const int i500 = 500;
          const int i1500 = 1500;
          const int i150 = 150;
          const int i100 = 100;

          TH1D * hleadingJetPt = new TH1D("leadingJetPt",";Leading p_{T}^{gen};#sigma (pb)",i50,i500,i1500);
          TH1D * hjetPtReco = new TH1D("JetPtReco",";p_{T}^{gen};#sigma (pb)",i150,i100,i3500);
          TH1D * hjetPtGen = new  TH1D("JetPtGen",";p_{T}^{gen};#sigma (pb)",i150,i100,i3500);



    TH1D* hBinDist_corrected = new TH1D("hBinDist_corrected","hBinDist_corrected",bin360,bin0, bin120);
    TH1D* hBinDist_uncorrected = new TH1D("hBinDist_uncorrected","hBinDist_uncorrected",bin360,bin0, bin120);



    TH1D* hjet_avg_numerator_two = new TH1D("hjet_avg_numerator_two", "hjet_avg_numerator_two" ,trackbin , trackbinEdge);
    TH1D* hjet_avg_denominat_two = new TH1D("hjet_avg_denominat_two", "hjet_avg_denominat_two" ,trackbin , trackbinEdge);
    TH1D* hjet_avg_numerator_four = new TH1D("hjet_avg_numerator_four", "hjet_avg_numerator_four" ,trackbin , trackbinEdge);
    TH1D* hjet_avg_denominat_four = new TH1D("hjet_avg_denominat_four", "hjet_avg_denominat_four" ,trackbin , trackbinEdge);

    
   
    TH1D* hBinDist_gen[trackbin];
    TH1D* hBinDist_gen_corrected[trackbin];
    TH1D* hPhiDrawA[trackbin];
    TH1D* hPhiDrawT[trackbin];
    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        hBinDist_gen[wtrk-1]    = new TH1D(Form("hBinDist_gen_%d",wtrk),Form("hBinDist_gen_%d",wtrk), bin360, bin0, bin120);
        hBinDist_gen_corrected[wtrk-1] = new TH1D(Form("hBinDist_gen_corrected_%d",wtrk),Form("hBinDist_gen_corrected_%d",wtrk), bin360, bin0, bin120);
        hPhiDrawA[wtrk-1] = new TH1D(Form("hPhiDrawA_%d",wtrk),Form("hPhiDrawA_%d",wtrk), 1000, -TMath::Pi(), TMath::Pi());
        hPhiDrawT[wtrk-1] = new TH1D(Form("hPhiDrawT_%d",wtrk),Form("hPhiDrawT_%d",wtrk), 1000, -TMath::Pi(), TMath::Pi());
        
        
    }

    double n_harm = 2.0;


  

    


  
    


    TH2D* hReco2D[fileList.size()];
    TH2D* hGen2D[fileList.size()];
    TH1D* hdid500;
    TH1D* hdid400;

    TFile *f_jet_HLT_lookup = new TFile("did400500_v2_all.root");
    hdid500 = (TH1D*)f_jet_HLT_lookup->Get("d500")->Clone("did500");
    hdid400 = (TH1D*)f_jet_HLT_lookup->Get("d400")->Clone("did400");
    hdid500->Divide(hdid400);




       std::cout << "Starting event loop" << std::endl;
    std::cout << "Total Number of Files in this Job: " << fileList.size() << std::endl;
    for(int f = 0; f<fileList.size(); f++){
        int f_from_file = f;
        fFile = TFile::Open(fileList.at(f).c_str(),"read");
        TTree *tree = (TTree*)fFile->Get("analyzer500/trackTree");
        Init(tree);

        std::cout << "File " << f+1 << " out of " << fileList.size() << std::endl;
        Long64_t nbytes = 0, nb = 0;
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"Total Entries is:"<<endl;
        cout<< nentries <<endl;

        vector<string> era_vec;
        vector<string> matched_cor_table_vec;
        era_vec.push_back("2018D_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2018C_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2018B_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2018A_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2017B_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017C_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017D_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017E_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017F_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2016Bv2_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016C_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016D_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016E_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016F-HIPM_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016F_"); matched_cor_table_vec.push_back("16_v2");
        era_vec.push_back("2016G_"); matched_cor_table_vec.push_back("16_v2");
        era_vec.push_back("2016H_"); matched_cor_table_vec.push_back("16_v2");
        int i_keep=999;
        for(int i=0; i<era_vec.size(); i++){
            if(fileList.at(f).find(era_vec[i]) != std::string::npos){
                i_keep=i;
                continue;
            }
        }

        string filename = "~/"+matched_cor_table_vec[i_keep]+".root";
        TFile *f_pt_eta_DCA_lookup = new TFile(filename.c_str());

        hReco2D[f] = (TH2D*)f_pt_eta_DCA_lookup->Get("h2_MAIN_DCA_Dau_Reco_Pt_Eta_Lab_All")->Clone(Form("h2_DCA_Dau_Reco_Pt_Eta_Lab_All_%d",f));
        hGen2D[f] = (TH2D*)f_pt_eta_DCA_lookup->Get("h2_Dau_Gen_Pt_Eta_Lab_All")->Clone(Form("h2_Dau_Gen_Pt_Eta_Lab_All_%d",f));
        hReco2D[f]->Divide(hGen2D[f]);
        int thisEffTable =f_from_file;
    


        // ENTERING EVENT LOOP
        // for(int f = 0; f<fileList.size(); f++){


        for (Long64_t ievent=0; ievent <nentries; ievent ++){
                  Long64_t jevent = LoadTree(ievent);
                  nb = fChain->GetEntry(ievent);   nbytes += nb;

//what is the def of genDau_pt, genJetPt
                  if(jetEta->size()==0) continue;
                  if(chargedMultiplicity->size()==0) continue;


                  if(!F_eventpass(jetPt, jetPt->size(), jetPtCut_Event)){
                      continue;
                  }
                  int gjN = jetPhi->size();


            // hEvent_Pass->Fill(1);


            
            
            
            
            
            //ENTERING JET LOOP

            //in this first loop I choose the jetA and count the trks of A 
            for(int kjet=0; kjet < jetPt->size(); kjet++){

                int ijet = kjet;
                int Gjet = kjet;
                long int NNtrk1 = (dau_pt->at(ijet)).size();

                TVector3 JetA;
                JetA.SetPtEtaPhi((*jetPt)[ijet],(*jetEta)[ijet],(*jetPhi)[ijet]);
                // TLorentzVector JetA_4 (JetA, JetA.Mag());
                
                if( fabs(JetA.Eta()) > jetEtaCut ) continue;
                if( JetA.Perp() < jetPtCut_Jet   ) continue;


                double jet_HLT_weight = 1.0;
                if((*jetPt)[ijet] < 880 && (*jetPt)[ijet] > 550){
                    jet_HLT_weight = 1.0/ (hdid500->GetBinContent(hdid500->FindBin((*jetPt)[ijet]) ) );
                }
                //(hdid500->GetBinContent(hdid500->FindBin((*jetPt)[ijet]) ) )

                 //    count the trks of A  
                int n_G_ChargeMult_count  = 0;
                int n_G_ChargeMult_count_corrected = 0;
                
                
                for(int  G_trk=0; G_trk < NNtrk1; G_trk++ ){
                    if((*dau_chg)[Gjet][G_trk] == 0) continue;
                    if(fabs((*dau_pt)[Gjet][G_trk])  < 0.3)     continue;
                    if(fabs((*dau_eta)[Gjet][G_trk]) > 2.4)     continue;
                    n_G_ChargeMult_count += 1;
                    double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][G_trk] , (*dau_eta)[ijet][    G_trk] )));
                    n_G_ChargeMult_count_corrected += (1.0/nUnc_weight);
                }

                hBinDist_corrected           ->Fill(n_G_ChargeMult_count_corrected, 1.0*jet_HLT_weight);
                hBinDist_uncorrected        ->Fill(n_G_ChargeMult_count, 1.0*jet_HLT_weight);
                
                int tkBool[trackbin] = {0};
                

                for(int i = 0; i < trackbin; i++){
                    if(n_G_ChargeMult_count >= trackbinbounds[i] && n_G_ChargeMult_count < trackbinboundsUpper[i]){
                        tkBool[i] = 1;
                        // hJet_Pass           ->Fill(i);
                        hBinDist_gen[i]         ->Fill(n_G_ChargeMult_count);
                    }
                }

                for(int i = 0; i < trackbin; i++){
                    if(n_G_ChargeMult_count_corrected >= trackbinbounds[i] && n_G_ChargeMult_count_corrected < trackbinboundsUpper[i]){
                        // tkBool[i] = 1;
                        // hJet_Pass           ->Fill(i);
                        hBinDist_gen_corrected[i]         ->Fill(n_G_ChargeMult_count_corrected);
                    }
                }



                std::complex<double> Q_all2A (0, 0);
                std::complex<double> Q_all4A (0, 0);
                std::complex<double> Q_all2T (0, 0);
                std::complex<double> Q_all4T (0, 0);

                double M = 0.0;
                double N = 0.0;

                // calculate the A_ptbool pile up Ntrig in jetA first,
                //and then we do this in jetB so that we can get the complete Ntrig
                for(int  A_trk=0; A_trk < NNtrk1; A_trk++ ){
                
                    TVector3 dau_A0;
                    dau_A0.SetPtEtaPhi((double)(*dau_pt)[ijet][A_trk],(double)(*dau_eta)[ijet][A_trk],(double)(*dau_phi)[ijet][A_trk]);
                    // TLorentzVector dau_A0_4(dau_A0,dau_A0.Mag());
                    
                    if((*dau_chg)[ijet][A_trk] == 0) continue;
                    if(fabs(dau_A0.Eta()) > 2.4)        continue;
                    if(fabs(dau_A0.Perp())  < 0.3)      continue;

                    //     daughter pt with respect to the jet axis                 pt With Respect To Jet 
                    double jet_dau_pt    =  ptWRTJet(JetA, dau_A0);

                    if(jet_dau_pt >3.0) continue;
                    // if(jet_dau_pt0 <0.3) continue;

                    // TLorentzVector dau_A_4 = BeamBoost(Boost_to_CM,dau_A0_4);
                    // TVector3       dau_A   = dau_A_4.Vect();

                    double jet_dau_eta   = etaWRTJet(JetA, dau_A0);
                    //     daughter phi with respect to the jet axis                 phi With Respect To Jet 
                    double jet_dau_phi   = phiWRTJet(JetA, dau_A0) ;

                    // double jet_dau_pt    =  ptWRTJet(JetA, dau_A0);



                    // if(jet_dau_pt >3.0) continue;
                    if(jet_dau_pt  <0.3) continue;

                    double Atrk_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][    A_trk] )));

                    

                    double phi = jet_dau_phi;


                    
                    std::complex<double> Q_part2 (TMath::Cos(n_harm*phi) , TMath::Sin(n_harm*phi));
                    std::complex<double> Q_part4 (TMath::Cos(2.0*n_harm*phi) , TMath::Sin(2.0*n_harm*phi));
                    if ((jet_dau_eta>0.86) && (jet_dau_eta<1.75)){

                        Q_all2A = Q_all2A + Q_part2;
                        Q_all4A = Q_all4A + Q_part4;
                        double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][    A_trk] )));
                        M += (1.0* jet_HLT_weight/(nUnc_weight*Atrk_weight));
                        // M++;
                    }

                    else if ((jet_dau_eta>2.75) && (jet_dau_eta<5.00)){

                        Q_all2T = Q_all2T + Q_part2;
                        Q_all4T = Q_all4T + Q_part4;
                       double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][    A_trk] )));
                        N += (1.0* jet_HLT_weight/(nUnc_weight*Atrk_weight));
                        // N++;
                    }
                    else continue;
                        

                    for(int i = 0; i < trackbin; i++){
                    //if((*chargedMultiplicity)[indicesR[kjet]] >= trackbinbounds[i] && (*chargedMultiplicity)[indicesR[kjet]] < trackbinboundsUpper[i]){
                        if(tkBool[i] == 1){
                            if (jet_dau_eta>0.86 && jet_dau_eta<1.75) hPhiDrawA[i]->Fill(phi);
                            if (jet_dau_eta>2.75 && jet_dau_eta<5.00) hPhiDrawT[i]->Fill(phi);
                        }
                    }

                        
                }



                // int M = n_G_ChargeMult_count ;
                if (M<3) continue;
                if (N<3) continue;

                double Q_all_absAT = std::real(Q_all2A*std::conj(Q_all2T));
                double particle_twoAT = Q_all_absAT / (M*N);

                double weight_two = (M*N);


                for(int i = 0; i < trackbin; i++){
                    if(tkBool[i] == 1){

                    hjet_avg_numerator_two->Fill(i,((weight_two)*(particle_twoAT)));
                    hjet_avg_denominat_two->Fill(i,(weight_two));
                    }
                }

                double particle_four = std::real( (Q_all2A*Q_all2A-Q_all4A)*std::conj((Q_all2T*Q_all2T-Q_all4T)) ) / (M*(M-1)*N*(N-1));
                double weight_four = M*(M-1)*N*(N-1);

                for(int i = 0; i < trackbin; i++){
                    if(tkBool[i] == 1){

                    hjet_avg_numerator_four->Fill(i,((weight_four)*(particle_four)) );
                    hjet_avg_denominat_four->Fill(i,(weight_four));
                    }
                }


                
             
            }//kjet
                    }//Event
                    fFile->Close();
                    }//File




    
    
                    
    TH1D* hRand_jet_avg_numerator_two = new TH1D("hRand_jet_avg_numerator_two", "hRand_jet_avg_numerator_two" ,trackbin , trackbinEdge);
    TH1D* hRand_jet_avg_denominat_two = new TH1D("hRand_jet_avg_denominat_two", "hRand_jet_avg_denominat_two" ,trackbin , trackbinEdge);
    TH1D* hRand_jet_avg_numerator_four = new TH1D("hRand_jet_avg_numerator_four", "hRand_jet_avg_numerator_four" ,trackbin , trackbinEdge);
    TH1D* hRand_jet_avg_denominat_four = new TH1D("hRand_jet_avg_denominat_four", "hRand_jet_avg_denominat_four" ,trackbin , trackbinEdge);


    for(int f = 0; f<fileList.size(); f++){
        int f_from_file = f;
        fFile = TFile::Open(fileList.at(f).c_str(),"read");
        TTree *tree = (TTree*)fFile->Get("analyzer500/trackTree");
        Init(tree);

        std::cout << "File " << f+1 << " out of " << fileList.size() << std::endl;
        Long64_t nbytes = 0, nb = 0;
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"Total Entries is:"<<endl;
        cout<< nentries <<endl;

        vector<string> era_vec;
        vector<string> matched_cor_table_vec;
        era_vec.push_back("2018D_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2018C_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2018B_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2018A_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2017B_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017C_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017D_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017E_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017F_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2016Bv2_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016C_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016D_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016E_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016F-HIPM_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016F_"); matched_cor_table_vec.push_back("16_v2");
        era_vec.push_back("2016G_"); matched_cor_table_vec.push_back("16_v2");
        era_vec.push_back("2016H_"); matched_cor_table_vec.push_back("16_v2");
        int i_keep=999;
        for(int i=0; i<era_vec.size(); i++){
            if(fileList.at(f).find(era_vec[i]) != std::string::npos){
                i_keep=i;
                continue;
            }
        }

        string filename = "~/"+matched_cor_table_vec[i_keep]+".root";
        TFile *f_pt_eta_DCA_lookup = new TFile(filename.c_str());

        hReco2D[f] = (TH2D*)f_pt_eta_DCA_lookup->Get("h2_MAIN_DCA_Dau_Reco_Pt_Eta_Lab_All")->Clone(Form("h2_DCA_Dau_Reco_Pt_Eta_Lab_All_%d",f));
        hGen2D[f] = (TH2D*)f_pt_eta_DCA_lookup->Get("h2_Dau_Gen_Pt_Eta_Lab_All")->Clone(Form("h2_Dau_Gen_Pt_Eta_Lab_All_%d",f));
        hReco2D[f]->Divide(hGen2D[f]);
        int thisEffTable =f_from_file;


        // ENTERING EVENT LOOP
        //for(int f = 0; f<fileList.size(); f++){


        for (Long64_t ievent=0; ievent <nentries; ievent ++){
                  Long64_t jevent = LoadTree(ievent);
                  nb = fChain->GetEntry(ievent);   nbytes += nb;

//what is the def of genDau_pt, genJetPt
                  if(jetEta->size()==0) continue;
                  if(chargedMultiplicity->size()==0) continue;


                  if(!F_eventpass(jetPt, jetPt->size(), jetPtCut_Event)){
                      continue;
                  }
                  int gjN = jetPhi->size();


            // hEvent_Pass->Fill(1);


            
            
            
            
            
            //ENTERING JET LOOP

            //in this first loop I choose the jetA and count the trks of A 
            for(int kjet=0; kjet < jetPt->size(); kjet++){

                int ijet = kjet;
                int Gjet = kjet;
                long int NNtrk1 = (dau_pt->at(ijet)).size();

                TVector3 JetA;
                JetA.SetPtEtaPhi((*jetPt)[ijet],(*jetEta)[ijet],(*jetPhi)[ijet]);
                // TLorentzVector JetA_4 (JetA, JetA.Mag());
                
                if( fabs(JetA.Eta()) > jetEtaCut ) continue;
                if( JetA.Perp() < jetPtCut_Jet   ) continue;


                double jet_HLT_weight = 1.0;
                if((*jetPt)[ijet] < 880 && (*jetPt)[ijet] > 550){
                    jet_HLT_weight = 1.0/(hdid500->GetBinContent(hdid500->FindBin((*jetPt)[ijet]) ) );
                }
                // /(hdid500->GetBinContent(hdid500->FindBin((*jetPt)[ijet]) ) )

                 //    count the trks of A  
                int n_G_ChargeMult_count  = 0;
                int n_G_ChargeMult_count_corrected = 0;
                
                
                for(int  G_trk=0; G_trk < NNtrk1; G_trk++ ){
                    if((*dau_chg)[Gjet][G_trk] == 0) continue;
                    if(fabs((*dau_pt)[Gjet][G_trk])  < 0.3)     continue;
                    if(fabs((*dau_eta)[Gjet][G_trk]) > 2.4)     continue;
                    n_G_ChargeMult_count += 1;
                    double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][G_trk] , (*dau_eta)[ijet][    G_trk] )));
                    n_G_ChargeMult_count_corrected += (1.0/nUnc_weight);
                }

                // hBinDist_corrected           ->Fill(n_G_ChargeMult_count_corrected, 1.0*jet_HLT_weight);
                // hBinDist_uncorrected        ->Fill(n_G_ChargeMult_count, 1.0*jet_HLT_weight);
                
                int tkBool[trackbin] = {0};
                

                for(int i = 0; i < trackbin; i++){
                    if(n_G_ChargeMult_count >= trackbinbounds[i] && n_G_ChargeMult_count < trackbinboundsUpper[i]){
                        tkBool[i] = 1;
                        // hJet_Pass           ->Fill(i);
                        // hBinDist_gen[i]         ->Fill(n_G_ChargeMult_count);
                    }
                }




                std::complex<double> Q_all2A (0, 0);
                std::complex<double> Q_all4A (0, 0);
                std::complex<double> Q_all2T (0, 0);
                std::complex<double> Q_all4T (0, 0);

                double M = 0.0;
                double N = 0.0;

                // calculate the A_ptbool pile up Ntrig in jetA first,
                //and then we do this in jetB so that we can get the complete Ntrig
                for(int  A_trk=0; A_trk < NNtrk1; A_trk++ ){
                
                    TVector3 dau_A0;
                    dau_A0.SetPtEtaPhi((double)(*dau_pt)[ijet][A_trk],(double)(*dau_eta)[ijet][A_trk],(double)(*dau_phi)[ijet][A_trk]);
                    // TLorentzVector dau_A0_4(dau_A0,dau_A0.Mag());
                    
                    if((*dau_chg)[ijet][A_trk] == 0) continue;
                    if(fabs(dau_A0.Eta()) > 2.4)        continue;
                    if(fabs(dau_A0.Perp())  < 0.3)      continue;

                    //     daughter pt with respect to the jet axis                 pt With Respect To Jet 
                    double jet_dau_pt    =  ptWRTJet(JetA, dau_A0);

                    if(jet_dau_pt >3.0) continue;
                    // if(jet_dau_pt0 <0.3) continue;

                    // TLorentzVector dau_A_4 = BeamBoost(Boost_to_CM,dau_A0_4);
                    // TVector3       dau_A   = dau_A_4.Vect();

                    double jet_dau_eta   = etaWRTJet(JetA, dau_A0);
                    //     daughter phi with respect to the jet axis                 phi With Respect To Jet 
                    double jet_dau_phi   = phiWRTJet(JetA, dau_A0) ;

                    // double jet_dau_pt    =  ptWRTJet(JetA, dau_A0);



                    // if(jet_dau_pt >3.0) continue;
                    if(jet_dau_pt  <0.3) continue;

                    double Atrk_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][    A_trk] )));

                    gRandom->SetSeed(0);
                    double phi;
                    if (jet_dau_eta>0.86 && jet_dau_eta<1.75){

                        for(int i = 0; i < trackbin; i++){
                            if(tkBool[i] == 1){
                                phi = hPhiDrawA[i]->GetRandom();
                            }
                        }

                        
                    }

                    else if (jet_dau_eta>2.75 && jet_dau_eta<5.00){

                        for(int i = 0; i < trackbin; i++){
                            if(tkBool[i] == 1){
                                phi = hPhiDrawT[i]->GetRandom();
                            }
                        }
                    }
                    else continue;


                    
                    std::complex<double> Q_part2 (TMath::Cos(n_harm*phi) , TMath::Sin(n_harm*phi));
                    std::complex<double> Q_part4 (TMath::Cos(2.0*n_harm*phi) , TMath::Sin(2.0*n_harm*phi));
                    if ((jet_dau_eta>0.86) && (jet_dau_eta<1.75)){

                        Q_all2A = Q_all2A + Q_part2;
                        Q_all4A = Q_all4A + Q_part4;
                        double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][    A_trk] )));
                        M += (1.0* jet_HLT_weight/(nUnc_weight*Atrk_weight));
                        // M++;
                    }

                    else if ((jet_dau_eta>2.75) && (jet_dau_eta<5.00)){

                        Q_all2T = Q_all2T + Q_part2;
                        Q_all4T = Q_all4T + Q_part4;
                        double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][    A_trk] )));
                        N += (1.0* jet_HLT_weight/(nUnc_weight*Atrk_weight));
                        // N++
                    }
                    else continue;
                           
                }



                // int M = n_G_ChargeMult_count ;
                if (M<3) continue;
                if (N<3) continue;

                double Q_all_absAT = std::real(Q_all2A*std::conj(Q_all2T));
                double particle_twoAT = Q_all_absAT / (M*N);

                double weight_two = (M*N);


                for(int i = 0; i < trackbin; i++){
                    if(tkBool[i] == 1){

                    hRand_jet_avg_numerator_two -> Fill(i,((weight_two)*(particle_twoAT))); 
                    hRand_jet_avg_denominat_two -> Fill(i,(weight_two));
                    }
                }

                double particle_four = std::real( (Q_all2A*Q_all2A-Q_all4A)*std::conj((Q_all2T*Q_all2T-Q_all4T)) ) / (M*(M-1)*N*(N-1));
                double weight_four = M*(M-1)*N*(N-1);

                for(int i = 0; i < trackbin; i++){
                    if(tkBool[i] == 1){

                    hRand_jet_avg_numerator_four -> Fill (i,(weight_four)*(particle_four));
                    hRand_jet_avg_denominat_four -> Fill (i,weight_four);
                    }
                }


                
             
            }//kjet
                    }//Event
                    fFile->Close();
                    }//File

    string subList = fList.substr(fList.size() - 3);
    TFile* fS_tempA = new TFile(Form("pythia_batch_data_output/root_out_copy_2/dijob_%s.root",subList.c_str()), "recreate");
    
    hBinDist_corrected->Write();
    hBinDist_uncorrected->Write();
    hjet_avg_numerator_two ->Write();
    hjet_avg_denominat_two ->Write();
    hjet_avg_numerator_four->Write();
    hjet_avg_denominat_four ->Write();
    
    hRand_jet_avg_numerator_two ->Write();
    hRand_jet_avg_denominat_two ->Write();
    hRand_jet_avg_numerator_four->Write();
    hRand_jet_avg_denominat_four ->Write();

    for(int i=0; i<trackbin; i++){
        hBinDist_gen[i]->Write();
        hBinDist_gen_corrected[i]->Write();
    }
    fS_tempA->Close();
    // fS_2_temp->Close();

}


                    //Code enters execution here
                    int main(int argc, const char* argv[])
                    {
                        if(argc != 4)
                        {
                            std::cout << "Usage: Z_mumu_Channel <fileList> <jobNumber> <nJobs>" << std::endl;
                            return 1;
                        }


                        //read input parameters
                        std::string fList = argv[1];
                        std::string buffer;
                        std::vector<std::string> listOfFiles;
                        std::ifstream inFile(fList.data());

                        int job = (int)std::atoi(argv[2]);
                        int nJobs = (int)std::atoi(argv[3]);


                        //read the file list and spit it into a vector of strings based on how the parallelization is to be done
                        //each vector is a separate subset of the fileList based on the job number
                        if(!inFile.is_open())
                        {
                            std::cout << "Error opening jet file. Exiting." <<std::endl;
                            return 1;
                        }
                        else
                        {
                            int line = 0;
                            while(true)
                            {
                                inFile >> buffer;
                                if(inFile.eof()) break;
                                if( line%nJobs == job) listOfFiles.push_back(buffer);
                                line++;
                            }
                        }

                        //create the MyClass Object
                        MyClass m = MyClass(listOfFiles);
                        m.Loop(job, fList);

                        return 0;
                    }