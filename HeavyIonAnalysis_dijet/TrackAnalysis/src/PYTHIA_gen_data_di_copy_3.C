#define MyClass_cxx
// touch



#include "include/TrimPythia_data.h"
//below allows for rotation into jet frame
#include "include/coordinateTools.h"
//defines constants
#include "include/sherpa_constants_2_2.h"

#include <iostream>
#include <iomanip>

#include <vector>
#include "math.h"
#include <numeric>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>

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

      //}}}

    //Initializing Histograms

    TH3D* hPairs        = new TH3D("hPairs","hPairs",  trackbin, bin0,trackbin, ptbin, bin0, ptbin, ptbin, bin0, ptbin);

    TH2D* hNtrig        = new TH2D("hNtrig","hNtrig",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);


    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", trackbin,bin0,trackbin);
    TH1D* hJet_Pass     = new TH1D("hJet_Pass"  ,"hJet_Pass"  , trackbin,bin0,trackbin);
    TH1D* hJetPt = new TH1D("hJetPt"  ,"hJetPt" ,i150,i100,i3500);
    TH1D* hBinDist_gen_single = new TH1D("hBinDist_gen_single","hBinDist_gen_single",bin360,bin0, bin120);
    TH1D* hBinDist_reco_single = new TH1D("hBinDist_reco_single","hBinDist_reco_single",bin360,bin0,bin120);

    TH1D* hBinDist_corrected = new TH1D("hBinDist_corrected","hBinDist_corrected",bin360,bin0, bin120);
    TH1D* hBinDist_uncorrected = new TH1D("hBinDist_uncorrected","hBinDist_uncorrected",bin360,bin0, bin120);

    


    TH2D* hEtaPhiA = new TH2D("hEtaPhiA","hEtaPhiA", 2*EPD_xb   , -EPD_xhi, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
    TH1D* hEtaA = new TH1D("hEtaA","hEtaA", 2*EPD_xb   , -EPD_xhi, EPD_xhi );
    TH1D* hPhiA = new TH1D("hPhiA","hPhiA",EPD_yb      , EPD_ylo    , EPD_yhi);
    TH1D* hJtA = new TH1D("hJtA","hJtA",60,0,3);
    // TH1D* hEtaJetA =new TH1D("hEtaJetA","hEtaJetA",50,-1,1);
    
   

    TH2D* hSignalShifted[trackbin][ptbin][ptbin];
    TH2D* hBckrndShifted[trackbin][ptbin][ptbin];
    TH2D* hEPDrawA[trackbin][ptbin];

    TH1D* hBinDist_gen[trackbin];
    TH1D* hBinDist_reco[trackbin];
    TH1D* hJt[trackbin][ptbin]; 
    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        hBinDist_gen[wtrk-1]    = new TH1D(Form("hBinDist_gen_%d",wtrk),Form("hBinDist_gen_%d",wtrk), bin360, bin0, bin120);
        hBinDist_reco[wtrk-1]   = new TH1D(Form("hBinDist_reco_%d",wtrk),Form("hBinewnDist_reco_%d",wtrk), bin360, bin0, bin120);
        
        for(int wppt = 1; wppt<ptbin+1; wppt++){
            hEPDrawA[wtrk-1][wppt-1]            = new TH2D(Form("hEPDrawA_trk_%d_ppt_%d",wtrk,wppt) ,Form( "hEPDrawA_trk_%d_ppt_%d",wtrk,wppt) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
            hJt[wtrk-1][wppt-1]     = new TH1D(Form("hJt_trk_%d_ppt_%d",wtrk,wppt),Form("hJt_trk_%d_ppt_%d",wtrk,wppt),60,0,3);
            for(int wppt2 = wppt; wppt2<ptbin+1; wppt2++){

                hBckrndShifted[wtrk-1][wppt-1][wppt2-1]      = new TH2D(Form("hBckrndS_trk_%d_ppt_%d_pt_%d",wtrk,wppt,wppt2) ,Form("hBckrndS_trk_%d_ppt_%d_pt_%d",wtrk,wppt,wppt2) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                hSignalShifted[wtrk-1][wppt-1][wppt2-1]      = new TH2D(Form("hSignalS_trk_%d_ppt_%d_pt_%d",wtrk,wppt,wppt2) ,Form("hSignalS_trk_%d_ppt_%d_pt_%d",wtrk,wppt,wppt2) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                
                

            }
        }
    }


  

    


  
    


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
        // TTree *tree = (TTree*)_file0->Get("analyzer500/trackTree");
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


            hEvent_Pass->Fill(1);


            
            
            //ENTERING JET LOOP

            //in this first loop I choose the jetA and count the trks of A 
            for(int kjet=0; kjet < jetPt->size(); kjet++){

                int ijet = kjet;
                int Gjet = kjet;
                long int NNtrk1 = (dau_pt->at(ijet)).size();

                TVector3 JetA;
                JetA.SetPtEtaPhi((*jetPt)[ijet],(*jetEta)[ijet],(*jetPhi)[ijet]);
                
                if( fabs(JetA.Eta()) > jetEtaCut ) continue;
                if( JetA.Perp() < jetPtCut_Jet   ) continue;


                 //    count the trks of A  
                // int n_G_ChargeMult_count  = 0;
                
                
                double jet_HLT_weight = 1.0;
                if((*jetPt)[ijet] < 880 && (*jetPt)[ijet] > 550){
                    jet_HLT_weight = 1.0/ (hdid500->GetBinContent(hdid500->FindBin((*jetPt)[ijet]) ) );
                }
                //(hdid500->GetBinContent(hdid500->FindBin((*jetPt)[ijet]) ) )

                 //    count the trks of A  
                double n_G_ChargeMult_count  = 0.0;
                // double n_G_ChargeMult_count1 = 0.0;
                int n_G_ChargeMult_count_corrected1 = 0.0;
                
                
                for(int  G_trk=0; G_trk < NNtrk1; G_trk++ ){
                    if((*dau_chg)[Gjet][G_trk] == 0) continue;
                    if(fabs((*dau_pt)[Gjet][G_trk])  < 0.3)     continue;
                    if(fabs((*dau_eta)[Gjet][G_trk]) > 2.4)     continue;
                    n_G_ChargeMult_count += 1;
                    double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][G_trk] , (*dau_eta)[ijet][    G_trk] )));
                    n_G_ChargeMult_count_corrected1 += (1.0/nUnc_weight);
                }

                hBinDist_corrected           ->Fill(n_G_ChargeMult_count_corrected1, 1.0*jet_HLT_weight);
                hBinDist_uncorrected        ->Fill(n_G_ChargeMult_count, 1.0*jet_HLT_weight);
                
                
                
            //    hBinDist_gen_single            ->Fill(n_G_ChargeMult_count);

                //some useful bools 
                // We have to do this here because we need to get correct mult first
                // and then get tktool, Ntrig

                int tkBool[trackbin] = {0};
                int Ntrig[trackbin][ptbin] = {0};
                int NtrigM[trackbin][ptbin] = {0};
                int NtrigP[trackbin][ptbin] = {0};
                int A_ptBool[NNtrk1][ptbin] = {0};    

                for(int i = 0; i < trackbin; i++){
                //if((*chargedMultiplicity)[indicesR[kjet]] >= trackbinbounds[i] && (*chargedMultiplicity)[indicesR[kjet]] < trackbinboundsUpper[i]){
                    if(n_G_ChargeMult_count >= trackbinbounds[i] && n_G_ChargeMult_count < trackbinboundsUpper[i]){
                        tkBool[i] = 1;
                        hJet_Pass           ->Fill(i);
                        hBinDist_gen[i]         ->Fill(n_G_ChargeMult_count);
                        
                    }
                }





                // calculate the A_ptbool pile up Ntrig in jetA first,
                //and then we do this in jetB so that we can get the complete Ntrig
                for(int  A_trk=0; A_trk < NNtrk1; A_trk++ ){
                
                    TVector3 dau_A;
                    dau_A.SetPtEtaPhi((double)(*dau_pt)[ijet][A_trk],(double)(*dau_eta)[ijet][A_trk],(double)(*dau_phi)[ijet][A_trk]);
                    
                    
                    if((*dau_chg)[ijet][A_trk] == 0) continue;
                    if(fabs(dau_A.Eta()) > 2.4)        continue;
                    if(fabs(dau_A.Perp())  < 0.3)      continue;

                    //     daughter pt with respect to the jet axis                 pt With Respect To Jet 
                    double jet_dau_pt    =  ptWRTJet(JetA, dau_A);

                    if(jet_dau_pt >3.0) continue;

                    double jet_dau_eta   = etaWRTJet(JetA, dau_A);
                    if(jet_dau_eta > track_eta_lim) continue;


                    // boosted:
                    // double jet_dau_pt    =  ptWRTJet(JetA, dau_A);

                    // find A_trk in i ptbin
                    for(int i = 0; i < ptbin; i++){
                        if(jet_dau_pt >= ptbinbounds_lo[i] && jet_dau_pt < ptbinbounds_hi[i]){
                            A_ptBool[A_trk][i] = 1;
                        }
                    }

                    double Atrk_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][    A_trk] )));

                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin; j++){
                            if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][    A_trk] )));
                                Ntrig[i][j] += (1.0/(nUnc_weight));
                                // Ntrig[i][j] += 1;
                                if((*dau_chg)[ijet][A_trk] > 0){
                                    NtrigP[i][j] += (1.0/(nUnc_weight));
                                }
                                if((*dau_chg)[ijet][A_trk] < 0){
                                    NtrigM[i][j] += (1.0/(nUnc_weight));
                                }
                            }
                        }
                    }

                        
                }

                // Here should be the final Ntrig for jetAB

                for(int i = 0; i < trackbin; i++){
                    for(int j = 0; j < ptbin; j++){
                        hNtrig->Fill(i,j,Ntrig[i][j]);
                    }
                }


                for(int  A_trk=0; A_trk < NNtrk1; A_trk++ ){
                    
                    TVector3 dau_A;
                    dau_A.SetPtEtaPhi((double)(*dau_pt)[ijet][A_trk],(double)(*dau_eta)[ijet][A_trk],(double)(*dau_phi)[ijet][A_trk]);
                    
                    if((*dau_chg)[ijet][A_trk] == 0) continue;
                    if(fabs(dau_A.Eta()) > 2.4)        continue;
                    if(fabs(dau_A.Perp())  < 0.3)      continue;

                    //     daughter pt with respect to the jet axis                 pt With Respect To Jet 
                    double jet_dau_pt    =  ptWRTJet(JetA, dau_A);

                    if(jet_dau_pt >3.0) continue;

        

                    //     daughter eta with respect to the jet axis                 eta With Respect To Jet 
                    double jet_dau_eta   = etaWRTJet(JetA, dau_A);
                    //     daughter phi with respect to the jet axis                 phi With Respect To Jet 
                    double jet_dau_phi   = phiWRTJet(JetA, dau_A) ;

                    if(jet_dau_eta > track_eta_lim) continue;

                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin; j++){
                            if(tkBool[i] + (jet_dau_pt >= ptbinbounds_lo[j] && jet_dau_pt < ptbinbounds_hi[j]) == 2){
                                hJt[i][j]->Fill(jet_dau_pt);
                            }
                        }
                    }




                    hEtaPhiA->Fill(jet_dau_eta, jet_dau_phi, 1);
                    hEtaA -> Fill(jet_dau_eta,1);
                    hPhiA -> Fill(jet_dau_phi,1);
                    hJtA  -> Fill(jet_dau_pt, 1);



                    double Atrk_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));

                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin; j++){
                            if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                int k_PU=0;
                                if ((Ntrig[i][j])==0) continue;
                                double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                                hEPDrawA[i][j]->Fill(jet_dau_eta, jet_dau_phi, 1.0* jet_HLT_weight/(Atrk_weight*Ntrig[i][j]));
                            }
                        }
                    }



                    

                    for(long int T_trk = A_trk+1; T_trk < NNtrk1; T_trk++ ){

                        
                        if((*dau_chg)[ijet][T_trk] == 0) continue;
                        if(fabs((*dau_eta)[ijet][T_trk]) > 2.4) continue;
                        if(fabs((*dau_pt)[ijet][T_trk])  < 0.3)      continue;

                            //Unboosted dau_T0
                        TVector3 dau_T;
                        dau_T.SetPtEtaPhi((double)(*dau_pt)[ijet][T_trk],(double)(*dau_eta)[ijet][T_trk],(double)(*dau_phi)[ijet][T_trk]);    
                        
                        double T_jet_dau_pt    =  ptWRTJet(JetA, dau_T);  
                        if(T_jet_dau_pt >3.0) continue;

                        double T_jet_dau_eta   = etaWRTJet(JetA, dau_T);
                        
                        if(T_jet_dau_eta > track_eta_lim) continue;




                        //boosted B wrt old A   
                        // double T_jet_dau_eta   = etaWRTJet(JetA, dau_T);
                        double T_jet_dau_phi   = phiWRTJet(JetA, dau_T);
                        // double T_jet_dau_pt    =  ptWRTJet(JetA, dau_T);
                        

                        //correlation function
                                            //A_trk(dau_A)  T_trk(dau_B)
                        double deltaEta = (jet_dau_eta - T_jet_dau_eta);
                                                                    //A_trk        T_trk
                        double deltaPhi = (TMath::ACos(TMath::Cos(jet_dau_phi - T_jet_dau_phi)));

                        double Ttrk_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][T_trk] , (*dau_eta)[ijet][    T_trk] )));

                        
                        // double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                        

                        // double deltaJt  = fabs(jet_dau_pt - T_jet_dau_pt);
                        
                        for(        int i = 0; i < trackbin; i++){
                            for(    int j = 0; j < ptbin;    j++){ 
                                for(int k = j; k < ptbin;    k++){ 
        


                                    if(tkBool[i] + A_ptBool[A_trk][j] + A_ptBool[T_trk][k] == 3){
                                            hPairs->Fill(i,j,k);
                                            if ((Ntrig[i][j])==0) continue;
                                            
                                            hSignalShifted[i][j][k]->Fill(deltaEta, deltaPhi,                 1.0* jet_HLT_weight/(Atrk_weight*Ttrk_weight*Ntrig[i][j]));
                                            hSignalShifted[i][j][k]->Fill(-deltaEta, deltaPhi,                1.0* jet_HLT_weight/(Atrk_weight*Ttrk_weight*Ntrig[i][j]));
                                            hSignalShifted[i][j][k]->Fill(deltaEta, -deltaPhi,                1.0* jet_HLT_weight/(Atrk_weight*Ttrk_weight*Ntrig[i][j]));
                                            hSignalShifted[i][j][k]->Fill(-deltaEta, -deltaPhi,               1.0* jet_HLT_weight/(Atrk_weight*Ttrk_weight*Ntrig[i][j]));
                                            hSignalShifted[i][j][k]->Fill( deltaEta,2*TMath::Pi() - deltaPhi, 1.0* jet_HLT_weight/(Atrk_weight*Ttrk_weight*Ntrig[i][j]));
                                            hSignalShifted[i][j][k]->Fill(-deltaEta,2*TMath::Pi() - deltaPhi, 1.0* jet_HLT_weight/(Atrk_weight*Ttrk_weight*Ntrig[i][j]));
                                            // hMomSignalShifted[i][j][k_PU]->Fill(deltaJt,                         1/(Ntrig[i][j]));

                                    }

                                }

                    
                            }
                        }

                    }//T_trk;  AA
                }// A_trk; AA
             
            }//kjet
                    }//Event
                    fFile->Close();
                    }//File
                    int backMult =10;
                    //filling background
                    //{{{
                    //

std::cout<< "made 4" << endl;
                    for(int wtrk = 1; wtrk < trackbin+1; wtrk++){
                              std::cout << wtrk << "/" << trackbin << std::endl;
                              for(int wppt = 1; wppt < ptbin+1; wppt++){
                                  std::cout << wppt << "/" << ptbin << std::endl;
                            

                                  
                                      
                                      //Nent is the number of pairs in the signal which we will try to 10x
                                      //Xent is the number of pseudoparticles requried such that when we build the pairs nCp = Xent CHOOSE 2 will give 
                                      //us 10 times as many pairs as we have in the signal histogrm.

                                      long int NENT0 =  hPairs->GetBinContent(wtrk, wppt, wppt);
                                      long int XENT0 =  ((1+floor(sqrt(1+(4*2*backMult*NENT0))))/2) ;
                                      double A_ETA0[XENT0] = {0};
                                      double A_PHI0[XENT0] = {0};
                                      


                                      
				                    //   float A_Jt[XENT]  = {0};

                                      for(int x = 0; x<XENT0; x++){
                                          gRandom->SetSeed(0);
                                          double WEtaA0, WPhiA0;//making the pseudoparticles
                                          
                                          hEPDrawA[wtrk-1][wppt-1]->GetRandom2(WEtaA0, WPhiA0);
                                          A_ETA0[x] = WEtaA0;
                                          A_PHI0[x] = WPhiA0; 
                                        //   A_Jt[x]  = WJt1;
                                      }
                                      
                                      for(long int i = 0; i < XENT0; i++){
                                          for(long int j = i+1; j < XENT0; j++){

                                              double WdeltaEta = (A_ETA0[i]-A_ETA0[j]);
                                              double WdeltaPhi = (TMath::ACos(TMath::Cos(A_PHI0[i]-A_PHI0[j])));
                    //                         //   double WdeltaJt  = fabs(A_Jt[i]-A_Jt[j]);

                                              hBckrndShifted[wtrk-1][wppt-1][wppt-1]->Fill(WdeltaEta, WdeltaPhi, 1);//./XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wppt-1]->Fill(-WdeltaEta, WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wppt-1]->Fill(WdeltaEta, -WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wppt-1]->Fill(-WdeltaEta, -WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wppt-1]->Fill(WdeltaEta, 2*TMath::Pi() - WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wppt-1]->Fill(-WdeltaEta,2*TMath::Pi() - WdeltaPhi, 1);//../XENT);
                    //                         //   hMomBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(WdeltaJt, 1);

                                          }
                                      }
                                  
                              




                                  for(int wppt2 = wppt+1; wppt2 < ptbin+1; wppt2++){
                                      std::cout << wppt2 << "/" << ptbin << std::endl;
                                      //Nent is the number of pairs in the signal which we will try to 10x
                                      //Xent is the number of pseudoparticles requried such that when we build the pairs nCp = Xent CHOOSE 2 will give 
                                      //us 10 times as many pairs as we have in the signal histogrm.

                                      long int NENT =  hPairs->GetBinContent(wtrk, wppt, wppt2);
                                      long int XENT =  floor(sqrt(10*NENT));
                                      double A_ETA[XENT] = {0};
                                      double A_PHI[XENT] = {0};
                                      double T_ETA[XENT] = {0};
                                      double T_PHI[XENT] = {0};


                                      
				                    //   float A_Jt[XENT]  = {0};

                                      for(int x = 0; x<XENT; x++){
                                          gRandom->SetSeed(0);
                                          double WEtaA, WPhiA;//making the pseudoparticles
                                          
                                          hEPDrawA[wtrk-1][wppt-1]->GetRandom2(WEtaA, WPhiA);
                                          A_ETA[x] = WEtaA;
                                          A_PHI[x] = WPhiA; 


                                          double WEtaT, WPhiT;  
                                          gRandom->SetSeed(gRandom->GetSeed() + 1);
                                          hEPDrawA[wtrk-1][wppt2-1]->GetRandom2(WEtaT, WPhiT);
                                          T_ETA[x] = WEtaT;
                                          T_PHI[x] = WPhiT;
                                        //   A_Jt[x]  = WJt1;
                                      }
                                      
                                      for(long int i = 0; i < XENT; i++){
                                          for(long int j = 0; j < XENT; j++){

                                              double WdeltaEta = (A_ETA[i]-T_ETA[j]);
                                              double WdeltaPhi = (TMath::ACos(TMath::Cos(A_PHI[i]-T_PHI[j])));
                    //                         //   double WdeltaJt  = fabs(A_Jt[i]-A_Jt[j]);

                                              hBckrndShifted[wtrk-1][wppt-1][wppt2-1]->Fill(WdeltaEta, WdeltaPhi, 1);//./XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wppt2-1]->Fill(-WdeltaEta, WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wppt2-1]->Fill(WdeltaEta, -WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wppt2-1]->Fill(-WdeltaEta, -WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wppt2-1]->Fill(WdeltaEta, 2*TMath::Pi() - WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wppt2-1]->Fill(-WdeltaEta,2*TMath::Pi() - WdeltaPhi, 1);//../XENT);
                    //                         //   hMomBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(WdeltaJt, 1);

                                          }
                                      }
                                  }
                         }
                    }
                        
                      //}}}


                      




                    string subList = fList.substr(fList.size() - 3);

                    TFile* fS_tempA = new TFile(Form("pythia_batch_data_output/root_out_copy_3/dijob_%s.root",subList.c_str()), "recreate");
                    
                    for(int wtrk =1; wtrk <trackbin+1; wtrk++){
                        hBinDist_gen[wtrk-1]         ->Write();
                        for(int wppt =1; wppt <ptbin+1; wppt++){
                            hEPDrawA                   [wtrk-1][wppt-1]->Write(Form("hEPDA_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1])));
                            hJt[wtrk-1][wppt-1] ->Write(Form("hJt_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1])));
                            for(int wppt2 = wppt; wppt2 <ptbin+1; wppt2++){

                                hSignalShifted             [wtrk-1][wppt-1][wppt2-1]->Write(Form("hSigS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),(int)(10*ptbinbounds_lo[wppt2-1]),(int)(10*ptbinbounds_hi[wppt2-1])    ));
                                hBckrndShifted             [wtrk-1][wppt-1][wppt2-1]->Write(Form("hBckS_%d_to_%d_pt_%d_to_%d_pt_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),(int)(10*ptbinbounds_lo[wppt2-1]),(int)(10*ptbinbounds_hi[wppt2-1])    ));
                                // hEPDraw                    [wtrk-1][wppt-1][wpPU-1]->Write(Form("hEPD_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU     ));
                                
                                                    
                            }
                        }
                    }


                    

                    hBinDist_gen_single->Write();
                    hEvent_Pass   ->Write();
                    hJet_Pass     ->Write();
                        
                    hEtaPhiA->Write();
                    hEtaA ->Write();
                    hPhiA ->Write();
                    hJtA ->Write();

                    fS_tempA->Close();
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