#define MyClass_cxx
// touch



#include "include/TrimPythia.h"
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

      //}}}

    //Initializing Histograms


    // int hist_bins = mult_bins.size();
    double mult_binsA[trackbin]; 

    long double jet_avg_numerator_two[trackbin] = {0};
    long double jet_avg_denominat_two[trackbin] = {0};
    double mult_bin_avg_two[trackbin] = {0};

    long double jet_avg_numerator_four[trackbin] = {0};
    long double jet_avg_denominat_four[trackbin] = {0};
    double mult_bin_avg_four[trackbin] = {0};

    TH1D* hPhiDrawA = new TH1D("hPhiDrawA", "hPhiDrawA", 1000, -TMath::Pi(), TMath::Pi());

    TH1D* hc22    = new TH1D("hc22", "hc22" ,trackbin-1 , mult_binsA);
    TH1D* hc24    = new TH1D("hc24", "hc24" ,trackbin-1 , mult_binsA);
    TH1D* hv22    = new TH1D("hv22", "hc22" ,trackbin-1 , mult_binsA);
    TH1D* hv24    = new TH1D("hv24", "hc24" ,trackbin-1 , mult_binsA);

    TH2D* hPairs        = new TH2D("hPairs","hPairs",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);

    TH2D* hNtrig        = new TH2D("hNtrig","hNtrig",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);


    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", trackbin,bin0,trackbin);
    TH1D* hJet_Pass     = new TH1D("hJet_Pass"  ,"hJet_Pass"  , trackbin,bin0,trackbin);
    TH1D* hJetPt = new TH1D("hJetPt"  ,"hJetPt" ,i150,i100,i3500);
    TH1D* hBinDist_gen_single = new TH1D("hBinDist_gen_single","hBinDist_gen_single",bin360,bin0, bin120);
    TH1D* hBinDist_reco_single = new TH1D("hBinDist_reco_single","hBinDist_reco_single",bin360,bin0,bin120);
    TH1D* hdeltaR = new TH1D("hdeltaR", "hdeltaR", 240,0,2.4);
    TH1D* hdeltaJetPhi = new TH1D("hdeltaJetPhi", "hdeltaJetPhi",50,3.0/4.0*M_PI,M_PI);
    TH1D* hdeltaJetTheta = new TH1D("hdeltaJetTheta", "hdeltaJetTheta",50,-M_PI,M_PI);
    TH1D* hdeltaJetEta = new TH1D("hdeltaJetEta", "hdeltaJetEta",50,-4,4);

    TH1D* hJJT1D                = new TH1D(Form("hJJT1D") ,Form("hJJT1D") , 120, 0, 1.2);
    TH2D* hJJT              = new TH2D(Form("hJJT") ,Form("hJJT") , 120, 0, 1.2, 140, 550,700);
    TH1D* hMult_ratio_AB        = new TH1D("hMult_ratio_AB", "hMult_ratio_AB", 250, 0, 5);
    TH2D* hMult_AB              = new TH2D("hMult_AB", "hMult_AB", 140, 0, 140, 140, 0, 140);
    TH2D* hJT_Mult_AB           = new TH2D("hJT_Mult_AB","hJT_Mult_AB",120, 0, 1.2, 200,0,2);




    TH1D* hEtaT = new TH1D("hEtaT","hEtaT",100,-5,5);
    TH1D* hPhiT = new TH1D("hPhiT","hPhiT",100,-M_PI,M_PI);
    TH1D* hJtT  = new TH1D("hJtT","hJtT",100,0,10);
    TH2D* hEtaPhiT = new TH2D("hEtaPhiT","hEtaPhiT", 2*EPD_xb   , -EPD_xhi, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);

    TH1D* hEtaTT = new TH1D("hEtaTT","hEtaTT",100,-5,5);
    TH1D* hPhiTT = new TH1D("hPhiTT","hPhiTT",100,-M_PI,M_PI);
    TH1D* hJtTT  = new TH1D("hJtTT","hJtTT",100,0,10);
    TH2D* hEtaPhiTT = new TH2D("hEtaPhiTT","hEtaPhiTT", 2*EPD_xb   , -EPD_xhi, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);

    TH2D* hEtaPhiA = new TH2D("hEtaPhiA","hEtaPhiA", 2*EPD_xb   , -EPD_xhi, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
    TH1D* hEtaA = new TH1D("hEtaA","hEtaA", 2*EPD_xb   , -EPD_xhi, EPD_xhi );
    TH1D* hPhiA = new TH1D("hPhiA","hPhiA",EPD_yb      , EPD_ylo    , EPD_yhi);
    TH1D* hJtA = new TH1D("hJtA","hJtA",100,0,10);
    // TH1D* hEtaJetA =new TH1D("hEtaJetA","hEtaJetA",50,-1,1);
    
   

    TH1D* hBinDist_gen[trackbin];
    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        hBinDist_gen[wtrk-1]    = new TH1D(Form("hBinDist_gen_%d",wtrk),Form("hBinDist_gen_%d",wtrk), bin360, bin0, bin120);
    }

    double n_harm = 2.0;


  

    


  
    


    std::cout << "Starting event loop" << std::endl;
    std::cout << "Total Number of Files in this Job: " << fileList.size() << std::endl;


//main loops
    for(int f = 0; f<fileList.size(); f++){
//processing data from CMS
        fFile = TFile::Open(fileList.at(f).c_str(),"read");
        TTree *tree = (TTree*)fFile->Get("trackTree");
        Init(tree);

        std::cout << "File " << f+1 << " out of " << fileList.size() << std::endl;
        Long64_t nbytes = 0, nb = 0;
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"Total Entries is:"<<endl;
        cout<< nentries <<endl;
std::cout << "File is " << fileList.at(f).c_str() << endl;


        // ENTERING EVENT LOOP
        //for(int f = 0; f<fileList.size(); f++){


        for (Long64_t ievent=0; ievent <nentries; ievent ++){
                  Long64_t jevent = LoadTree(ievent);
                  nb = fChain->GetEntry(ievent);   nbytes += nb;

//what is the def of genDau_pt, genJetPt
                  if(genJetPt->size()==0) continue;
                  if(genJetChargedMultiplicity->size()==0) continue;


                  if(!F_eventpass(genJetPt, genJetPt->size(), jetPtCut_Event)){
                      continue;
                  }
                  int gjN = genJetPhi->size();


            hEvent_Pass->Fill(1);


            
            
            
            
            
            //ENTERING JET LOOP

            //in this first loop I choose the jetA and count the trks of A 
            for(int kjet=0; kjet < genJetPt->size(); kjet++){

                int ijet = kjet;
                int Gjet = kjet;
                long int NNtrk1 = (genDau_pt->at(ijet)).size();

                TVector3 JetA;
                JetA.SetPtEtaPhi((*genJetPt)[ijet],(*genJetEta)[ijet],(*genJetPhi)[ijet]);
                TLorentzVector JetA_4 (JetA, JetA.Mag());
                
                if( fabs(JetA.Eta()) > jetEtaCut ) continue;
                if( JetA.Perp() < jetPtCut_Jet   ) continue;


                 //    count the trks of A  
                int n_G_ChargeMult_count  = 0;
                int n_G_ChargeMult_count1 = 0;
                
                
                for(int  G_trk=0; G_trk < NNtrk1; G_trk++ ){
                    if((*genDau_chg)[Gjet][G_trk] == 0) continue;
                    if(fabs((*genDau_pt)[Gjet][G_trk])  < 0.3)     continue;
                    if(fabs((*genDau_eta)[Gjet][G_trk]) > 2.4)     continue;
                    n_G_ChargeMult_count1 += 1;
                }
                // if (n_G_ChargeMult_count1<30) continue;
                
                
                for(int jjet=ijet+1; (jjet< genJetPt->size()); jjet++){


                    TVector3 JetB;
                    JetB.SetPtEtaPhi((*genJetPt)[jjet],(*genJetEta)[jjet],(*genJetPhi)[jjet]);
                    TLorentzVector JetB_4 (JetB, JetB.Mag());

                    if( fabs(JetB.Eta()) > jetEtaCut ) continue;
                    if( JetB.Perp() < jetPtCut_Jet-200 ) continue;
                    // if( JetB.Perp() < JetA.Perp()*0.95 ) continue;
              
                    TLorentzVector Boost_to_CM = JetA_4 + JetB_4;
                    TLorentzVector JetAA_4 = BeamBoost(Boost_to_CM,JetA_4);
                    TLorentzVector JetBB_4 = BeamBoost(Boost_to_CM,JetB_4);

                    TVector3 JetAA = JetAA_4.Vect();
                    TVector3 JetBB = JetBB_4.Vect();

                    // TVector3 JetAB = BeamBoost(JetA.Perp(),JetA.Eta(),JetA.Phi(),JetB.Perp(),JetB.Eta(),JetB.Phi());
                    double deltaJetEta = JetAA.Eta() + JetBB.Eta();
                    double deltaJetPhi =  fabs(JetA.Phi()-JetB.Phi());
                    hdeltaJetPhi -> Fill(deltaJetPhi);
                    hdeltaJetEta -> Fill(deltaJetEta);

                    if (fabs(M_PI-deltaJetPhi) > 0.1) continue;
                    // if (fabs(deltaJetEta)>0.15) continue;
                    long int NNtrk2 = (genDau_pt->at(jjet)).size();
                    // hdeltaR -> Fill(deltaR);

                   
                    
                    // Calculate the trks in jetB
                    int n_G_ChargeMult_count2 = 0;
                    for(int G_trk2=0; G_trk2 < NNtrk2; G_trk2++ ){
                        if((*genDau_chg)[jjet][G_trk2] == 0) continue;
                        if(fabs((*genDau_pt)[jjet][G_trk2])  < 0.3)     continue;
                        if(fabs((*genDau_eta)[jjet][G_trk2]) > 2.4)     continue;
                        n_G_ChargeMult_count2 += 1;
                    }

                    

                    // n_G_ChargeMult_count = n_G_ChargeMult_count1 + n_G_ChargeMult_count2 ;
                    // n_G_ChargeMult_count = ((1+floor(sqrt(1+(4*2*n_G_ChargeMult_count1*n_G_ChargeMult_count2))))/2) ;
                    
                    n_G_ChargeMult_count = n_G_ChargeMult_count1;
                    hBinDist_gen_single            ->Fill(n_G_ChargeMult_count1);


                    int tkBool[trackbin] = {0};
                    int Ntrig[trackbin][ptbin] = {0};
                    int NtrigM[trackbin][ptbin] = {0};
                    int NtrigP[trackbin][ptbin] = {0};
                    int A_ptBool[NNtrk1][ptbin] = {0};    //
                    int T_ptBool[NNtrk2][ptbin]     = {0};// This is for the AB

                    for(int i = 0; i < trackbin; i++){
                    //if((*chargedMultiplicity)[indicesR[kjet]] >= trackbinbounds[i] && (*chargedMultiplicity)[indicesR[kjet]] < trackbinboundsUpper[i]){
                        if(n_G_ChargeMult_count >= trackbinbounds[i] && n_G_ChargeMult_count < trackbinboundsUpper[i]){
                            tkBool[i] = 1;
                            // hJet_Pass           ->Fill(i);
                            hBinDist_gen[i]         ->Fill(n_G_ChargeMult_count);
                        }
                    }



                    std::complex<double> Q_all2 (0, 0);
                    std::complex<double> Q_all4 (0, 0);

                    // calculate the A_ptbool pile up Ntrig in jetA first,
                    //and then we do this in jetB so that we can get the complete Ntrig
                    for(int  A_trk=0; A_trk < NNtrk1; A_trk++ ){
                    
                        TVector3 dau_A0;
                        dau_A0.SetPtEtaPhi((double)(*genDau_pt)[ijet][A_trk],(double)(*genDau_eta)[ijet][A_trk],(double)(*genDau_phi)[ijet][A_trk]);
                        TLorentzVector dau_A0_4(dau_A0,dau_A0.Mag());
                        
                        if((*genDau_chg)[ijet][A_trk] == 0) continue;
                        if(fabs(dau_A0.Eta()) > 2.4)        continue;
                        if(fabs(dau_A0.Perp())  < 0.3)      continue;

                        //     daughter pt with respect to the jet axis                 pt With Respect To Jet 
                        double jet_dau_pt0    =  ptWRTJet(JetA, dau_A0);

                        if(jet_dau_pt0 >3.0) continue;
                        if(jet_dau_pt0 <0.3) continue;

                        TLorentzVector dau_A_4 = BeamBoost(Boost_to_CM,dau_A0_4);
                        TVector3       dau_A   = dau_A_4.Vect();

                        double phi = dau_A.Phi();


                        
                        std::complex<double> Q_part2 (TMath::Cos(n_harm*phi) , TMath::Sin(n_harm*phi));
                        Q_all2 = Q_all2 + Q_part2;

                        std::complex<double> Q_part4 (TMath::Cos(2.0*n_harm*phi) , TMath::Sin(2.0*n_harm*phi));
                        Q_all4 = Q_all4 + Q_part4;

                        hPhiDrawA->Fill(phi);



                            
                    }


                    int M = NNtrk1;

                    double Q_all_abs = std::abs(Q_all2);
                    double Q_all_sqr = Q_all_abs * Q_all_abs;
                    double particle_two = (Q_all_sqr - M) / (M*(M-1));
                    int weight_two = (M*(M-1));


                    for(int i = 0; i < trackbin; i++){
                        if(tkBool[i] == 1){

                        jet_avg_numerator_two[i] = jet_avg_numerator_two[i] + ((weight_two)*(particle_two)); 
                        jet_avg_denominat_two[i] = jet_avg_denominat_two[i] + (weight_two);
                        }
                    }

                    double Q_all_fourth = Q_all_sqr * Q_all_sqr;            
                    double Q_all_abs_2n = std::abs(Q_all4);   
                    double Q_all_sqr_2n = Q_all_abs_2n * Q_all_abs_2n;      
                    double Q_all_re_threetrm = std::real(
                            (2.0 * Q_all4) *
                            std::conj( Q_all2 ) *
                            std::conj( Q_all2 )
                            );                                              
                    double particle_four =
                        (   (Q_all_fourth + Q_all_sqr_2n - 2.0 * Q_all_re_threetrm)
                            / ( M*(M-1)*(M-2)*(M-3) )
                        )
                        -
                        (   2.0* ( (2.0*(M-2)*Q_all_sqr) - (M*(M-3)) )
                            / ( M*(M-1)*(M-2)*(M-3) )
                        );
                    int weight_four = M*(M-1)*(M-2)*(M-3);

                    for(int i = 0; i < trackbin; i++){
                        if(tkBool[i] == 1){

                        jet_avg_numerator_four[i] = jet_avg_numerator_four[i] + ((weight_four)*(particle_four));
                        jet_avg_denominat_four[i] = jet_avg_denominat_four[i] + (weight_four);
                        }
                    }


                }//jjet
             
            }//kjet
                    }//Event
                    fFile->Close();
                    }//File

    long double c22[trackbin] = {0};
    long double c24[trackbin] = {0};
    long double v22[trackbin] = {0};
    long double v24[trackbin] = {0};

    for(int i = 1; i < trackbin; i++){
        mult_bin_avg_two[i]  =jet_avg_numerator_two[i]/jet_avg_denominat_two[i];
        mult_bin_avg_four[i] =jet_avg_numerator_four[i]/jet_avg_denominat_four[i];

        c22[i] = mult_bin_avg_two[i];
        c24[i] = mult_bin_avg_four[i] - (2* mult_bin_avg_two[i] * mult_bin_avg_two[i]);

        v22[i] = std::real(std::pow(std::complex<double>(c22[i]),.50));
        v24[i] = std::real(std::pow(std::complex<double>(-c24[i]),.25));


        hc22->Fill(mult_bins[i]+.5, c22[i]);
        hc24->Fill(mult_bins[i]+.5, c24[i]);
        hv22->Fill(mult_bins[i]+.5, v22[i]);
        hv24->Fill(mult_bins[i]+.5, v24[i]);
    }
    TFile* fS_tempA = new TFile(Form("pythia_batch_output/root_out_2/dijob_%s.root",subList.c_str()), "recreate");
    hc22->Write();
    hc24->Write();
    hv22->Write();
    hv24->Write();
    fS_temp->Close();
                    
    long double Rand_jet_avg_numerator_two[trackbin] = {0};
    long double Rand_jet_avg_denominat_two[trackbin] = {0};
    double Rand_mult_bin_avg_two[trackbin] = {0};

    long double Rand_jet_avg_numerator_four[trackbin] = {0};
    long double Rand_jet_avg_denominat_four[trackbin] = {0};
    double Rand_mult_bin_avg_four[trackbin] = {0};

    TH1D* hRand_c22    = new TH1D("hRand_c22", "hRand_c22" ,trackbin-1 , mult_binsA);
    TH1D* hRand_c24    = new TH1D("hRand_c24", "hRand_c24" ,trackbin-1 , mult_binsA);
    TH1D* hRand_v22    = new TH1D("hRand_v22", "hRand_c22" ,trackbin-1 , mult_binsA);
    TH1D* hRand_v24    = new TH1D("hRand_v24", "hRand_c24" ,trackbin-1 , mult_binsA);

   for(int f = 0; f<fileList.size(); f++){
//processing data from CMS
        fFile = TFile::Open(fileList.at(f).c_str(),"read");
        TTree *tree = (TTree*)fFile->Get("trackTree");
        Init(tree);

        std::cout << "File " << f+1 << " out of " << fileList.size() << std::endl;
        Long64_t nbytes = 0, nb = 0;
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"Total Entries is:"<<endl;
        cout<< nentries <<endl;
std::cout << "File is " << fileList.at(f).c_str() << endl;


        // ENTERING EVENT LOOP
        //for(int f = 0; f<fileList.size(); f++){


        for (Long64_t ievent=0; ievent <nentries; ievent ++){
                  Long64_t jevent = LoadTree(ievent);
                  nb = fChain->GetEntry(ievent);   nbytes += nb;

//what is the def of genDau_pt, genJetPt
                  if(genJetPt->size()==0) continue;
                  if(genJetChargedMultiplicity->size()==0) continue;


                  if(!F_eventpass(genJetPt, genJetPt->size(), jetPtCut_Event)){
                      continue;
                  }
                  int gjN = genJetPhi->size();


            hEvent_Pass->Fill(1);


            
            
            
            
            
            //ENTERING JET LOOP

            //in this first loop I choose the jetA and count the trks of A 
            for(int kjet=0; kjet < genJetPt->size(); kjet++){

                int ijet = kjet;
                int Gjet = kjet;
                long int NNtrk1 = (genDau_pt->at(ijet)).size();

                TVector3 JetA;
                JetA.SetPtEtaPhi((*genJetPt)[ijet],(*genJetEta)[ijet],(*genJetPhi)[ijet]);
                TLorentzVector JetA_4 (JetA, JetA.Mag());
                
                if( fabs(JetA.Eta()) > jetEtaCut ) continue;
                if( JetA.Perp() < jetPtCut_Jet   ) continue;


                 //    count the trks of A  
                int n_G_ChargeMult_count  = 0;
                int n_G_ChargeMult_count1 = 0;
                
                
                for(int  G_trk=0; G_trk < NNtrk1; G_trk++ ){
                    if((*genDau_chg)[Gjet][G_trk] == 0) continue;
                    if(fabs((*genDau_pt)[Gjet][G_trk])  < 0.3)     continue;
                    if(fabs((*genDau_eta)[Gjet][G_trk]) > 2.4)     continue;
                    n_G_ChargeMult_count1 += 1;
                }
                // if (n_G_ChargeMult_count1<30) continue;
                
                
                for(int jjet=ijet+1; (jjet< genJetPt->size()); jjet++){


                    TVector3 JetB;
                    JetB.SetPtEtaPhi((*genJetPt)[jjet],(*genJetEta)[jjet],(*genJetPhi)[jjet]);
                    TLorentzVector JetB_4 (JetB, JetB.Mag());

                    if( fabs(JetB.Eta()) > jetEtaCut ) continue;
                    if( JetB.Perp() < jetPtCut_Jet-200 ) continue;
                    // if( JetB.Perp() < JetA.Perp()*0.95 ) continue;
              
                    TLorentzVector Boost_to_CM = JetA_4 + JetB_4;
                    TLorentzVector JetAA_4 = BeamBoost(Boost_to_CM,JetA_4);
                    TLorentzVector JetBB_4 = BeamBoost(Boost_to_CM,JetB_4);

                    TVector3 JetAA = JetAA_4.Vect();
                    TVector3 JetBB = JetBB_4.Vect();

                    // TVector3 JetAB = BeamBoost(JetA.Perp(),JetA.Eta(),JetA.Phi(),JetB.Perp(),JetB.Eta(),JetB.Phi());
                    double deltaJetEta = JetAA.Eta() + JetBB.Eta();
                    double deltaJetPhi =  fabs(JetA.Phi()-JetB.Phi());
                    hdeltaJetPhi -> Fill(deltaJetPhi);
                    hdeltaJetEta -> Fill(deltaJetEta);

                    if (fabs(M_PI-deltaJetPhi) > 0.1) continue;
                    // if (fabs(deltaJetEta)>0.15) continue;
                    long int NNtrk2 = (genDau_pt->at(jjet)).size();
                    // hdeltaR -> Fill(deltaR);

                   
                    
                    // Calculate the trks in jetB
                    int n_G_ChargeMult_count2 = 0;
                    for(int G_trk2=0; G_trk2 < NNtrk2; G_trk2++ ){
                        if((*genDau_chg)[jjet][G_trk2] == 0) continue;
                        if(fabs((*genDau_pt)[jjet][G_trk2])  < 0.3)     continue;
                        if(fabs((*genDau_eta)[jjet][G_trk2]) > 2.4)     continue;
                        n_G_ChargeMult_count2 += 1;
                    }

                    

                    // n_G_ChargeMult_count = n_G_ChargeMult_count1 + n_G_ChargeMult_count2 ;
                    // n_G_ChargeMult_count = ((1+floor(sqrt(1+(4*2*n_G_ChargeMult_count1*n_G_ChargeMult_count2))))/2) ;
                    
                    n_G_ChargeMult_count = n_G_ChargeMult_count1;
                    hBinDist_gen_single            ->Fill(n_G_ChargeMult_count1);


                    int tkBool[trackbin] = {0};
                    int Ntrig[trackbin][ptbin] = {0};
                    int NtrigM[trackbin][ptbin] = {0};
                    int NtrigP[trackbin][ptbin] = {0};
                    int A_ptBool[NNtrk1][ptbin] = {0};    //
                    int T_ptBool[NNtrk2][ptbin]     = {0};// This is for the AB

                    for(int i = 0; i < trackbin; i++){
                    //if((*chargedMultiplicity)[indicesR[kjet]] >= trackbinbounds[i] && (*chargedMultiplicity)[indicesR[kjet]] < trackbinboundsUpper[i]){
                        if(n_G_ChargeMult_count >= trackbinbounds[i] && n_G_ChargeMult_count < trackbinboundsUpper[i]){
                            tkBool[i] = 1;
                            // hJet_Pass           ->Fill(i);
                            hBinDist_gen[i]         ->Fill(n_G_ChargeMult_count);
                        }
                    }



                    std::complex<double> Q_all2 (0, 0);
                    std::complex<double> Q_all4 (0, 0);

                    // calculate the A_ptbool pile up Ntrig in jetA first,
                    //and then we do this in jetB so that we can get the complete Ntrig
                    for(int  A_trk=0; A_trk < NNtrk1; A_trk++ ){
                    

                        gRandom->SetSeed(0);
                        double phi;
                        phi = hPhiDraw->GetRandom();


                        
                        std::complex<double> Q_part2 (TMath::Cos(n_harm*phi) , TMath::Sin(n_harm*phi));
                        Q_all2 = Q_all2 + Q_part2;

                        std::complex<double> Q_part4 (TMath::Cos(2.0*n_harm*phi) , TMath::Sin(2.0*n_harm*phi));
                        Q_all4 = Q_all4 + Q_part4;

                        // hPhiDrawA->Fill(phi);



                            
                    }


                    int M = NNtrk1;

                    double Q_all_abs = std::abs(Q_all2);
                    double Q_all_sqr = Q_all_abs * Q_all_abs;
                    double particle_two = (Q_all_sqr - M) / (M*(M-1));
                    int weight_two = (M*(M-1));


                    for(int i = 0; i < trackbin; i++){
                        if(tkBool[i] == 1){

                        jet_avg_numerator_two[i] = jet_avg_numerator_two[i] + ((weight_two)*(particle_two)); 
                        jet_avg_denominat_two[i] = jet_avg_denominat_two[i] + (weight_two);
                        }
                    }

                    double Q_all_fourth = Q_all_sqr * Q_all_sqr;            
                    double Q_all_abs_2n = std::abs(Q_all4);   
                    double Q_all_sqr_2n = Q_all_abs_2n * Q_all_abs_2n;      
                    double Q_all_re_threetrm = std::real(
                            (2.0 * Q_all4) *
                            std::conj( Q_all2 ) *
                            std::conj( Q_all2 )
                            );                                              
                    double particle_four =
                        (   (Q_all_fourth + Q_all_sqr_2n - 2.0 * Q_all_re_threetrm)
                            / ( M*(M-1)*(M-2)*(M-3) )
                        )
                        -
                        (   2.0* ( (2.0*(M-2)*Q_all_sqr) - (M*(M-3)) )
                            / ( M*(M-1)*(M-2)*(M-3) )
                        );
                    int weight_four = M*(M-1)*(M-2)*(M-3);

                    for(int i = 0; i < trackbin; i++){
                        if(tkBool[i] == 1){

                        jet_avg_numerator_four[i] = jet_avg_numerator_four[i] + ((weight_four)*(particle_four));
                        jet_avg_denominat_four[i] = jet_avg_denominat_four[i] + (weight_four);
                        }
                    }


                }//jjet
             
            }//kjet
                    }//Event
                    fFile->Close();
                    }//File

    double Rand_c22[trackbin] = {0};
    double Rand_c24[trackbin] = {0};
    double Rand_v22[trackbin] = {0};
    double Rand_v24[trackbin] = {0};

    for(int i = 1; i < trackbin; i++){
        Rand_mult_bin_avg_two[i]  =Rand_jet_avg_numerator_two[i]/Rand_jet_avg_denominat_two[i];
        Rand_mult_bin_avg_four[i] =Rand_jet_avg_numerator_four[i]/Rand_jet_avg_denominat_four[i];

        Rand_c22[i] = Rand_mult_bin_avg_two[i];
        Rand_c24[i] = Rand_mult_bin_avg_four[i] - (2* Rand_mult_bin_avg_two[i] * Rand_mult_bin_avg_two[i]);
        Rand_v22[i] = std::real(std::pow(std::complex<double>(Rand_c22[i]),.50));
        Rand_v24[i] = std::real(std::pow(std::complex<double>(-Rand_c24[i]),.25));

//std::real(std::pow(std::complex<double>(-Rand_c24[i]),.25))

        hRand_c22->Fill(mult_bins[i]+.5, Rand_c22[i]);
        hRand_c24->Fill(mult_bins[i]+.5, Rand_c24[i]);
        hRand_v22->Fill(mult_bins[i]+.5, Rand_v22[i]);
        hRand_v24->Fill(mult_bins[i]+.5, Rand_v24[i]);
    }
    TFile* fS_2_temp = new TFile(Form("pythia_batch_output/root_out_2/rand_dijob_%s.root",subList.c_str()), "recreate");
    hRand_c22->Write();
    hRand_c24->Write();
    hRand_v22->Write();
    hRand_v24->Write();
    fS_2_temp->Close();

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