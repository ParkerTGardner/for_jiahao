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

    TH2D* hPairs        = new TH2D("hPairs","hPairs",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);

    TH2D* hNtrig        = new TH2D("hNtrig","hNtrig",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);


    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", trackbin,bin0,trackbin);
    TH1D* hJet_Pass     = new TH1D("hJet_Pass"  ,"hJet_Pass"  , trackbin,bin0,trackbin);
    TH1D* hJetPt = new TH1D("hJetPt"  ,"hJetPt" ,i150,i100,i3500);
    TH1D* hBinDist_gen_single = new TH1D("hBinDist_gen_single","hBinDist_gen_single",bin360,bin0, bin120);
    TH1D* hBinDist_reco_single = new TH1D("hBinDist_reco_single","hBinDist_reco_single",bin360,bin0,bin120);
    TH1D* hdeltaR = new TH1D("hdeltaR", "hdeltaR", 100,2.5,4);
    TH1D* hdeltaJetPhi = new TH1D("hdeltaJetPhi", "hdeltaJetPhi",50,0,5/4*M_PI);
    TH1D* hdeltaJetTheta = new TH1D("hdeltaJetTheta", "hdeltaJetTheta",50,-M_PI,M_PI);
    TH1D* hdeltaJetEta = new TH1D("hdeltaJetEta", "hdeltaJetEta",50,-4,4);
    TH1D* hT_jet_dau_eta = new TH1D("hT_jet_dau_eta","hT_jet_dau_eta",100,-5,5);
    TH1D* hT_jet_dau_phi = new TH1D("hT_jet_dau_phi","hT_jet_dau_phi",100,-M_PI,M_PI);
    TH1D* hJet_Eta_ave_cut1 = new TH1D("hJet_Eta_ave_cut1","hJet_Eta_ave_cut1",100,-2,2);
    TH1D* hJet_Ptw_Eta_ave_cut1 = new TH1D("hJet_Ptw_Eta_ave_cut1","hJet_Ptw_Eta_ave_cut1", 100,-2,2);
    TH1D* hJet_Eta_ave_cutR = new TH1D("hJet_Eta_ave_cutR","hJet_Eta_ave_cutR",100,-2,2);
    TH1D* hJet_Ptw_Eta_ave_cutR = new TH1D("hJet_Ptw_Eta_ave_cutR","hJet_Ptw_Eta_ave_cutR", 100,-2,2);
    TH1D* hJet_Eta_ave_nocut = new TH1D("hJet_Eta_ave_nocut","hJet_Eta_ave_nocut",100,-2,2);
    TH1D* hJet_Ptw_Eta_ave_nocut = new TH1D("hJet_Ptw_Eta_ave_nocut","hJet_Ptw_Eta_ave_nocut", 100,-2,2);
    TH1D* hT_jet_dau_pt  = new TH1D("hT_jet_dau_pt","hT_jet_dau_pt",100,0,10);
    // TH1D* htest_t_eta = new TH1D("htest_t_eta","htest_t_eta",100,-5,5);
    // TH1D* htest_t_phi = new TH1D("htest_t_phi","htest_t_phi",100,-M_PI,M_PI);
    // TH1D* htest_t_eta0 = new TH1D("htest_t_eta0","htest_t_eta0",100,-5,5);
    // TH1D* htest_t_phi0 = new TH1D("htest_t_phi0","htest_t_phi0",100,-M_PI,M_PI);
    // TH1D* htest_a_eta = new TH1D("htest_a_eta","htest_a_eta",100,-5,5);
    // TH1D* htest_a_phi = new TH1D("htest_a_phi","htest_a_phi",100,-M_PI,M_PI);
    // TH1D* htest_a_eta0 = new TH1D("htest_a_eta0","htest_a_eta0",100,-5,5);
    // TH1D* htest_a_phi0 = new TH1D("htest_a_phi0","htest_a_phi0",100,-M_PI,M_PI);
    // TH1D* htest_t_deltaeta = new TH1D("htest_t_deltaeta","htest_t_deltaeta",100,-5,5);
    // TH1D* htest_t_deltaphi = new TH1D("htest_t_deltaphi","htest_t_deltaphi",100,-M_PI,M_PI);
    // TH1D* htest_t_deltaeta0 = new TH1D("htest_t_deltaeta0","htest_t_deltaeta0",100,-5,5);
    // TH1D* htest_t_deltaphi0 = new TH1D("htest_t_deltaphi0","htest_t_deltaphi0",100,-M_PI,M_PI);

    
    //2D Corr histograms

    // TH3D* hEPDraw[trackbin][ptbin][PUbin];
    TH2D* hSignalShifted[trackbin][ptbin][PUbin];
    TH2D* hBckrndShifted[trackbin][ptbin][PUbin];
    TH2D* hEPDrawA[trackbin][ptbin][PUbin];
    TH2D* hEPDrawT[trackbin][ptbin][PUbin];
    // TH2D* hPPT[trackbin][ptbin][PUbin];
    TH2D* hJJT                  = new TH2D(Form("hJJT") ,Form("hJJT") , 100, 0, 1);
    TH2D* hJJT_Cut              = new TH2D(Form("hJJT_Cut") ,Form("hJJT_Cut") , 100, 0, 1);

  

    


    TH1D* hBinDist_gen[trackbin];
    TH1D* hBinDist_reco[trackbin];
    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        hBinDist_gen[wtrk-1]    = new TH1D(Form("hBinDist_gen_%d",wtrk),Form("hBinDist_gen_%d",wtrk), bin360, bin0, bin120);
        hBinDist_reco[wtrk-1]   = new TH1D(Form("hBinDist_reco_%d",wtrk),Form("hBinewnDist_reco_%d",wtrk), bin360, bin0, bin120);
        for(int wppt = 1; wppt<ptbin+1; wppt++){
            for(int wpPU = 1; wpPU<PUbin+1; wpPU++){

                //these are the main histograms
                // hMomBckrndShifted[wtrk-1][wppt-1][wpPU-1]   = new TH1D(Form("hjtw_MomBckrndS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form("hjtw_MomBckrndS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) , 30, 0, 1.2);
                // hMomSignalShifted[wtrk-1][wppt-1][wpPU-1]   = new TH1D(Form("hjtw_MomSignalS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form("hjtw_MomSignalS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) , 30, 0, 1.2);
                hBckrndShifted[wtrk-1][wppt-1][wpPU-1]      = new TH2D(Form("hBckrndS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form("hBckrndS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                hSignalShifted[wtrk-1][wppt-1][wpPU-1]      = new TH2D(Form("hSignalS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form("hSignalS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                // hEPDraw[wtrk-1][wppt-1][wpPU-1]             = new TH3D(Form("hjtw_EPDraw_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form( "hjtw_EPDraw_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
                hEPDrawA[wtrk-1][wppt-1][wpPU-1]             = new TH2D(Form("hEPDrawA_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form( "hEPDrawA_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
                hEPDrawT[wtrk-1][wppt-1][wpPU-1]             = new TH2D(Form("hEPDrawT_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form( "hEPDrawT_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) , 2*EPD_xb   , -EPD_xhi, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
                // hPPT[wtrk-1][wppt-1][wpPU-1]                  = new TH2D(Form("hPPT_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form( "hPPT_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) , 100, 0, 1);
                
                // h_jet_dau[wtrk-1][wppt-1][wpPU-1]             = new TH2D(Form("hjtw_BckrndS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form("hjtw_BckrndS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                // TH2D* h_jet_dau_A[trackbin][ptbin][PUbin];
                // TH2D* h_jet_dau_T[trackbin][ptbin][PUbin];

            }
        }
    }


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

                //defining how many particles are in the jet, set observable based on jet multiplicity
                //for example, what is the shape of jet(5,10,100 particles in jet)
                

                if( fabs((*genJetEta)[ijet]) > jetEtaCut ) continue;
                if( (*genJetPt)[ijet] < jetPtCut_Jet   ) continue;

                //    count the trks of A  
                int n_G_ChargeMult_count = 0;
                int n_G_ChargeMult_count1 =0;
                
                
                for(int  G_trk=0; G_trk < NNtrk1; G_trk++ ){
                    if((*genDau_chg)[Gjet][G_trk] == 0) continue;
                    if(fabs((*genDau_pt)[Gjet][G_trk])  < 0.3)     continue;
                    if(fabs((*genDau_eta)[Gjet][G_trk]) > 2.4)     continue;
                    n_G_ChargeMult_count1 += 1;
                }

                // Enter the jetB
                // Actually most of our calculation is in this loop
                for(int jjet=ijet+1; (jjet< genJetPt->size()); jjet++){

                    // hJet_Eta_ave_nocut->Fill(((*genJetEta)[ijet]+(*genJetEta)[jjet])/2.0);
                    // hJet_Ptw_Eta_ave_nocut->Fill(((*genJetEta)[ijet]*(*genJetPt)[ijet]+(*genJetEta)[jjet]*(*genJetPt)[jjet])/((*genJetPt)[ijet]+(*genJetPt)[jjet]));

                    // choose the jetB                    
                    // if (jjet == ijet) continue;
                    hJJT->Fill((*genJetPt)[ijet]/(*genJetPt)[jjet],(*genJetPt)[ijet]);
                    if( fabs((*genJetEta)[jjet]) > jetEtaCut ) continue;
                    if( (*genJetPt)[jjet] < jetPtCut_Jet-200   ) continue;    
                    hJJT_Cut->Fill((*genJetPt)[ijet]/(*genJetPt)[jjet],(*genJetPt)[ijet]);

                    hJet_Eta_ave_cut1->Fill(((*genJetEta)[jjet]-(*genJetEta)[ijet])/2.0);
                    hJet_Ptw_Eta_ave_cut1->Fill((((*genJetEta)[jjet]-(*genJetEta)[ijet])*(*genJetPt)[jjet])/((*genJetPt)[ijet]+(*genJetPt)[jjet]));  

                    // double deltaJetTheta = thetaWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genJetPt)[jjet], -(double)(*genJetEta)[jjet], (double)(*genJetPhi)[jjet]+M_PI);
                    // double deltaJetEta0   = etaWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genJetPt)[jjet], (double)(*genJetEta)[jjet], (double)(*genJetPhi)[jjet]);
                    // double deltaJetPhi    = phiWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genJetPt)[jjet], (double)(*genJetEta)[jjet], (double)(*genJetPhi)[jjet]);
                    // double deltaJetTheta0 = 2*TMath::ATan(TMath::Exp(-deltaJetEta0));
                    // double deltaJetTheta  = M_PI/2 - deltaJetTheta0;
                    // double deltaJetEta    = -TMath::Log( TMath::Tan( deltaJetTheta/2.0));
                    //if we take symmectric jet and deltaJetPhi<0, it means the jet just situates nearside
                    // deltaJetTheta = M_PI - deltaJetTheta;
                    
                    // double deltaJetR1   = sqrt( pow((deltaJetPhi),2) + pow(deltaJetEta,2));
                    // double deltaJetR2   = sqrt( pow((M_PI-deltaJetPhi),2) + pow(deltaJetEta,2));
                    double deltaJetEta =  -(double)(*genJetEta)[jjet]+(*genJetEta)[ijet];
                    // double deltaJetTheta = fabs(2*TMath::ATan(TMath::Exp(-(double)(*genJetEta)[ijet])) - 2*TMath::ATan(TMath::Exp(-(double)(*genJetEta)[jjet])));
                    double deltaJetPhi =  fabs((double)(*genJetPhi)[ijet]-(double)(*genJetPhi)[jjet]);
                    double deltaR = sqrt(pow(M_PI-deltaJetPhi,2)+pow(deltaJetEta,2));
                    hdeltaR -> Fill(deltaR);
                    hdeltaJetEta -> Fill(deltaJetEta);
                    hdeltaJetPhi -> Fill(deltaJetPhi);
                    // hdeltaJetTheta -> Fill(deltaJetTheta);


                    // We only consider the back-to-back system, so we choose deltaeta<0
                    // if (deltaJetEta < 0) continue;
                    if (fabs(M_PI-deltaJetPhi) > 0.1) continue;
                    if (deltaR > 0.8 )  continue;
                    // if (fabs(deltaJetEta)>1.6) continue;

                    hJet_Eta_ave_cutR->Fill(((*genJetEta)[jjet]-(*genJetEta)[ijet])/2.0);
                    hJet_Ptw_Eta_ave_cutR->Fill((((*genJetEta)[jjet]-(*genJetEta)[ijet])*(*genJetPt)[jjet])/((*genJetPt)[ijet]+(*genJetPt)[jjet]));  

                    long int NNtrk2 = (genDau_pt->at(jjet)).size();


                    // Calculate the trks in jetB
                    int n_G_ChargeMult_count2 = 0;
                    for(int G_trk2=0; G_trk2 < NNtrk2; G_trk2++ ){
                        if((*genDau_chg)[jjet][G_trk2] == 0) continue;
                        if(fabs((*genDau_pt)[jjet][G_trk2])  < 0.3)     continue;
                        if(fabs((*genDau_eta)[jjet][G_trk2]) > 2.4)     continue;
                        n_G_ChargeMult_count2 += 1;

                    }
                    // Here we define the mult of AB as max(A,B)
                    // n_G_ChargeMult_count = max(n_G_ChargeMult_count1 , n_G_ChargeMult_count2) ;
                    n_G_ChargeMult_count = n_G_ChargeMult_count1 + n_G_ChargeMult_count2 ;




                    hBinDist_gen_single            ->Fill(n_G_ChargeMult_count);


                    //some useful bools 
                    // We have to do this here because we need to get correct mult first
                    // and then get tktool, Ntrig

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
                            hJet_Pass           ->Fill(i);
                            hBinDist_gen[i]         ->Fill(n_G_ChargeMult_count);     
                        }
                    }
                

                    
                    
                    

                    // calculate the A_ptbool pile up Ntrig in jetA first,
                    //and then we do this in jetB so that we can get the complete Ntrig
                    for(int  A_trk=0; A_trk < NNtrk1; A_trk++ ){
                        
                    
                        if((*genDau_chg)[ijet][A_trk] == 0) continue;
                        if(fabs((*genDau_eta)[ijet][A_trk]) > 2.4) continue;
                        if(fabs((*genDau_pt)[ijet][A_trk])  < 0.3)     continue;

                            //     daughter pt with respect to the jetA axis                 pt With Respect To JetA
                        double jet_dau_pt    =  ptWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);

                        if(jet_dau_pt >3.0) continue;// why we drop this

                        // find A_trk in i ptbin
                        for(int i = 0; i < ptbin; i++){
                            if(jet_dau_pt >= ptbinbounds_lo[i] && jet_dau_pt < ptbinbounds_hi[i]){
                                A_ptBool[A_trk][i] = 1;
                            }
                        }

                        // for ith trkbin and jth ptbin, pile up the Ntrig of jetA

                        for(int i = 0; i < trackbin; i++){
                            for(int j = 0; j < ptbin; j++){
                                if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                    Ntrig[i][j] += 1;
                                    if((*genDau_chg)[ijet][A_trk] > 0){
                                        NtrigP[i][j] += 1;
                                    }
                                    if((*genDau_chg)[ijet][A_trk] < 0){
                                        NtrigM[i][j] += 1;
                                    }
                                }
                            }
                        }    
                    }

                    
                    // Multiplicity = max(jetA,jetB) ;


                    
                
                
                    // Get the T_ptBool
                    
                    for(long int T_trk=0; T_trk< NNtrk2; T_trk++ ){

                        if((*genDau_chg)[jjet][T_trk] == 0) continue;
                        if(fabs((*genDau_eta)[jjet][T_trk]) > 2.4) continue;
                        if(fabs((*genDau_pt)[jjet][T_trk])  < 0.3)      continue;

                        double T_jet_dau_pt  = 0;
                        double T_jet_dau_eta = 0;
                        double T_jet_dau_phi = 0;

                        // pt, eta phi wrt Jet A, here we consider eta as -eta (reflect)
                        
                        double T_jet_dau_pt0    =  ptWRTJet((double)(*genJetPt)[jjet], (double)(*genJetEta)[jjet], (double)(*genJetPhi)[jjet], (double)(*genDau_pt)[jjet][T_trk], (double)(*genDau_eta)[jjet][T_trk], (double)(*genDau_phi)[jjet][T_trk]);  
                        if(T_jet_dau_pt0 >3.0) continue;

                        // double T_dau_Theta0 = 2*TMath::ATan(TMath::Exp(-(double)(*genDau_eta)[jjet][T_trk]))+deltaJetTheta;
                        // double T_dau_eta = -TMath::Log( TMath::Tan( T_dau_Theta0/2.0));

                        //     daughter eta with respect to the jetA axis                 eta With Respect To JetA 
                        T_jet_dau_eta   = etaWRTJet((double)(*genJetPt)[ijet], 0, (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[jjet][T_trk], (double)(*genDau_eta)[jjet][T_trk]-(double)(*genJetEta)[ijet], (double)(*genDau_phi)[jjet][T_trk]);
                        //     daughter phi with respect to the jetA axis                 phi With Respect To JetA
                        T_jet_dau_phi   = phiWRTJet((double)(*genJetPt)[ijet], 0, (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[jjet][T_trk], (double)(*genDau_eta)[jjet][T_trk]-(double)(*genJetEta)[ijet], (double)(*genDau_phi)[jjet][T_trk]); 
                        
                        T_jet_dau_pt    =  ptWRTJet((double)(*genJetPt)[ijet], 0, (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[jjet][T_trk], (double)(*genDau_eta)[jjet][T_trk]-(double)(*genJetEta)[ijet], (double)(*genDau_phi)[jjet][T_trk]);

                        hT_jet_dau_eta->Fill(T_jet_dau_eta);
                        hT_jet_dau_phi->Fill(T_jet_dau_phi);
                        hT_jet_dau_pt ->Fill(T_jet_dau_pt);
                        
                        
                        
                        for(int i = 0; i < ptbin; i++){
                            if(T_jet_dau_pt0 >= ptbinbounds_lo[i] && T_jet_dau_pt0 < ptbinbounds_hi[i]){
                                T_ptBool[T_trk][i] = 1;
                            }
                        }

                        // The same as for jetA
                        
                        for(int i = 0; i < trackbin; i++){
                            for(int j = 0; j < ptbin; j++){
                                if(tkBool[i] + T_ptBool[T_trk][j] == 2){
                                    
                                    Ntrig[i][j] += 1;
                                    if((*genDau_chg)[jjet][T_trk] > 0){
                                        NtrigP[i][j] += 1;
                                    }
                                    if((*genDau_chg)[jjet][T_trk] < 0){
                                        NtrigM[i][j] += 1;
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


                    //Main loop for corr func

                    // first particle loop, here first particle is only from jetA, and the second just from jetB
                    // Actually if we go over all the jet and get jet pairs, and we will get jetAB and jetBA
                    // Here we only calculate corr(dau_A,dau_A) and corr(dau_A,dau_B)
                    // when kjet == jjet and jjet == kjet must suitable 
                    // and we can get corr(dau_B, dau_A) and corr(dau_B, dau_B) this time
                    // for example, here mulA=mulB=5, if all daus suitable we will get 2C10=45 pairs
                    // and for the first loop we wiil get 2C5+ 5*5 pairs, 
                    // and we should weight the corr(dau_A,dau_B) 1/2 to avoid double-counting


                    for(int  A_trk=0; A_trk < NNtrk1; A_trk++ ){
                        
                        if((*genDau_chg)[ijet][A_trk] == 0) continue;
                        if(fabs((*genDau_eta)[ijet][A_trk]) > 2.4) continue;
                        if(fabs((*genDau_pt)[ijet][A_trk])  < 0.3)     continue;

                            //     daughter pt with respect to the jet axis                 pt With Respect To Jet 
                        double jet_dau_pt    =  ptWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);


                            //if(jet_dau_pt >3.0) continue;

                        if(jet_dau_pt >3.0) continue;// why we drop this

                        //     daughter eta with respect to the jet axis                 eta With Respect To Jet 
                        double jet_dau_eta   = etaWRTJet((double)(*genJetPt)[ijet], 0, (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk]-(double)(*genJetEta)[ijet], (double)(*genDau_phi)[ijet][A_trk]);
                        //     daughter phi with respect to the jet axis                 phi With Respect To Jet 
                        double jet_dau_phi   = phiWRTJet((double)(*genJetPt)[ijet], 0, (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk]-(double)(*genJetEta)[ijet], (double)(*genDau_phi)[ijet][A_trk]);

                        double jet_dau_theta = 2*ATan(Exp(-(jet_dau_eta)));


                        // the background should be just from this single jet loop
                        // Actully because of the jjet loop, it will be overcounting the const times( the number of jets)
                        // And it has no effect at all

                        // for(int i = 0; i < trackbin; i++){
                        //     for(int j = 0; j < ptbin; j++){
                        //         if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                        //             int k_PU=0;

                        //             hEPDraw[i][j][k_PU]->Fill(jet_dau_eta, jet_dau_phi, jet_dau_pt, (float)((float)(1.0)/Ntrig[i][j]));;
                        //         }
                        //     }
                        // }

                        for(int i = 0; i < trackbin; i++){
                            for(int j = 0; j < ptbin; j++){
                                if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                    int k_PU=0;

                                    hEPDrawA[i][j][k_PU]->Fill(jet_dau_eta, jet_dau_phi, ((double)(1.0)/Ntrig[i][j]));
                                }
                            }
                        }



                        

                        //A_trk is the first track from the first loop
                        //T_trk is the second loop


                        // This is for corr(dau_A, dau_B)
                        for(long int T_trk = 0; T_trk < NNtrk2; T_trk++ ){

                            if((*genDau_chg)[jjet][T_trk] == 0) continue;
                            if(fabs((*genDau_eta)[jjet][T_trk]) > 2.4) continue;
                            if(fabs((*genDau_pt)[jjet][T_trk])  < 0.3)      continue;



                            double T_jet_dau_pt  = 0;
                            double T_jet_dau_eta = 0;
                            double T_jet_dau_phi = 0;


                            
                            T_jet_dau_pt    =  ptWRTJet((double)(*genJetPt)[jjet], (double)(*genJetEta)[jjet], (double)(*genJetPhi)[jjet], (double)(*genDau_pt)[jjet][T_trk], (double)(*genDau_eta)[jjet][T_trk], (double)(*genDau_phi)[jjet][T_trk]);  
                            if(T_jet_dau_pt >3.0) continue;

                            // double T_dau_Theta0 = 2*TMath::ATan(TMath::Exp(-(double)(*genDau_eta)[jjet][T_trk]))+deltaJetTheta;
                            // double T_dau_eta = -TMath::Log( TMath::Tan( T_dau_Theta0/2.0));

                            //     daughter eta with respect to the jetA axis                 eta With Respect To JetA 
                            T_jet_dau_eta   = etaWRTJet((double)(*genJetPt)[ijet], 0, (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[jjet][T_trk], (double)(*genDau_eta)[jjet][T_trk]-(double)(*genJetEta)[ijet], (double)(*genDau_phi)[jjet][T_trk]);
                            //     daughter phi with respect to the jetA axis                 phi With Respect To JetA
                            T_jet_dau_phi   = phiWRTJet((double)(*genJetPt)[ijet], 0, (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[jjet][T_trk], (double)(*genDau_eta)[jjet][T_trk]-(double)(*genJetEta)[ijet], (double)(*genDau_phi)[jjet][T_trk]); 

                            T_jet_dau_pt    =  ptWRTJet((double)(*genJetPt)[ijet], 0, (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[jjet][T_trk], (double)(*genDau_eta)[jjet][T_trk]-(double)(*genJetEta)[ijet], (double)(*genDau_phi)[jjet][T_trk]);

                            double T_jet_dau_eta0   = etaWRTJet((double)(*genJetPt)[jjet], (double)(*genJetEta)[jjet], (double)(*genJetPhi)[jjet], (double)(*genDau_pt)[jjet][T_trk], (double)(*genDau_eta)[jjet][T_trk], (double)(*genDau_phi)[jjet][T_trk]);
                            
                            if(T_jet_dau_eta0 > track_eta_lim) continue;

                            // for(int i = 0; i < trackbin; i++){
                            //     for(int j = 0; j < ptbin; j++){
                            //         if(tkBool[i] + T_ptBool[T_trk][j] == 2){
                            //             int k_PU=0;

                            //             hEPDraw[i][j][k_PU]->Fill(T_jet_dau_eta, T_jet_dau_phi , T_jet_dau_pt, (float)((float)(1.0)/Ntrig[i][j]));
                            //         }
                            //     }
                            // }

                            for(int i = 0; i < trackbin; i++){
                                for(int j = 0; j < ptbin; j++){
                                    if(tkBool[i] + T_ptBool[T_trk][j] == 2){
                                        int k_PU=0;

                                        hEPDrawT[i][j][k_PU]->Fill(T_jet_dau_eta, T_jet_dau_phi ,  ((double)(1.0)/Ntrig[i][j]));
                                    }
                                }
                            }
                            

                            //correlation function
                                             //A_trk(dau_A)  T_trk(dau_B)
                            double deltaEta = (jet_dau_eta - T_jet_dau_eta);
                                                                        //A_trk        T_trk
                            double deltaPhi = (TMath::ACos(TMath::Cos(jet_dau_phi - T_jet_dau_phi)));

                            // htest_t_deltaeta0->Fill(deltaEta);
                            // htest_t_deltaphi0->Fill(deltaPhi);
                            

                            // double deltaJt  = fabs(jet_dau_pt - T_jet_dau_pt);
                            
                            for(        int i = 0; i < trackbin; i++){
                                for(    int j = 0; j < ptbin;    j++){ 
            


                                    if(tkBool[i] + A_ptBool[A_trk][j] + T_ptBool[T_trk][j] == 3){
                                            hPairs->Fill(i,j);
                                            int k_PU=0;
                                            hSignalShifted[i][j][k_PU]->Fill(deltaEta, deltaPhi,                 ((double)(1.0)/Ntrig[i][j]));
                                            hSignalShifted[i][j][k_PU]->Fill(-deltaEta, deltaPhi,                ((double)(1.0)/Ntrig[i][j]));
                                            hSignalShifted[i][j][k_PU]->Fill(deltaEta, -deltaPhi,                ((double)(1.0)/Ntrig[i][j]));
                                            hSignalShifted[i][j][k_PU]->Fill(-deltaEta, -deltaPhi,               ((double)(1.0)/Ntrig[i][j]));
                                            hSignalShifted[i][j][k_PU]->Fill( deltaEta,2*TMath::Pi() - deltaPhi, ((double)(1.0)/Ntrig[i][j]));
                                            hSignalShifted[i][j][k_PU]->Fill(-deltaEta,2*TMath::Pi() - deltaPhi, ((double)(1.0)/Ntrig[i][j]));
                                            // hMomSignalShifted[i][j][k_PU]->Fill(deltaJt,                         1/(Ntrig[i][j]));




                                        //}}}
                                        // This is the mixed charge signal. Each duaghter will serve as a trigger so regular Ntrig suffices.
                                        // EPdraw will constitute sampling from EPD_P and EPD_M for the pseudo particles. 
                                        // So that I take random phi and eta from P and random Phi and Eta from M and make a opposite sign difference.
                                    }

                        
                                }
                            }

                
                        }//T_trk;  AB


                    }





                }//jjet
             
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
                                  for(int wpPU = 1; wpPU < PUbin+1; wpPU++){
                                      std::cout << wpPU << "/" << PUbin << std::endl;
                                      //Nent is the number of pairs in the signal which we will try to 10x
                                      //Xent is the number of pseudoparticles requried such that when we build the pairs nCp = Xent CHOOSE 2 will give 
                                      //us 10 times as many pairs as we have in the signal histogrm.

                                      long int NENT =  hPairs->GetBinContent(wtrk, wppt);
                                      long int XENT =  ((1+floor(sqrt(1+(4*2*backMult*NENT))))/2) ;
                                      float A_ETA[XENT] = {0};
                                      float A_PHI[XENT] = {0};
                                      float T_ETA[XENT] = {0};
                                      float T_PHI[XENT] = {0};
				                    //   float A_Jt[XENT]  = {0};

                                      for(int x = 0; x<XENT; x++){
                                          gRandom->SetSeed(0);
                                          double WEtaA, WPhiA;//making the pseudoparticles
                                          double WEtaT, WPhiT;
                                          hEPDrawA[wtrk-1][wppt-1][wpPU-1]->GetRandom2(WEtaA, WPhiA);
                                          hEPDrawT[wtrk-1][wppt-1][wpPU-1]->GetRandom2(WEtaT, WPhiT);
                                          A_ETA[x] = WEtaA;
                                          A_PHI[x] = WPhiA;
                                          T_ETA[x] = WEtaT;
                                          T_PHI[x] = WPhiT;
                                        //   A_Jt[x]  = WJt1;
                                      }
                                      for(long int i = 0; i < (XENT-1); i++){
                                          for(long int j = (i+1); j < XENT; j++){

                                              double WdeltaEta = (A_ETA[i]-T_ETA[j]);
                                              double WdeltaPhi = (TMath::ACos(TMath::Cos(A_PHI[i]-T_PHI[j])));
                                            //   double WdeltaJt  = fabs(A_Jt[i]-A_Jt[j]);

                                              hBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(WdeltaEta, WdeltaPhi, 1);//./XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(-WdeltaEta, WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(WdeltaEta, -WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(-WdeltaEta, -WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(WdeltaEta, 2*TMath::Pi() - WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(-WdeltaEta,2*TMath::Pi() - WdeltaPhi, 1);//../XENT);
                                            //   hMomBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(WdeltaJt, 1);

                                          }
                                      }
                                  }
                              }
                          }
                      //}}}




                    string subList = fList.substr(fList.size() - 3);

                    TFile* fS_tempA = new TFile(Form("pythia_batch_output/root_out_3/dijob_%s.root",subList.c_str()), "recreate");
                    for(int wtrk =1; wtrk <trackbin+1; wtrk++){
                        hBinDist_gen[wtrk-1]         ->Write();
                        for(int wppt =1; wppt <ptbin+1; wppt++){
                            for(int wpPU =1; wpPU<PUbin+1; wpPU++){

                                hSignalShifted             [wtrk-1][wppt-1][wpPU-1]->Write(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU    ));
                                hBckrndShifted             [wtrk-1][wppt-1][wpPU-1]->Write(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU    ));
                                // hEPDraw                    [wtrk-1][wppt-1][wpPU-1]->Write(Form("hEPD_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU     ));
                                hEPDrawA                   [wtrk-1][wppt-1][wpPU-1]->Write(Form("hEPDA_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU     ));
                                hEPDrawT                   [wtrk-1][wppt-1][wpPU-1]->Write(Form("hEPDT_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU     ));
                                
                            }
                        }
                    }

                    hBinDist_gen_single->Write();
                    hEvent_Pass   ->Write();
                    hJet_Pass     ->Write();
                    hdeltaR       ->Write();
                    hdeltaJetEta->Write();
                    hdeltaJetPhi->Write();
                    hdeltaJetTheta->Write();
                    hT_jet_dau_eta->Write();
                    hT_jet_dau_phi->Write();
                    hT_jet_dau_pt->Write();
                    hNtrig->Write();
                    hJJT->Write();
                    hJJT_Cut->Write();
                    // hJet_Ptw_Eta_ave_nocut->Write();
                    hJet_Eta_ave_cut1->Write();
                    hJet_Eta_ave_cutR->Write();
                    hJet_Ptw_Eta_ave_cut1->Write();
                    // hJet_Ptw_Eta_ave_nocut->Write();
                    hJet_Ptw_Eta_ave_cutR->Write();
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

