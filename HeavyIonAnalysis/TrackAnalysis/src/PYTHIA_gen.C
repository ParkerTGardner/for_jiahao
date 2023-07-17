#define MyClass_cxx
// touch
#include "include/TrimPythia.h"
#include "include/coordinateTools.h"
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


void MyClass::Loop(int job, std::string fList){

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
    TH1D* hBinDist_gen_single = new TH1D("hBinDist_gen_single","hBinDist_gen_single",bin360,bin0,bin120);
    TH1D* hBinDist_reco_single = new TH1D("hBinDist_reco_single","hBinDist_reco_single",bin360,bin0,bin120);

    //2D Corr histograms

    TH2D* hEPDraw[trackbin][ptbin][PUbin];
    TH2D* hSignalShifted[trackbin][ptbin][PUbin];
    TH2D* hBckrndShifted[trackbin][ptbin][PUbin];


    TH1D* hBinDist_gen[trackbin];
    TH1D* hBinDist_reco[trackbin];
    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        hBinDist_gen[wtrk-1]    = new TH1D(Form("hBinDist_gen_%d",wtrk),Form("hBinDist_gen_%d",wtrk), bin360, bin0, bin120);
        hBinDist_reco[wtrk-1]    = new TH1D(Form("hBinDist_reco_%d",wtrk),Form("hBinDist_reco_%d",wtrk), bin360, bin0, bin120);
        for(int wppt = 1; wppt<ptbin+1; wppt++){
            for(int wpPU = 1; wpPU<PUbin+1; wpPU++){
                hBckrndShifted[wtrk-1][wppt-1][wpPU-1]      = new TH2D(Form("hBckrndS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form("hBckrndS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                hSignalShifted[wtrk-1][wppt-1][wpPU-1]      = new TH2D(Form("hSignalS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form("hSignalS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                hEPDraw[wtrk-1][wppt-1][wpPU-1]             = new TH2D(Form( "hEPDraw_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form( "hEPDraw_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);


            }
        }
    }


    std::cout << "Starting event loop" << std::endl;
    std::cout << "Total Number of Files in this Job: " << fileList.size() << std::endl;
    for(int f = 0; f<fileList.size(); f++){
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


                  if(genJetPt->size()==0) continue;
                  if(genJetChargedMultiplicity->size()==0) continue;


                  if(!F_eventpass(genJetPt, genJetPt->size(), jetPtCut_Event)){
                      continue;
                  }
                  int gjN = genJetPhi->size();

////if(made_bo) std::cout<< "made 3.4" << endl;

            hEvent_Pass->Fill(1);

////if(made_bo) std::cout<< "made 3.6" << endl;
            //ENTERING JET LOOP

            //for(int f = 0; f<fileList.size(); f++){
            //for (Long64_t ievent=0; ievent <nentries; ievent ++){
            for(int kjet=0; kjet < genJetPt->size(); kjet++){

                int ijet = kjet;
                int Gjet = kjet;
                long int NNtrk = (genDau_pt->at(ijet)).size();


                if( fabs((*genJetEta)[ijet]) > jetEtaCut ) continue;
                if( (*genJetPt)[ijet] < jetPtCut_Jet   ) continue;

                // filling distrivutions within track bins
                // ALSO VERY IMPORTANLTY changing the tkBool to 1 for this particular jet. This will be usefull later wen I create conditons for filling other historgams.
                int tkBool[trackbin] = {0};

                int n_G_ChargeMult_count =0;
                for(int  G_trk=0; G_trk < NNtrk; G_trk++ ){
                    if((*genDau_chg)[Gjet][G_trk] == 0) continue;
                    if(fabs((*genDau_pt)[Gjet][G_trk])  < 0.3)     continue;
                    if(fabs((*genDau_eta)[Gjet][G_trk]) > 2.4)     continue;
                    n_G_ChargeMult_count += 1;
                }//G_trk end

                hBinDist_gen_single            ->Fill(n_G_ChargeMult_count);

                //****** REMOVE REMOVE REMOVE REMOVE REMOVE ****
                //****** REMOVE REMOVE REMOVE REMOVE REMOVE ****
                //****** REMOVE REMOVE REMOVE REMOVE REMOVE ****

                //****** REMOVE REMOVE REMOVE REMOVE REMOVE ****
                //****** REMOVE REMOVE REMOVE REMOVE REMOVE ****
                //****** REMOVE REMOVE REMOVE REMOVE REMOVE ****

                continue;

                //****** REMOVE REMOVE REMOVE REMOVE REMOVE ****
                //****** REMOVE REMOVE REMOVE REMOVE REMOVE ****
                //****** REMOVE REMOVE REMOVE REMOVE REMOVE ****

                //****** REMOVE REMOVE REMOVE REMOVE REMOVE ****
                //****** REMOVE REMOVE REMOVE REMOVE REMOVE ****
                //****** REMOVE REMOVE REMOVE REMOVE REMOVE ****





                for(int i = 0; i < trackbin; i++){
                    //if((*chargedMultiplicity)[indicesR[kjet]] >= trackbinbounds[i] && (*chargedMultiplicity)[indicesR[kjet]] < trackbinboundsUpper[i]){
                    if(n_G_ChargeMult_count >= trackbinbounds[i] && n_G_ChargeMult_count < trackbinboundsUpper[i]){
                        tkBool[i] = 1;
                        hJet_Pass           ->Fill(i);
                        hBinDist_gen[i]         ->Fill(n_G_ChargeMult_count);     
                    }
                }

                int Ntrig[trackbin][ptbin] = {0};
                int NtrigM[trackbin][ptbin] = {0};
                int NtrigP[trackbin][ptbin] = {0};
                int A_ptBool[NNtrk][ptbin] = {0};

                // VERY IMPORTANT calculating daughter pt wrt to jet axis.
                // So this needs to be 2d vector, for pt bin and for daughter index.
                // In each case I create a true falsse for that daughter falling in to the specific pt bin.
                // Note this is NOT jet x daughter. It's pt bin x daughtr

                //for(int f = 0; f<fileList.size(); f++){
                //for (Long64_t ievent=0; ievent <nentries; ievent ++){
                //for(int ijet=0; ijet < jetCounter; ijet++){
                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    if((*genDau_chg)[ijet][A_trk] == 0) continue;
                    if(fabs((*genDau_pt)[ijet][A_trk])  < 0.3)     continue;
                    if(fabs((*genDau_eta)[ijet][A_trk]) > 2.4) continue;

                    double jet_dau_pt    =  ptWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);

                    //excluding outside outermost limits.
                    //if( jet_dau_pt < ptbinmin || jet_dau_pt > ptbinmax) continue;
                    //loop through pt bins and fill the boolean array
                    for(int i = 0; i < ptbin; i++){
                        if(jet_dau_pt >= ptbinbounds_lo[i] && jet_dau_pt < ptbinbounds_hi[i]){
                            A_ptBool[A_trk][i] = 1;
                        }
                    }
                    //in the case where the total jet mult and the individual daughter pt is acceptable for this track bin and pt bin, we increase the Ntrig count.
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

                }//here ends the boolean array creations.

                //continuation of main loops. Here is where the 2D Corr plots are created using the above booleans and 

                for(int i = 0; i < trackbin; i++){
                    for(int j = 0; j < ptbin; j++){
                        hNtrig->Fill(i,j,Ntrig[i][j]);
                    }
                }

                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    //this should be redundant if it passes the bools above? i guess it helps skip daughters faster. maybe i can reindex and run through the daughters quickly by aranging all the charged dauhghter sat the front.
                    if((*genDau_chg)[ijet][A_trk] == 0) continue;
                    if(fabs((*genDau_eta)[ijet][A_trk]) > 2.4) continue;
                    if(fabs((*genDau_pt)[ijet][A_trk])  < 0.3)     continue;

                    double jet_dau_pt    =  ptWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);

                    if(jet_dau_pt >3.0) continue;
                    //if(jet_dau_pt <0.5) continue;



                    double jet_dau_eta   = etaWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    double jet_dau_phi   = phiWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);

                    double jet_dau_theta = 2*ATan(Exp(-(jet_dau_eta)));
                    if(jet_dau_eta > track_eta_lim) continue;
                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin; j++){
                            if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                int k_PU=0;

                                hEPDraw[i][j][k_PU]->Fill(jet_dau_eta, jet_dau_phi, 1.0/( Ntrig[i][j] ));
                            }
                        }
                    }

                    if(A_trk == NNtrk - 1) continue;

                    for(long int T_trk=A_trk+1; T_trk< NNtrk; T_trk++ ){

                        if((*genDau_chg)[ijet][T_trk] == 0) continue;
                        if(fabs((*genDau_eta)[ijet][T_trk]) > 2.4) continue;
                        if(fabs((*genDau_pt)[ijet][T_trk])  < 0.3)     continue;
                        //if((*highPurity)[ijet][A_trk] == 0) continue;

                        double T_jet_dau_pt    =  ptWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][T_trk], (double)(*genDau_eta)[ijet][T_trk], (double)(*genDau_phi)[ijet][T_trk]);

                    if(T_jet_dau_pt >3.0) continue;
                    //if(T_jet_dau_pt <0.5) continue;

                        double T_jet_dau_eta   = etaWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][T_trk], (double)(*genDau_eta)[ijet][T_trk], (double)(*genDau_phi)[ijet][T_trk]);
                        double T_jet_dau_phi   = phiWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][T_trk], (double)(*genDau_eta)[ijet][T_trk], (double)(*genDau_phi)[ijet][T_trk]);
                        if(T_jet_dau_eta > track_eta_lim) continue;
                        double deltaEta = (jet_dau_eta - T_jet_dau_eta);
                        double deltaPhi = (TMath::ACos(TMath::Cos(jet_dau_phi - T_jet_dau_phi)));

                        for(        int i = 0; i < trackbin; i++){
                            for(    int j = 0; j < ptbin;    j++){ 
        


                                if(tkBool[i] + A_ptBool[A_trk][j] + A_ptBool[T_trk][j] == 3){
                                          hPairs->Fill(i,j);
                                          int k_PU=0;
                                          hSignalShifted[i][j][k_PU]->Fill(deltaEta, deltaPhi,                 ( 1./Ntrig[i][j]));
                                          hSignalShifted[i][j][k_PU]->Fill(-deltaEta, deltaPhi,                ( 1./Ntrig[i][j]));
                                          hSignalShifted[i][j][k_PU]->Fill(deltaEta, -deltaPhi,                ( 1./Ntrig[i][j]));
                                          hSignalShifted[i][j][k_PU]->Fill(-deltaEta, -deltaPhi,               ( 1./Ntrig[i][j])); 
                                          hSignalShifted[i][j][k_PU]->Fill( deltaEta,2*TMath::Pi() - deltaPhi, ( 1./Ntrig[i][j]));
                                          hSignalShifted[i][j][k_PU]->Fill(-deltaEta,2*TMath::Pi() - deltaPhi, ( 1./Ntrig[i][j]));

                                      //}}}
                                    // This is the mixed charge signal. Each duaghter will serve as a trigger so regular Ntrig suffices.
                                    // EPdraw will constitute sampling from EPD_P and EPD_M for the pseudo particles. 
                                    // So that I take random phi and eta from P and random Phi and Eta from M and make a opposite sign difference.
                                }
                            }
                        }
                    }//T_trk
                    }//A_trk
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


                                      for(int x = 0; x<XENT; x++){
                                          gRandom->SetSeed(0);
                                          double WEta1, WPhi1;//making the pseudoparticles
                                          hEPDraw[wtrk-1][wppt-1][wpPU-1]->GetRandom2(WEta1, WPhi1);
                                          A_ETA[x] = WEta1;
                                          A_PHI[x] = WPhi1;
                                      }
                                      for(long int i = 0; i < (XENT-1); i++){
                                          for(long int j = (i+1); j < XENT; j++){

                                              double WdeltaEta = (A_ETA[i]-A_ETA[j]);
                                              double WdeltaPhi = (TMath::ACos(TMath::Cos(A_PHI[i]-A_PHI[j])));

                                              hBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(WdeltaEta, WdeltaPhi, 1);//./XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(-WdeltaEta, WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(WdeltaEta, -WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(-WdeltaEta, -WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(WdeltaEta, 2*TMath::Pi() - WdeltaPhi, 1);//../XENT);
                                              hBckrndShifted[wtrk-1][wppt-1][wpPU-1]->Fill(-WdeltaEta,2*TMath::Pi() - WdeltaPhi, 1);//../XENT);

                                          }
                                      }
                                  }
                              }
                          }
                      //}}}




                    string subList = fList.substr(fList.size() - 3);

                    TFile* fS_tempA = new TFile(Form("pythia_batch_output/root_out/job_%s.root",subList.c_str()), "recreate");
                    for(int wtrk =1; wtrk <trackbin+1; wtrk++){
                        hBinDist_gen[wtrk-1]         ->Write();
                        for(int wppt =1; wppt <ptbin+1; wppt++){
                            for(int wpPU =1; wpPU<PUbin+1; wpPU++){

                                hSignalShifted             [wtrk-1][wppt-1][wpPU-1]->Write(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU    ));
                                hBckrndShifted             [wtrk-1][wppt-1][wpPU-1]->Write(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU    ));
                                hEPDraw                    [wtrk-1][wppt-1][wpPU-1]->Write(Form("hEPD_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU     ));
                            }
                        }
                    }

                    hBinDist_gen_single->Write();
                    hEvent_Pass   ->Write();
                    hJet_Pass     ->Write();

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

