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

    TH2D* hPairs        = new TH2D("hPairs","hPairs",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);

    TH2D* hNtrig        = new TH2D("hNtrig","hNtrig",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);


    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", trackbin,bin0,trackbin);
    TH1D* hJet_Pass     = new TH1D("hJet_Pass"  ,"hJet_Pass"  , trackbin,bin0,trackbin);
    TH1D* hJetPt = new TH1D("hJetPt"  ,"hJetPt" ,i150,i100,i3500);
    TH1D* hBinDist_gen_single = new TH1D("hBinDist_gen_single","hBinDist_gen_single",bin360,bin0, bin120);
    TH1D* hBinDist_reco_single = new TH1D("hBinDist_reco_single","hBinDist_reco_single",bin360,bin0,bin120);

    


    
    TH1D* hJt = new TH1D("hJt","hJt", 60,0,3);
    TH1D* hJtA[trackbin][ptbin][PUbin]; 
    // TH1D* hEtaJetA =new TH1D("hEtaJetA","hEtaJetA",50,-1,1);
    
   

    TH1D* hBinDist_gen[trackbin];
    TH1D* hBinDist_reco[trackbin];
    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        hBinDist_gen[wtrk-1]    = new TH1D(Form("hBinDist_gen_%d",wtrk),Form("hBinDist_gen_%d",wtrk), bin360, bin0, bin120);
        hBinDist_reco[wtrk-1]   = new TH1D(Form("hBinDist_reco_%d",wtrk),Form("hBinewnDist_reco_%d",wtrk), bin360, bin0, bin120);
        
        for(int wppt = 1; wppt<ptbin+1; wppt++){
            for(int wpPU = 1; wpPU<PUbin+1; wpPU++){

                hJtA[wtrk-1][wppt-1][wpPU-1]     = new TH1D(Form("hJtA_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU),Form("hJtA_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU),60,0,3);
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

                TVector3 JetA;
                JetA.SetPtEtaPhi((*genJetPt)[ijet],(*genJetEta)[ijet],(*genJetPhi)[ijet]);
                
                if( fabs(JetA.Eta()) > jetEtaCut ) continue;
                if( JetA.Perp() < jetPtCut_Jet   ) continue;


                 //    count the trks of A  
                int n_G_ChargeMult_count  = 0;
                
                
                for(int  G_trk=0; G_trk < NNtrk1; G_trk++ ){
                    if((*genDau_chg)[Gjet][G_trk] == 0) continue;
                    if(fabs((*genDau_pt)[Gjet][G_trk])  < 0.3)     continue;
                    if(fabs((*genDau_eta)[Gjet][G_trk]) > 2.4)     continue;
                    n_G_ChargeMult_count += 1;
                }
                // if (n_G_ChargeMult_count1<30) continue;
                
                
               hBinDist_gen_single            ->Fill(n_G_ChargeMult_count);

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
                    dau_A.SetPtEtaPhi((double)(*genDau_pt)[ijet][A_trk],(double)(*genDau_eta)[ijet][A_trk],(double)(*genDau_phi)[ijet][A_trk]);
                    
                    
                    if((*genDau_chg)[ijet][A_trk] == 0) continue;
                    if(fabs(dau_A.Eta()) > 2.4)        continue;
                    if(fabs(dau_A.Perp())  < 0.3)      continue;

                    //     daughter pt with respect to the jet axis                 pt With Respect To Jet 
                    double jet_dau_pt    =  ptWRTJet(JetA, dau_A);

                    if(jet_dau_pt >3.0) continue;


                    // boosted:
                    // double jet_dau_pt    =  ptWRTJet(JetA, dau_A)

                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin; j++){
                            if(tkBool[i] + (jet_dau_pt >= ptbinbounds_lo[j] && jet_dau_pt < ptbinbounds_hi[j]) == 2){
                                hJtA[i][j]->Fill(jet_dau_pt);
                            }
                        }
                    }

                    hJt->Fill(jet_dau_pt);
                        
                }

                // Here should be the final Ntrig for jetA
             
            }//kjet
                    }//Event
                    fFile->Close();
                    }//File
                    int backMult =10;
                    //filling background
                    //{{{
                    //


                      //}}}




                    string subList = fList.substr(fList.size() - 3);

                    TFile* fS_tempA = new TFile(Form("pythia_batch_output/root_out/dijob_%s.root",subList.c_str()), "recreate");
                    
                    for(int wtrk =1; wtrk <trackbin+1; wtrk++){
                        hBinDist_gen[wtrk-1]         ->Write();
                        for(int wppt =1; wppt <ptbin+1; wppt++){
                            for(int wpPU =1; wpPU<PUbin+1; wpPU++){

                                hJtA[wtrk-1][wppt-1][wpPU-1] ->Write(Form("hJtA_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU));
                            }
                        }
                    }


                    

                    hBinDist_gen_single->Write();
                    hEvent_Pass   ->Write();
                    hJet_Pass     ->Write();
                    hJt->Write();

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