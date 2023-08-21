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

    // TH2D* hPairs        = new TH2D("hPairs","hPairs",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);

    // TH2D* hNtrig        = new TH2D("hNtrig","hNtrig",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);

    

   
    TH1D* hkt_Pairs     = new TH1D("hkt_Pairs", "hkt_Pairs", ktbin, bin0, ktbin);

    TH1D* hkt_Ntrig_single     = new TH1D("hkt_Ntrig", "hkt_Ntrig", ktbin, bin0, ktbin);


    TH1D* hkt_Ntrig[ktbin];    
    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", trackbin,bin0,trackbin);
    TH1D* hJet_Pass     = new TH1D("hJet_Pass"  ,"hJet_Pass"  , trackbin,bin0,trackbin);
    TH1D* hJetPt = new TH1D("hJetPt"  ,"hJetPt" ,i150,i100,i3500);
    // TH1D* hdeltaR = new TH1D("hdeltaR", "hdeltaR", 100,0,5);
    
        
    TH1D* hMomSignalShifted[ktbin];
    TH1D* hMomBckrndShifted[ktbin];
    TH1D* kt_EPDraw[ktbin];
    TH3D* hkt_EPDraw = new TH3D(Form("hkt_EPDraw_kt"), Form("hkt_EPDraw_kt") , EPD_xb, EPD_xlo, EPD_xhi, EPD_yb, EPD_ylo, EPD_yhi, 100*ptbin, -ptbin, ptbin);


    for(int wkt =1; wkt < ktbin+1; wkt++){
        hMomSignalShifted[wkt-1] = new TH1D(Form("hMomSignalS_kt_%d",wkt), Form("hMomSignalS_kt_%d",wkt) , 30, 0, 1.2);
        hMomBckrndShifted[wkt-1] = new TH1D(Form("hMomBckrndS_kt_%d",wkt), Form("hMomBckrndS_kt_%d",wkt) , 30, 0, 1.2);
        kt_EPDraw[wkt-1] = new TH1D(Form("kt_EPDraw_kt_%d",wkt), Form("kt_EPDraw_kt_%d",wkt) , 30,0,1.2); 
        hkt_Ntrig[wkt-1] = new TH1D(Form("hkt_Ntrig_%d",wkt), Form("hkt_Ntrig_%d",wkt), ktbin, bin0, ktbin);
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


                
                int kt_tkbool[ktbin] = {0}; // This is ktbin
                int kt_Ntrig[ktbin] = {0};
            
                // Fill kt_Ntrig in jetA
                for(int  A_trk=0; A_trk < NNtrk1; A_trk++ ){
                
                    if((*genDau_chg)[ijet][A_trk] == 0) continue;
                    if(fabs((*genDau_eta)[ijet][A_trk]) > 2.4) continue;
                    if(fabs((*genDau_pt)[ijet][A_trk])  < 0.3)     continue;

                        //     daughter pt with respect to the jetA axis                 pt With Respect To JetA
                    double jet_dau_pt    =  ptWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);

                    if(jet_dau_pt >3.0) continue;// why we drop this

                    //     daughter eta with respect to the jet axis                 eta With Respect To Jet 
                    double jet_dau_eta   = etaWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    //     daughter phi with respect to the jet axis                 phi With Respect To Jet 
                    double jet_dau_phi   = phiWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);

                    hkt_EPDraw->Fill(jet_dau_eta, jet_dau_phi, jet_dau_pt, 1);


                    for(long int T_trk=A_trk+1; T_trk< NNtrk1; T_trk++ ){

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
                        for(int i = 0; i < ktbin; i++){

                            if((ktbinbounds_lo[i] < 0.5*fabs(jet_dau_pt+T_jet_dau_pt)) && (ktbinbounds_hi[i] >= 0.5*fabs(jet_dau_pt+T_jet_dau_pt))){
                                kt_Ntrig[i] += 1;
                            }
                        }
                        
                    }//T_trk AA


                }//A_trk

                for(int i = 0; i < ktbin; i++){
                    hkt_Ntrig[i]->Fill(kt_Ntrig[i]);
                    hJet_Pass           ->Fill(i);
                }

                //Main loop for corr func

                for(int  A_trk=0; A_trk < NNtrk1; A_trk++ ){
                    
                    if((*genDau_chg)[ijet][A_trk] == 0) continue;
                    if(fabs((*genDau_eta)[ijet][A_trk]) > 2.4) continue;
                    if(fabs((*genDau_pt)[ijet][A_trk])  < 0.3)     continue;

                        //     daughter pt with respect to the jet axis                 pt With Respect To Jet 
                    double jet_dau_pt    =  ptWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);


                        //if(jet_dau_pt >3.0) continue;

                    if(jet_dau_pt >3.0) continue;// why we drop this

                    //     daughter eta with respect to the jet axis                 eta With Respect To Jet 
                    double jet_dau_eta   = etaWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    //     daughter phi with respect to the jet axis                 phi With Respect To Jet 
                    double jet_dau_phi   = phiWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);

                    double jet_dau_theta = 2*ATan(Exp(-(jet_dau_eta)));



                    // This is for background
                    // For every new jet, we get kt_Ntrig for nomalization
                    // For every new daughter in this jet, we calculate every deltap and fill the corresponding ktbin
                    // Finally we get all particle's deltap for each ktbin---- It's really time-comsuming
                    // Then we choose 10 times pairs than signal's for each ktbin

                    int nBinsX = hkt_EPDraw->GetNbinsX();
                    int nBinsY = hkt_EPDraw->GetNbinsY();
                    int nBinsZ = hkt_EPDraw->GetNbinsZ();

                    for (int i = 1; i <= nBinsX; i++) {
                        for (int j = 1; j <= nBinsY; j++) {
                            for (int k = 1; k <= nBinsZ; k++) {
                                double binContent = hkt_EPDraw->GetBinContent(i, j, k);
                                double dau_eta_lo = hkt_EPDraw->GetXaxis()->GetBinLowEdge(i);
                                double dau_phi_lo = hkt_EPDraw->GetYaxis()->GetBinLowEdge(j);
                                double dau_pt_lo = hkt_EPDraw->GetZaxis()->GetBinLowEdge(k);
                                double dau_eta_wi = hkt_EPDraw->GetXaxis()->GetBinWidth(i);
                                double dau_phi_wi = hkt_EPDraw->GetYaxis()->GetBinWidth(j);
                                double dau_pt_wi = hkt_EPDraw->GetZaxis()->GetBinWidth(k);
                                double dau_eta_ave = dau_eta_lo+0.5*(dau_eta_wi);
                                double dau_phi_ave = dau_phi_lo+0.5*(dau_phi_wi);
                                double dau_pt_ave = dau_pt_lo+0.5*(dau_pt_wi);
                                
                                for (int l=0; l<ktbin; l++){

                                    if((ktbinbounds_lo[l] < 0.5*fabs(dau_pt_ave+jet_dau_pt)) && (ktbinbounds_hi[l] >= 0.5*fabs(dau_pt_ave+jet_dau_pt))){

                                        kt_EPDraw[l]->Fill(diffvec(dau_eta_ave,dau_phi_ave,dau_pt_ave,jet_dau_eta,jet_dau_phi,jet_dau_pt), (double)((double)(binContent)/kt_Ntrig[l]));

                                    }

                                }

                                
                                    
                                
                            }
                        }
                    }

                    // corr(dau_A,dau_A)
                    for(long int T_trk=A_trk+1; T_trk< NNtrk1; T_trk++ ){

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
                        
                        double deltaPt  = diffvec(jet_dau_pt, jet_dau_eta, jet_dau_phi, T_jet_dau_pt, T_jet_dau_eta, T_jet_dau_phi);

                        for(int i = 0; i < ktbin; i++){
                            if((ktbinbounds_lo[i] < 0.5*fabs(jet_dau_pt+T_jet_dau_pt)) && (ktbinbounds_hi[i] >= 0.5*fabs(jet_dau_pt+T_jet_dau_pt))){
                                hkt_Pairs->Fill(i);
                                hMomSignalShifted[i]->Fill(deltaPt,                         (double)((double)(1.0)/(kt_Ntrig[i])));
                            }
                        }

                        
                    }//T_trk AA


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
                        
                        


                        for(int wkt = 0; wkt < ktbin; wkt++){
                            std::cout << wkt  << "/" << ktbin << std::endl;
                            //Nent is the number of pairs in the signal which we will try to 10x
                            //Xent is the number of pseudoparticles requried such that when we build the pairs nCp = Xent CHOOSE 2 will give 
                            //us 10 times as many pairs as we have in the signal histogrm.

                            long int NENT =  hkt_Pairs->GetBinContent(wkt);
                            long int XENT =  ((1+floor(sqrt(1+(4*2*backMult*NENT))))/2) ;
                            

                            for(int x = 0; x<XENT; x++){
                                gRandom->SetSeed(0);
                                double WdeltaPt = kt_EPDraw[wkt]->GetRandom();//making the pseudoparticles
                                
                                hMomBckrndShifted[wkt]->Fill(WdeltaPt, 1);
                            }

                        }
                      //}}}




                    string subList = fList.substr(fList.size() - 3);

                    TFile* fS_tempA = new TFile(Form("pythia_batch_output/root_out_copy/dijob_%s.root",subList.c_str()), "recreate");
                   

                    for(int wkt =1; wkt<ktbin+1; wkt++){
                        
                        hMomSignalShifted [wkt-1]-> Write(Form("hMomSigS_%d_to_%d", (int)(10*ktbinbounds_lo[wkt-1]), (int)(10*ktbinbounds_hi[wkt-1])));
                        hMomBckrndShifted [wkt-1]-> Write(Form("hMomBckS_%d_to_%d", (int)(10*ktbinbounds_lo[wkt-1]), (int)(10*ktbinbounds_hi[wkt-1])));
                        hkt_Ntrig[wkt-1]-> Write(Form("hkt_Ntrig_%d_to_%d", (int)(10*ktbinbounds_lo[wkt-1]), (int)(10*ktbinbounds_hi[wkt-1])));
                    }

                    
                    hEvent_Pass   ->Write();
                    hJet_Pass     ->Write();
                    hkt_Ntrig_single     ->Write();
                    // hdeltaR       ->Write();
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

