#define MyClass_cxx
// touch



#include "include/TrimPythia_data.h"
//below allows for rotation into jet frame
#include "include/coordinateTools.h"
//defines constants
#include "include/sherpa_constants_cumulant.h"

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

    // TH1D* hBinDist_corrected = new TH1D("hBinDist_corrected","hBinDist_corrected",bin360,bin0, bin120);
    // TH1D* hBinDist_uncorrected = new TH1D("hBinDist_uncorrected","hBinDist_uncorrected",bin360,bin0, bin120);

    TH1D* hjet_avg_numerator_two = new TH1D("hjet_avg_numerator_two", "hjet_avg_numerator_two" ,trackbin , trackbinEdge);
    TH1D* hjet_avg_denominat_two = new TH1D("hjet_avg_denominat_two", "hjet_avg_denominat_two" ,trackbin , trackbinEdge);
    TH1D* hjet_avg_numerator_four = new TH1D("hjet_avg_numerator_four", "hjet_avg_numerator_four" ,trackbin , trackbinEdge);
    TH1D* hjet_avg_denominat_four = new TH1D("hjet_avg_denominat_four", "hjet_avg_denominat_four" ,trackbin , trackbinEdge);

    
   

    TH1D* hBinDist_gen[trackbin];
    TH1D* hPhiDrawA[trackbin];
    TH1D* hPhiDrawT[trackbin];
    TH1D* hEta[trackbin];
    TH1D* hPhi_EP[trackbin];
    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        hBinDist_gen[wtrk-1]    = new TH1D(Form("hBinDist_gen_%d",wtrk),Form("hBinDist_gen_%d",wtrk), bin360, bin0, bin120);
        hPhiDrawA[wtrk-1] = new TH1D(Form("hPhiDrawA_%d",wtrk),Form("hPhiDrawA_%d",wtrk), 1000, -TMath::Pi(), TMath::Pi());
        hPhiDrawT[wtrk-1] = new TH1D(Form("hPhiDrawT_%d",wtrk),Form("hPhiDrawT_%d",wtrk), 1000, -TMath::Pi(), TMath::Pi());
        hEta[wtrk-1] = new TH1D(Form("hEta_%d",wtrk),Form("hEta_%d",wtrk), 1000, 0, 8);
        hPhi_EP[wtrk-1] = new TH1D(Form("hPhi_EP_%d",wtrk),Form("hPhi_EP_%d",wtrk), 1000, -2*TMath::Pi(), 2*TMath::Pi());
        
    }

    double n_harm = 2.0;
    double v2 = 0.15;


  

    


  
    


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



            
            
            
            
            
            // def random Event plane to calculate v2 for each jet
            // phi*' =  phi* - Psi
            // weight = (1+2v2*cos(phi*'))
            TRandom3 randGenerator(0); 
            double Psi = randGenerator.Uniform(-TMath::Pi(), TMath::Pi());
            
            //ENTERING JET LOOP
            for(int kjet=0; kjet < jetPt->size(); kjet++){

                int ijet = kjet;
                int Gjet = kjet;
                long int NNtrk1 = (dau_pt->at(ijet)).size();

                TVector3 JetA;
                JetA.SetPtEtaPhi((*jetPt)[ijet],(*jetEta)[ijet],(*jetPhi)[ijet]);
                // TLorentzVector JetA_4 (JetA, JetA.Mag());
                
                if( fabs(JetA.Eta()) > jetEtaCut ) continue;
                if( JetA.Perp() < jetPtCut_Jet   ) continue;


                 //    count the trks of A  
                double n_G_ChargeMult_count  = 0.0;


                double jet_HLT_weight = 1.0;
                if((*jetPt)[ijet] < 880 && (*jetPt)[ijet] > 550){
                    jet_HLT_weight = 1.0/ (hdid500->GetBinContent(hdid500->FindBin((*jetPt)[ijet]) ) );
                }

                 //    count the trks of A  
                double n_G_ChargeMult_count_corrected = 0.0;
                
                
                for(int  G_trk=0; G_trk < NNtrk1; G_trk++ ){
                    if((*dau_chg)[Gjet][G_trk] == 0) continue;
                    if(fabs((*dau_pt)[Gjet][G_trk])  < 0.3)     continue;
                    if(fabs((*dau_eta)[Gjet][G_trk]) > 2.4)     continue;
                    n_G_ChargeMult_count += 1;
                    double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][G_trk] , (*dau_eta)[ijet][    G_trk] )));
                    n_G_ChargeMult_count_corrected += (1.0/nUnc_weight);
                }

                


                int tkBool[trackbin] = {0};
                

                for(int i = 0; i < trackbin; i++){
                //if((*chargedMultiplicity)[indicesR[kjet]] >= trackbinbounds[i] && (*chargedMultiplicity)[indicesR[kjet]] < trackbinboundsUpper[i]){
                    if(n_G_ChargeMult_count >= trackbinbounds[i] && n_G_ChargeMult_count < trackbinboundsUpper[i]){
                        tkBool[i] = 1;
                        // hJet_Pass           ->Fill(i);
                        hBinDist_gen[i]         ->Fill(n_G_ChargeMult_count);
                    }
                }



                std::complex<double> Q_all2A (0, 0);
                std::complex<double> Q_all4A (0, 0);
                std::complex<double> Q_all2T (0, 0);
                std::complex<double> Q_all4T (0, 0);
                // std::complex<double> Q_all4T (0, 0);

                int M = 0;
                int N = 0;
                double weight_sum_A = 0.0;
                double weight_sum_T = 0.0;
                double weight_sum_A_sqr = 0.0;
                double weight_sum_T_sqr = 0.0;

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
                    

                    double jet_dau_eta   = etaWRTJet(JetA, dau_A0);
                    //     daughter phi with respect to the jet axis                 phi With Respect To Jet 
                    double jet_dau_phi   = phiWRTJet(JetA, dau_A0) ;




                    if(jet_dau_pt  <0.7) continue;

                    TVector3 EP;
                    EP.SetXYZ(TMath::Cos(Psi),TMath::Sin(Psi),0);
                    TVector3 dau;
                    dau.SetXYZ(TMath::Cos(jet_dau_phi),TMath::Sin(jet_dau_phi),0);
                    TVector3 z;
                    z.SetXYZ(0.,0.,1.);
                    double phi_EP0 = TMath::ACos(EP*dau);
                    double phi_EP;
                    if((EP.Cross(dau))*z >= 0) phi_EP = phi_EP0;
                    else phi_EP = -phi_EP0;

                    double weight = 1.0 / (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                        
                    double phi = jet_dau_phi;


                    
                    std::complex<double> Q_part2 (weight*TMath::Cos(n_harm*phi) , weight*TMath::Sin(n_harm*phi));
                    std::complex<double> Q_part4 (weight*weight*TMath::Cos(2.0*n_harm*phi) , weight*weight*TMath::Sin(2.0*n_harm*phi));
                    

                    if ((jet_dau_eta>0.86) && (jet_dau_eta<2.25)){

                        Q_all2A = Q_all2A + Q_part2;
                        Q_all4A = Q_all4A + Q_part4;
                        M++;
                        weight_sum_A += weight;
                        weight_sum_A_sqr += weight*weight;
                    }

                    else if ((jet_dau_eta>2.25) && (jet_dau_eta<5.00)){

                        Q_all2T = Q_all2T + Q_part2;
                        Q_all4T = Q_all4T + Q_part4;
                        N++;
                        weight_sum_T += weight;
                        weight_sum_T_sqr += weight*weight;
                    }
                    else continue;

                    for(int i = 0; i < trackbin; i++){
                    //if((*chargedMultiplicity)[indicesR[kjet]] >= trackbinbounds[i] && (*chargedMultiplicity)[indicesR[kjet]] < trackbinboundsUpper[i]){
                        if(tkBool[i] == 1){
                            hEta[i]->Fill(jet_dau_eta, jet_HLT_weight*weight);
                            hPhi_EP[i]->Fill(phi_EP, jet_HLT_weight*weight);
                        }
                    }
                        

                        
                }


                // mult of subevent A, B >3
                if (M<3) continue;
                if (N<3) continue;

                
                // double particle_twoAT = Q_all_absAT / (weight_sum_A*weight_sum_T);
                double Q_all_absAT = std::real(Q_all2A*std::conj(Q_all2T));
                double weight_two_AT = (weight_sum_A*weight_sum_T);



                for(int i = 0; i < trackbin; i++){
                    if(tkBool[i] == 1){

                    
                    hjet_avg_numerator_two->Fill(i, jet_HLT_weight*Q_all_absAT);
                    hjet_avg_denominat_two->Fill(i,(jet_HLT_weight*weight_two_AT));

                    }
                }

               

                double S_A = weight_sum_A*weight_sum_A - weight_sum_A_sqr;
                double S_T = weight_sum_T*weight_sum_T - weight_sum_T_sqr;
                if (S_A*S_T == 0) continue;
                double four_numerator = std::real( (Q_all2A*Q_all2A-Q_all4A)*std::conj((Q_all2T*Q_all2T-Q_all4T)) );
                double particle_four = std::real( (Q_all2A*Q_all2A-Q_all4A)*std::conj((Q_all2T*Q_all2T-Q_all4T)) ) / (S_A * S_T);
                double weight_four = (S_A * S_T);

                for(int i = 0; i < trackbin; i++){
                    if(tkBool[i] == 1){

                    hjet_avg_numerator_four->Fill(i,jet_HLT_weight*four_numerator );
                    hjet_avg_denominat_four->Fill(i,(jet_HLT_weight*weight_four));

                    }
                }


                
             
            }//kjet
                    }//Event
                    fFile->Close();
                    }//File

    
    string subList = fList.substr(fList.size() - 3);
    TFile* fS_tempA = new TFile(Form("pythia_batch_data_output/root_out_3/dijob_%s.root",subList.c_str()), "recreate");
    hjet_avg_numerator_two ->Write();
    hjet_avg_denominat_two ->Write();
    hjet_avg_numerator_four->Write();
    hjet_avg_denominat_four ->Write();
    
    // hRand_jet_avg_numerator_two ->Write();
    // hRand_jet_avg_denominat_two ->Write();
    // hRand_jet_avg_numerator_four->Write();
    // hRand_jet_avg_denominat_four ->Write();

    for(int i=0; i<trackbin; i++){
        hBinDist_gen[i]->Write();
        hEta[i] -> Write();
        hPhi_EP[i] -> Write();
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