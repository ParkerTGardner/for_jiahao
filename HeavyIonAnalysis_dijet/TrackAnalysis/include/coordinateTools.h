#ifndef COORDTOOLS
#define COORDTOOLS

#include "TVector3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>

//int getJetTrkIndx( int jetIndx, int trkIndx, int nTrk){
//  return jetIndx * nTrk + trkIndx;
//}

double ptWRTJet(TVector3 jet, TVector3 p){
  return p.Perp(jet);
}

double ptWRTJet( double jetPt, double jetEta, double jetPhi, double trkPt, double trkEta, double trkPhi){
  TVector3 j = TVector3(0,0,0);
  j.SetPtEtaPhi( jetPt, jetEta, jetPhi);

  TVector3 trk = TVector3(0,0,0);
  trk.SetPtEtaPhi( trkPt, trkEta, trkPhi);
  return ptWRTJet( j, trk);
}

double thetaWRTJet(TVector3 jet, TVector3 p){
  return p.Angle(jet);
}

double thetaWRTJet( double jetPt, double jetEta, double jetPhi, double trkPt, double trkEta, double trkPhi){
  TVector3 j = TVector3(0,0,0);
  j.SetPtEtaPhi( jetPt, jetEta, jetPhi);

  TVector3 trk = TVector3(0,0,0);
  trk.SetPtEtaPhi( trkPt, trkEta, trkPhi);
  return thetaWRTJet( j, trk);
}


double etaWRTJet(TVector3 jet, TVector3 p){
  double eta = -TMath::Log( TMath::Tan( thetaWRTJet(jet,p)/2.0));
  if(eta > 999.) eta = 999.;
  else if(eta < -999.) eta = -999.;
  return eta;
}

double etaWRTJet( double jetPt, double jetEta, double jetPhi, double trkPt, double trkEta, double trkPhi){
  TVector3 j = TVector3(0,0,0);
  j.SetPtEtaPhi( jetPt, jetEta, jetPhi);

  TVector3 trk = TVector3(0,0,0);
  trk.SetPtEtaPhi( trkPt, trkEta, trkPhi);
  return etaWRTJet( j, trk);
}


double phiWRTJet(TVector3 jet, TVector3 p){
  TVector3 pt = p-((p*jet.Unit())*(jet.Unit()));//pt vector
  TVector3 z = TVector3(0,0,1);
  TVector3 phiOrigin = jet.Unit().Cross((jet.Unit().Cross(z)));//vector that will be phi=0 (in plane of jet and beam line)
  double phi = pt.Angle(phiOrigin);//get phi from 0 to pi

  //determine sign of phi based on cross product of pt and origin
  if( (phiOrigin.Cross(pt.Unit()))*jet >= 0) return phi;
  else return -phi;
}

double phiWRTJet( double jetPt, double jetEta, double jetPhi, double trkPt, double trkEta, double trkPhi){
  TVector3 j = TVector3(0,0,0);
  j.SetPtEtaPhi( jetPt, jetEta, jetPhi);

  TVector3 trk = TVector3(0,0,0);
  trk.SetPtEtaPhi( trkPt, trkEta, trkPhi);
  return phiWRTJet( j, trk);
}

double diffvec(TVector3 jet, TVector3 p){
  TVector3 diff = jet - p;
  return diff.Mag();
}

double diffvec( double jetPt, double jetEta, double jetPhi, double trkPt, double trkEta, double trkPhi){
  TVector3 j = TVector3(0,0,0);
  j.SetPtEtaPhi( jetPt, jetEta, jetPhi);

  TVector3 trk = TVector3(0,0,0);
  trk.SetPtEtaPhi( trkPt, trkEta, trkPhi);
  return diffvec( j, trk);
}

TVector3 BeamBoost( double jetPt, double jetEta, double jetPhi, double trkPt, double trkEta, double trkPhi){
  // boost with Jet.pz
  TVector3 j0 = TVector3(0.0,0.0,0.0);
  j0.SetPtEtaPhi( jetPt, jetEta, jetPhi);
  //the boost func in root is inversed boost, so here should use -v
  TVector3 j (0.,0.,-j0.Z()/j0.Mag());

  TVector3 trk0 = TVector3(0.0,0.0,0.0);
  trk0.SetPtEtaPhi( trkPt, trkEta, trkPhi);
  TLorentzVector trk(0.,0.,0.,0.);
  trk.SetPtEtaPhiE( trkPt, trkEta, trkPhi, trk0.Mag());
  trk.Boost(j);
  return trk.Vect();
}

#endif

