#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TChain.h"
#include "TObject.h"
#include "TClonesArray.h"
// #include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include <TLorentzVector.h>
#ifndef __CINT__
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "Riostream.h"
#include <cstdlib>
#include "TH3F.h"
#include "TH2F.h"
#include "THn.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TLegend.h"
// #include "StJetTreeStruct.h"
#include <vector>

#endif


TH1D *ProcessSpectraHistogram(TH1D *h, TH1D *k = NULL){
    TH1D *R = (TH1D *)h->Clone();
    for(int i = 1;i<h->GetNbinsX()+1;i++){
        double val = R->GetBinContent(i);
        double er = R->GetBinError(i);
        double width = R->GetBinWidth(i);
        double center = fabs(R->GetBinCenter(i));

        R->SetBinContent(i,val/width/2./1.2/TMath::Pi()/center/0.035);
        R->SetBinError(i,er/width/2./1.2/TMath::Pi()/center/0.035);
    }

    if (k) {R->Divide(k); cout << "Corrected by gen efficiency" << endl;}

    return R;
}

TH1 *SetName(TH1 *h, TString Name){
  h->SetNameTitle(Name.Data(), Name.Data());
  return h;
}

TH2 *SetName(TH2 *h, TString Name){
  h->SetNameTitle(Name.Data(), Name.Data());
  return h;
}

TH1 *SetAxisTitles(TH1 *h, TString xaxis = "", TString yaxis = ""){
  h->GetXaxis()->SetTitle(xaxis.Data());
  h->GetYaxis()->SetTitle(yaxis.Data());
  return h;
}

TH2 *SetAxisTitles(TH2 *h, TString xaxis = "", TString yaxis = ""){
  h->GetXaxis()->SetTitle(xaxis.Data());
  h->GetYaxis()->SetTitle(yaxis.Data());
  return h;
}

TH1 *SetColor(TH1 *h, int color, int marker = 20){
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  return h;
}

TH1D *Rebin(TH1D *h, TString Name, int nbins, double *bins){
  TH1D *k = (TH1D *)h->Rebin(nbins, Name.Data(), bins);  
  return k;
}

double CalculateChi2(TH1 *Unfolded, TH1 *MC){
    if (Unfolded->GetNbinsX() != MC->GetNbinsX()) {throw invalid_argument( "Histograms do not have the same binning." );}

    double chi2 = 0;

    for (int i = 1; i <= MC->GetNbinsX(); i++){
        if (MC->GetBinContent(i) != 0) chi2 += pow((Unfolded->GetBinContent(i) - MC->GetBinContent(i)), 2)/MC->GetBinContent(i);
    }

    return chi2; 
}

double CalculateChi2(TH2 *Unfolded, TH2 *MC){
    if (Unfolded->GetNbinsX() != MC->GetNbinsX() || Unfolded->GetNbinsY() != MC->GetNbinsY()) {throw invalid_argument( "Histograms do not have the same binning." );}

    double chi2 = 0;

    for (int i = 1; i <= MC->GetNbinsX(); i++){
      for (int j = 1; j <= MC->GetNbinsY(); j++){
        if (MC->GetBinContent(i,j) != 0) chi2 += pow((Unfolded->GetBinContent(i,j) - MC->GetBinContent(i,j)), 2)/MC->GetBinContent(i,j);
      }
    }

    return chi2; 
}


const double R = 0.4;
const double deltar = 0.05;
const int numberofbins = R/deltar;


// Important functions that are used repeatedly in the calculations

// function to calculate relative phi between 2 objects and shift between 0 and 2pi 
//___________________________________________________________________________________________
Double_t dPhi(Double_t phi1, Double_t phi2) {
  Double_t deltaPhi;
  deltaPhi = abs(phi1 - phi2); //TODO absolute values
  if (deltaPhi>(2*TMath::Pi()))  deltaPhi-=2*(TMath::Pi());
  if (deltaPhi<(0*TMath::Pi())) deltaPhi+=2*(TMath::Pi()); 

  if (deltaPhi > TMath::Pi()) deltaPhi= 2*(TMath::Pi()) - deltaPhi;
  return deltaPhi;   // dphi in [0, 2Pi]
}

Double_t standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}

// function to calculate relative eta between 2 objects
//___________________________________________________________________________________________
Double_t dEta(Double_t eta1, Double_t eta2) {
  Double_t deltaEta;
  deltaEta = eta1 - eta2;

  return deltaEta;
}

// function to calculate relative eta between 2 objects
//___________________________________________________________________________________________
Double_t dR(Double_t delphi, Double_t deleta) {
  Double_t dRad;
  dRad = TMath::Sqrt(pow(delphi,2) + pow(deleta,2));

  return dRad;
}

// function to calculate invariant mass of a pair of objects
//___________________________________________________________________________________________
Double_t Mass(Double_t m1, Double_t m2, Double_t E1, Double_t E2, Double_t px1, Double_t px2, Double_t py1, Double_t py2, Double_t pz1, Double_t pz2) {
  Double_t m;
  m = TMath::Sqrt(pow(m1,2) + pow(m2,2) + 2*(E1*E2 - px1*px2 - py1*py2 - pz1*pz2));

  return m;
}

// function to calculate invariant mass of a pair of objects
//___________________________________________________________________________________________
Double_t pTforD0(Double_t px1, Double_t px2, Double_t py1, Double_t py2) {
  Double_t m;
  m = TMath::Sqrt(pow(px1 + px2, 2) + pow(py1 + py2, 2));

  return m;
}

Double_t pX(Double_t pT, Double_t phi){
  Double_t px;
  px = pT*TMath::Cos(phi);

  return px;
}

Double_t pY(Double_t pT, Double_t phi){
  Double_t py;
  py = pT*TMath::Sin(phi);

  return py;
}

Double_t p(Double_t E, Double_t m){
   
  Double_t pmag;
  pmag = TMath::Sqrt(pow(E,2) - pow(m,2));

  return pmag;
}

Double_t p(Double_t px, Double_t py, Double_t pz){
   
  Double_t pmag;
  pmag = TMath::Sqrt(pow(px,2) + pow(py,2) + pow(pz, 2));

  return pmag;
}

Double_t pZ(Double_t E, Double_t m, Double_t eta){
   
  Double_t pz;
  pz = p(E,m)*(TMath::Exp(2*eta)-1)/(TMath::Exp(2*eta)+1);

  return pz;
}



