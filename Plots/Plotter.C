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
#include "THnSparse.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
// #include "StJetTreeStruct.h"
#include <vector>

#pragma link C++ class vector<int> +;

#include "Bindef.h"

using namespace std;

#endif

void Plotter(){
	
	TFile *f[2];
	f[0] = new TFile("Output_NoQuenching.root");
	f[1] = new TFile("Output_Quenching.root");

	TH1D *LeadPt[2];
	TH1D *SubLeadPt[2];

	TH1D *LeadJetLeadPt[2];
	TH1D *LeadJetSubLeadPt[2];

	TH1D *SubLeadJetLeadPt[2];
	TH1D *SubLeadJetSubLeadPt[2];

	TString Name[2] = {"Unquenched", "Quenched"};
	TString Color[2] = {kBlack, kGreen-2};

	for (int i = 0; i <= 1; i++){
		LeadPt[i] = (TH1D *)f[i]->Get("LeadingParticlePt");
		SetName(LeadPt[i], TString("Lead p_{T} ") + Name[i]);
		SetColor(LeadPt[i], Color[i], 20);
		LeadPt[i]->Scale(1./LeadPt[i]->Integral());
		SubLeadPt[i] = (TH1D *)f[i]->Get("SubLeadingParticlePt");
		SetName(SubLeadPt[i], TString("Sub Lead p_{T} ") + Name[i]);
		SetColor(SubLeadPt[i], Color[i], 20);
		SubLeadPt[i]->Scale(1./LeadPt[i]->Integral());

		LeadJetLeadPt[i] = (TH1D *)f[i]->Get("LeadJetLeadParticle");
		SetName(LeadJetLeadPt[i], TString("Lead Jet Lead p_{T} ") + Name[i]);
		SetColor(LeadJetLeadPt[i], Color[i], 25);
		LeadJetLeadPt[i]->Scale(1./LeadPt[i]->Integral());
		LeadJetSubLeadPt[i] = (TH1D *)f[i]->Get("LeadJetSubLeadParticle");
		SetName(LeadJetSubLeadPt[i], TString("Lead Jet Sub Lead p_{T} ") + Name[i]);
		SetColor(LeadJetSubLeadPt[i], Color[i], 24);
		LeadJetSubLeadPt[i]->Scale(1./LeadPt[i]->Integral());

		SubLeadJetLeadPt[i] = (TH1D *)f[i]->Get("SubLeadJetLeadParticle");
		SetName(SubLeadJetLeadPt[i], TString("Sub Lead Jet Lead p_{T} ") + Name[i]);
		SetColor(SubLeadJetLeadPt[i], Color[i], 33);
		SubLeadJetLeadPt[i]->Scale(1./LeadPt[i]->Integral());
		SubLeadJetSubLeadPt[i] = (TH1D *)f[i]->Get("SubLeadJetSubLeadParticle");
		SetName(SubLeadJetSubLeadPt[i], TString("Sub Lead Jet Sub Lead p_{T} ") + Name[i]);
		SetColor(SubLeadJetSubLeadPt[i], Color[i], 34);
		SubLeadJetSubLeadPt[i]->Scale(1./LeadPt[i]->Integral());
	}

	TCanvas *c[3];
	for (int i = 0; i < 3; i++){
		c[i] = new TCanvas(Form("c_%i", i), Form("c_%i", i), 1000, 1000);
	}

	c[0]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	LeadPt[0]->Draw("EP");
	LeadPt[1]->Draw("EP SAME");

	SubLeadPt[0]->Draw("EP");
	SubLeadPt[1]->Draw("EP SAME");

}