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

void Plotter(){

	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	
	TFile *f[2];
	f[0] = new TFile("Output_NoQuenching_pthat15.root");
	f[1] = new TFile("Output_Quenching_pthat15.root");

	TH1D *LeadJetPt[2];
	TH1D *SubLeadJetPt[2];

	TH1D *LeadPt[2];
	TH1D *SubLeadPt[2];

	TH1D *LeadJetLeadPt[2];
	TH1D *LeadJetSubLeadPt[2];

	TH1D *SubLeadJetLeadPt[2];
	TH1D *SubLeadJetSubLeadPt[2];

	TH1D *LeadJetPtW[2];
	TH1D *SubLeadJetPtW[2];

	TH1D *LeadPtW[2];
	TH1D *SubLeadPtW[2];

	TH1D *LeadJetLeadPtW[2];
	TH1D *LeadJetSubLeadPtW[2];

	TH1D *SubLeadJetLeadPtW[2];
	TH1D *SubLeadJetSubLeadPtW[2];

	TString Name[2] = {"Unquenched", "Quenched"};
	int Color[2] = {kBlack, kGreen-2};

	for (int i = 0; i <= 1; i++){
		LeadPt[i] = (TH1D *)f[i]->Get("LeadingParticlePt");
		LeadPt[i]->Rebin(4);
		SetName(LeadPt[i], TString("Lead p_{T} ") + Name[i]);
		SetColor(LeadPt[i], Color[i], 20);

		double scalefactor = LeadPt[i]->Integral();
		LeadPt[i]->Scale(1./scalefactor);

		SubLeadPt[i] = (TH1D *)f[i]->Get("SubLeadingParticlePt");
		SubLeadPt[i]->Rebin(4);
		SetName(SubLeadPt[i], TString("Sub Lead p_{T} ") + Name[i]);
		SetColor(SubLeadPt[i], Color[i], 21);
		SubLeadPt[i]->Scale(1./scalefactor);

		LeadJetPt[i] = (TH1D *)f[i]->Get("LeadingJetPt");
		LeadJetPt[i]->Rebin(4);
		SetName(LeadJetPt[i], TString("Lead Jet p_{T} ") + Name[i]);
		SetColor(LeadJetPt[i], Color[i], 45);
		LeadJetPt[i]->Scale(1./scalefactor);

		SubLeadJetPt[i] = (TH1D *)f[i]->Get("SubLeadingJetPt");
		SubLeadJetPt[i]->Rebin(4);
		SetName(SubLeadJetPt[i], TString("Sub Lead Jet p_{T} ") + Name[i]);
		SetColor(SubLeadJetPt[i], Color[i], 47);
		SubLeadJetPt[i]->Scale(1./scalefactor);

		LeadJetLeadPt[i] = (TH1D *)f[i]->Get("LeadJetLeadParticle");
		LeadJetLeadPt[i]->Rebin(4);
		SetName(LeadJetLeadPt[i], TString("Lead Jet Lead p_{T} ") + Name[i]);
		SetColor(LeadJetLeadPt[i], Color[i], 45);
		LeadJetLeadPt[i]->Scale(1./scalefactor);

		LeadJetSubLeadPt[i] = (TH1D *)f[i]->Get("LeadJetSubLeadParticle");
		LeadJetSubLeadPt[i]->Rebin(4);
		SetName(LeadJetSubLeadPt[i], TString("Lead Jet Sub Lead p_{T} ") + Name[i]);
		SetColor(LeadJetSubLeadPt[i], Color[i], 47);
		LeadJetSubLeadPt[i]->Scale(1./scalefactor);

		SubLeadJetLeadPt[i] = (TH1D *)f[i]->Get("SubLeadJetLeadParticle");
		SubLeadJetLeadPt[i]->Rebin(4);
		SetName(SubLeadJetLeadPt[i], TString("Sub Lead Jet Lead p_{T} ") + Name[i]);
		SetColor(SubLeadJetLeadPt[i], Color[i], 22);
		SubLeadJetLeadPt[i]->Scale(1./scalefactor);

		SubLeadJetSubLeadPt[i] = (TH1D *)f[i]->Get("SubLeadJetSubLeadParticle");
		SubLeadJetSubLeadPt[i]->Rebin(4);
		SetName(SubLeadJetSubLeadPt[i], TString("Sub Lead Jet Sub Lead p_{T} ") + Name[i]);
		SetColor(SubLeadJetSubLeadPt[i], Color[i], 23);
		SubLeadJetSubLeadPt[i]->Scale(1./scalefactor);

		LeadPtW[i] = (TH1D *)f[i]->Get("LeadingParticlePtW");
		LeadPtW[i]->Rebin(4);
		SetName(LeadPtW[i], TString("Lead p_{T} W") + Name[i]);
		SetColor(LeadPtW[i], Color[i], 20);

		double scalefactorW = LeadPtW[i]->Integral();
		LeadPtW[i]->Scale(1./scalefactorW);

		SubLeadPtW[i] = (TH1D *)f[i]->Get("SubLeadingParticlePtW");
		SubLeadPtW[i]->Rebin(4);
		SetName(SubLeadPtW[i], TString("Sub Lead p_{T} W") + Name[i]);
		SetColor(SubLeadPtW[i], Color[i], 21);
		SubLeadPtW[i]->Scale(1./scalefactorW);

		LeadJetPtW[i] = (TH1D *)f[i]->Get("LeadingJetPtW");
		LeadJetPtW[i]->Rebin(4);
		SetName(LeadJetPtW[i], TString("Lead Jet p_{T} W") + Name[i]);
		SetColor(LeadJetPtW[i], Color[i], 45);
		LeadJetPtW[i]->Scale(1./scalefactorW);

		SubLeadJetPtW[i] = (TH1D *)f[i]->Get("SubLeadingJetPtW");
		SubLeadJetPtW[i]->Rebin(4);
		SetName(SubLeadJetPtW[i], TString("Sub Lead Jet p_{T} W") + Name[i]);
		SetColor(SubLeadJetPtW[i], Color[i], 47);
		SubLeadJetPtW[i]->Scale(1./scalefactorW);

		LeadJetLeadPtW[i] = (TH1D *)f[i]->Get("LeadJetLeadParticleW");
		LeadJetLeadPtW[i]->Rebin(4);
		SetName(LeadJetLeadPtW[i], TString("Lead Jet Lead p_{T} W") + Name[i]);
		SetColor(LeadJetLeadPtW[i], Color[i], 45);
		LeadJetLeadPtW[i]->Scale(1./scalefactorW);

		LeadJetSubLeadPtW[i] = (TH1D *)f[i]->Get("LeadJetSubLeadParticleW");
		LeadJetSubLeadPtW[i]->Rebin(4);
		SetName(LeadJetSubLeadPtW[i], TString("Lead Jet Sub Lead p_{T} W") + Name[i]);
		SetColor(LeadJetSubLeadPtW[i], Color[i], 47);
		LeadJetSubLeadPtW[i]->Scale(1./scalefactorW);

		SubLeadJetLeadPtW[i] = (TH1D *)f[i]->Get("SubLeadJetLeadParticleW");
		SubLeadJetLeadPtW[i]->Rebin(4);
		SetName(SubLeadJetLeadPtW[i], TString("Sub Lead Jet Lead p_{T} W") + Name[i]);
		SetColor(SubLeadJetLeadPtW[i], Color[i], 22);
		SubLeadJetLeadPtW[i]->Scale(1./scalefactorW);

		SubLeadJetSubLeadPtW[i] = (TH1D *)f[i]->Get("SubLeadJetSubLeadParticleW");
		SubLeadJetSubLeadPtW[i]->Rebin(4);
		SetName(SubLeadJetSubLeadPtW[i], TString("Sub Lead Jet Sub Lead p_{T} W") + Name[i]);
		SetColor(SubLeadJetSubLeadPtW[i], Color[i], 23);
		SubLeadJetSubLeadPtW[i]->Scale(1./scalefactorW);
	}

	TCanvas *c[8];
	for (int i = 0; i < 8; i++){
		c[i] = new TCanvas(Form("c_%i", i), Form("c_%i", i), 1000, 1000);
	}

	c[0]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	LeadPt[0]->GetXaxis()->SetRangeUser(5, 30);
	LeadPt[0]->GetYaxis()->SetRangeUser(pow(10, -7), 1);
	SetAxisTitles(LeadPt[0], "p_{T} [GeV/#it{c}]", "arb. units");
	LeadPt[0]->Draw("EP");
	LeadPt[1]->Draw("EP SAME");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	c[0]->SaveAs("Plots/LeadPt.pdf");

	c[1]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	SubLeadPt[0]->GetXaxis()->SetRangeUser(5, 30);
	SubLeadPt[0]->GetYaxis()->SetRangeUser(pow(10, -7), 1);
	SetAxisTitles(SubLeadPt[0], "p_{T} [GeV/#it{c}]", "arb. units");
	SubLeadPt[0]->Draw("EP SAME");
	SubLeadPt[1]->Draw("EP SAME");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	c[1]->SaveAs("Plots/SubLeadPt.pdf");

	c[2]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	LeadJetLeadPt[0]->GetXaxis()->SetRangeUser(5, 30);
	LeadJetLeadPt[0]->GetYaxis()->SetRangeUser(pow(10, -7), 1);
	SetAxisTitles(LeadJetLeadPt[0], "p_{T} [GeV/#it{c}]", "arb. units");
	LeadJetLeadPt[0]->Draw("EP SAME");
	LeadJetLeadPt[1]->Draw("EP SAME");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	c[2]->SaveAs("Plots/LeadJetLeadPt.pdf");

	c[3]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	LeadJetSubLeadPt[0]->GetXaxis()->SetRangeUser(5, 30);
	LeadJetSubLeadPt[0]->GetYaxis()->SetRangeUser(pow(10, -7), 1);
	SetAxisTitles(LeadJetSubLeadPt[0], "p_{T} [GeV/#it{c}]", "arb. units");
	LeadJetSubLeadPt[0]->Draw("EP SAME");
	LeadJetSubLeadPt[1]->Draw("EP SAME");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	c[3]->SaveAs("Plots/LeadJetSubLeadPt.pdf");

	c[4]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	SubLeadJetLeadPt[0]->GetXaxis()->SetRangeUser(5, 30);
	SubLeadJetLeadPt[0]->GetYaxis()->SetRangeUser(pow(10, -7), 1);
	SetAxisTitles(SubLeadJetLeadPt[0], "p_{T} [GeV/#it{c}]", "arb. units");
	SubLeadJetLeadPt[0]->Draw("EP SAME");
	SubLeadJetLeadPt[1]->Draw("EP SAME");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	c[4]->SaveAs("Plots/SubLeadJetLeadPt.pdf");

	c[5]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	SubLeadJetSubLeadPt[0]->GetXaxis()->SetRangeUser(5, 30);
	SubLeadJetSubLeadPt[0]->GetYaxis()->SetRangeUser(pow(10, -7), 1);
	SetAxisTitles(SubLeadJetSubLeadPt[0], "p_{T} [GeV/#it{c}]", "arb. units");
	SubLeadJetSubLeadPt[0]->Draw("EP SAME");
	SubLeadJetSubLeadPt[1]->Draw("EP SAME");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	c[5]->SaveAs("Plots/SubLeadJetSubLeadPt.pdf");

	c[6]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	LeadJetPt[0]->GetXaxis()->SetRangeUser(5, 30);
	LeadJetPt[0]->GetYaxis()->SetRangeUser(pow(10, -7), 1);
	SetAxisTitles(LeadJetPt[0], "p_{T} [GeV/#it{c}]", "arb. units");
	LeadJetPt[0]->Draw("EP SAME");
	LeadJetPt[1]->Draw("EP SAME");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	c[6]->SaveAs("Plots/LeadJetPt.pdf");

	c[7]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	SubLeadJetPt[0]->GetXaxis()->SetRangeUser(5, 30);
	SubLeadJetPt[0]->GetYaxis()->SetRangeUser(pow(10, -7), 1);
	SetAxisTitles(SubLeadJetPt[0], "p_{T} [GeV/#it{c}]", "arb. units");
	SubLeadJetPt[0]->Draw("EP SAME");
	SubLeadJetPt[1]->Draw("EP SAME");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	c[7]->SaveAs("Plots/SubLeadJetPt.pdf");

	TCanvas *d[8];
	for (int i = 0; i < 8; i++){
		d[i] = new TCanvas(Form("d_%i", i), Form("d_%i", i), 1000, 1000);
	}

	TH1D *Ratio[8];

	d[0]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	Ratio[0] = (TH1D *)LeadPtW[1]->Clone();
	SetName(Ratio[0], "Lead p_{T}");
	Ratio[0]->Divide(LeadPtW[0]);
	Ratio[0]->GetXaxis()->SetRangeUser(5, 30);
	Ratio[0]->GetYaxis()->SetRangeUser(6*pow(10, -2), 5);
	SetAxisTitles(Ratio[0], "p_{T} [GeV/#it{c}]", "R_{CP}");
	Ratio[0]->Draw("EP");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	d[0]->SaveAs("Plots/LeadPt_Ratio.pdf");

	d[1]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	Ratio[1] = (TH1D *)SubLeadPtW[1]->Clone();
	SetName(Ratio[1], "Sub Lead p_{T}");
	Ratio[1]->Divide(SubLeadPtW[0]);
	Ratio[1]->GetXaxis()->SetRangeUser(5, 30);
	Ratio[1]->GetYaxis()->SetRangeUser(6*pow(10, -2), 5);
	SetAxisTitles(Ratio[1], "p_{T} [GeV/#it{c}]", "R_{CP}");
	Ratio[1]->Draw("EP");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	d[1]->SaveAs("Plots/SubLeadPt_Ratio.pdf");

	d[2]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	Ratio[2] = (TH1D *)LeadJetLeadPtW[1]->Clone();
	SetName(Ratio[2], "Lead Jet Lead p_{T}");
	Ratio[2]->Divide(LeadJetLeadPtW[0]);
	Ratio[2]->GetXaxis()->SetRangeUser(5, 30);
	Ratio[2]->GetYaxis()->SetRangeUser(6*pow(10, -2), 5);
	SetAxisTitles(Ratio[2], "p_{T} [GeV/#it{c}]", "R_{CP}");
	Ratio[2]->Draw("EP");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	d[2]->SaveAs("Plots/LeadJetLeadPt_Ratio.pdf");

	d[3]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	Ratio[3] = (TH1D *)LeadJetSubLeadPtW[1]->Clone();
	SetName(Ratio[3], "Lead Jet Sub Lead p_{T}");
	Ratio[3]->Divide(LeadJetSubLeadPtW[0]);
	Ratio[3]->GetXaxis()->SetRangeUser(5, 30);
	Ratio[3]->GetYaxis()->SetRangeUser(6*pow(10, -2), 5);
	SetAxisTitles(Ratio[3], "p_{T} [GeV/#it{c}]", "R_{CP}");
	Ratio[3]->Draw("EP");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	d[3]->SaveAs("Plots/LeadJetSubLeadPt_Ratio.pdf");

	d[4]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	Ratio[4] = (TH1D *)SubLeadJetLeadPtW[1]->Clone();
	SetName(Ratio[4], "Sublead Jet Lead p_{T}");
	Ratio[4]->Divide(SubLeadJetLeadPtW[0]);
	Ratio[4]->GetXaxis()->SetRangeUser(5, 30);
	Ratio[4]->GetYaxis()->SetRangeUser(6*pow(10, -2), 5);
	SetAxisTitles(Ratio[4], "p_{T} [GeV/#it{c}]", "R_{CP}");
	Ratio[4]->Draw("EP");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	d[4]->SaveAs("Plots/SubLeadJetLeadPt_Ratio.pdf");

	d[5]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	Ratio[5] = (TH1D *)SubLeadJetSubLeadPtW[1]->Clone();
	SetName(Ratio[5], "Sublead Jet Sub Lead p_{T}");
	Ratio[5]->Divide(SubLeadJetSubLeadPtW[0]);
	Ratio[5]->GetXaxis()->SetRangeUser(5, 30);
	Ratio[5]->GetYaxis()->SetRangeUser(6*pow(10, -2), 5);
	SetAxisTitles(Ratio[5], "p_{T} [GeV/#it{c}]", "R_{CP}");
	Ratio[5]->Draw("EP");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	d[5]->SaveAs("Plots/SubLeadJetSubLeadPt_Ratio.pdf");

	d[6]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	Ratio[6] = (TH1D *)LeadJetPtW[1]->Clone();
	SetName(Ratio[6], "Lead Jet p_{T}");
	Ratio[6]->Divide(LeadJetPtW[0]);
	Ratio[6]->GetXaxis()->SetRangeUser(5, 30);
	Ratio[6]->GetYaxis()->SetRangeUser(6*pow(10, -2), 5);
	SetAxisTitles(Ratio[6], "p_{T} [GeV/#it{c}]", "R_{CP}");
	Ratio[6]->Draw("EP");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	d[6]->SaveAs("Plots/LeadJetPt_Ratio.pdf");

	d[7]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	Ratio[7] = (TH1D *)SubLeadJetPtW[1]->Clone();
	SetName(Ratio[7], "Sub Lead Jet p_{T}");
	Ratio[7]->Divide(SubLeadJetPtW[0]);
	Ratio[7]->GetXaxis()->SetRangeUser(5, 30);
	Ratio[7]->GetYaxis()->SetRangeUser(6*pow(10, -2), 5);
	SetAxisTitles(Ratio[7], "p_{T} [GeV/#it{c}]", "R_{CP}");
	Ratio[7]->Draw("EP");
	gPad->BuildLegend(0.5,0.82,0.88,0.9);
	d[7]->SaveAs("Plots/SubLeadJetPt_Ratio.pdf");

}