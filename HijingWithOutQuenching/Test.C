#define Test_cxx
#include "Test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <algorithm>

void Test::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Test.C
//      root> Test t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   TH1D *LeadingParticlePt = new TH1D("LeadingParticlePt", "Leading Particle pT", 100, 0, 50);
   TH1D *SubLeadingParticlePt = new TH1D("SubLeadingParticlePt", "Sub Leading Particle pT", 100, 0, 50);

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      double mParticles_mPt[kMaxmParticles] = {};

      int goodparticles = 0;

      for (int i = 0; i < mParticles_; i++){
         if (mParticles_mStatus[i] <= 0) continue;

         TVector3 p(mParticles_mPx[i], mParticles_mPy[i], mParticles_mPz[i]);
         if (abs(p.Eta()) > 1.) continue;
         if (p.Pt() < 0.2) continue;
         mParticles_mPt[goodparticles] = p.Pt();

         // cout << p.Pt() << endl;
         goodparticles++;
      }

      sort(mParticles_mPt, mParticles_mPt + goodparticles, greater<double>());

      LeadingParticlePt->Fill(mParticles_mPt[0]);
      SubLeadingParticlePt->Fill(mParticles_mPt[1]);
   }

   TCanvas *c = new TCanvas("c", "c", 1000, 1000);
   gPad->SetLogy();
   LeadingParticlePt->Draw("EP");
   SubLeadingParticlePt->Draw("EP SAME");
   SubLeadingParticlePt->SetLineColor(kBlack);
   SubLeadingParticlePt->SetMarkerColor(kBlack);
}
