#define Test_cxx
#include "Test.h"

using namespace std;

// Include Fastjet Classes
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/Selector.hh>
#include "Bindef.h"

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

   fastjet::Strategy strategy = fastjet::Best;
   fastjet::RecombinationScheme recombScheme = fastjet::E_scheme; //Change as you need
   fastjet::JetDefinition *jetdefinition = NULL;
   fastjet::JetAlgorithm algorithm;

   algorithm = fastjet::antikt_algorithm; // Will use this one for now

   double R = 0.4;

   jetdefinition = new fastjet::JetDefinition(algorithm, R, recombScheme, strategy);

   // Defining Vectors for FastJet inputs
   std::vector <fastjet::PseudoJet> fjInputs; // Will store px, py, pz, E info

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   TH1D *LeadingJetPt = new TH1D("LeadingJetPt", "Leading Jet pT", 100, 0, 50);
   TH1D *SubLeadingJetPt = new TH1D("SubLeadingJetPt", "Subleading Jet pT", 100, 0, 50);

   TH1D *LeadingParticlePt = new TH1D("LeadingParticlePt", "Leading Particle pT", 100, 0, 50);
   TH1D *SubLeadingParticlePt = new TH1D("SubLeadingParticlePt", "Sub Leading Particle pT", 100, 0, 50);

   TH1D *LeadJetLeadParticle = new TH1D("LeadJetLeadParticle", "Lead Jet Leading Particle pT", 100, 0, 50);
   TH1D *SubLeadJetLeadParticle = new TH1D("SubLeadJetLeadParticle", "Sub Lead Jet Leading Particle pT", 100, 0, 50);

   TH1D *LeadJetSubLeadParticle = new TH1D("LeadJetSubLeadParticle", "Lead Jet Sub Leading Particle pT", 100, 0, 50);
   TH1D *SubLeadJetSubLeadParticle = new TH1D("SubLeadJetSubLeadParticle", "Sub Lead Jet Sub Leading Particle pT", 100, 0, 50);

   TH1D *LeadingJetPtW = new TH1D("LeadingJetPtW", "Leading Jet pT W", 100, 0, 50);
   TH1D *SubLeadingJetPtW = new TH1D("SubLeadingJetPtW", "Subleading Jet pT W", 100, 0, 50);

   TH1D *LeadingParticlePtW = new TH1D("LeadingParticlePtW", "Leading Particle pT W", 100, 0, 50);
   TH1D *SubLeadingParticlePtW = new TH1D("SubLeadingParticlePtW", "Sub Leading Particle pT W", 100, 0, 50);

   TH1D *LeadJetLeadParticleW = new TH1D("LeadJetLeadParticleW", "Lead Jet Leading Particle pT W", 100, 0, 50);
   TH1D *SubLeadJetLeadParticleW = new TH1D("SubLeadJetLeadParticleW", "Sub Lead Jet Leading Particle pT W", 100, 0, 50);

   TH1D *LeadJetSubLeadParticleW = new TH1D("LeadJetSubLeadParticleW", "Lead Jet Sub Leading Particle pT W", 100, 0, 50);
   TH1D *SubLeadJetSubLeadParticleW = new TH1D("SubLeadJetSubLeadParticleW", "Sub Lead Jet Sub Leading Particle pT W", 100, 0, 50);

   TH1D *LeadingJetPhotonCorr = new TH1D("LeadingJetPhotonCorr", "LeadingJetPhotonCorr", 100, -10, 10);
   TH1D *SubLeadingJetPhotonCorr = new TH1D("SubLeadingJetPhotonCorr", "SubLeadingJetPhotonCorr", 100, -10, 10);

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (jentry%1000==0) cout << jentry << " events finished." << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      double mParticles_mPt[kMaxmParticles] = {};

      int goodparticles = 0;

      fjInputs.clear();

      for (int i = 0; i < mParticles_; i++){
         TVector3 p(mParticles_mPx[i], mParticles_mPy[i], mParticles_mPz[i]);

         // if (mParticles_mId[i] == 22) cout << mParticles_mStatus[i] << "\t" << p.Pt() << "\t" << p.Eta() << "\t" << p.Phi() << endl;
         if (mParticles_mStatus[i] <= 0) continue;

         
         if (abs(p.Eta()) > 1.) continue;
         if (p.Pt() < 0.2) continue;
         mParticles_mPt[goodparticles] = p.Pt();

         if (p.Pt() < 2.0) continue;

         fastjet::PseudoJet k(mParticles_mPx[i], mParticles_mPy[i], mParticles_mPz[i], mParticles_mEnergy[i]);
         // Store as input to Fastjet
         fjInputs.push_back(k);   
         goodparticles++;
      }

      sort(mParticles_mPt, mParticles_mPt + goodparticles, greater<double>());

      double w = 1.0/(1.0*numberOfBinary);

      LeadingParticlePt->Fill(mParticles_mPt[0]);
      SubLeadingParticlePt->Fill(mParticles_mPt[1]);

      LeadingParticlePtW->Fill(mParticles_mPt[0], w);
      SubLeadingParticlePtW->Fill(mParticles_mPt[1], w);

      vector <fastjet::PseudoJet> inclusiveJets, sortedJets, selectedJets;
      fastjet::ClusterSequence clustSeq(fjInputs, *jetdefinition);

      // Extract Inclusive Jets sorted by pT
      inclusiveJets = clustSeq.inclusive_jets(10.0);
      sortedJets = sorted_by_pt(inclusiveJets);

      fastjet::Selector eta_selector = fastjet::SelectorEtaRange(-0.6, 0.6); // Selects |eta| < 0.6
      // fastjet::Selector pt_selector = fastjet::SelectorPtMin(10); // Selects jets with pT > 20 GeV/c
      fastjet::Selector selector = eta_selector;
      selectedJets = eta_selector(sortedJets);

      if (selectedJets.size() == 0) continue;
      // cout << "Selected Jets = " << selectedJets.size() << endl;

      for(unsigned jetid = 0; jetid < selectedJets.size(); jetid++){
         if (jetid > 1) continue;
         vector <fastjet::PseudoJet> constituents = selectedJets[jetid].constituents();
         vector <fastjet::PseudoJet> sortedconstituents = sorted_by_pt(constituents);

         if (sortedconstituents.size() <= 0) continue;

         // cout << "Jet # " << jetid << " Constituent Size = " << constituents.size() << "\t" << selectedJets[jetid].pt() << "\t" << sortedconstituents[0].pt() << "\t" << sortedconstituents[1].pt() << endl;

         if (jetid == 0) LeadingJetPt->Fill(selectedJets[jetid].pt());
         if (jetid == 1) SubLeadingJetPt->Fill(selectedJets[jetid].pt());

         if (jetid == 0) LeadingJetPtW->Fill(selectedJets[jetid].pt(), w);
         if (jetid == 1) SubLeadingJetPtW->Fill(selectedJets[jetid].pt(), w);

         if (jetid == 0) LeadJetLeadParticle->Fill(sortedconstituents[0].pt());
         if (jetid == 1) SubLeadJetLeadParticle->Fill(sortedconstituents[0].pt());

         if (jetid == 0) LeadJetLeadParticleW->Fill(sortedconstituents[0].pt(), w);
         if (jetid == 1) SubLeadJetLeadParticleW->Fill(sortedconstituents[0].pt(), w);

         if (jetid == 0 && sortedconstituents.size() > 1) LeadJetSubLeadParticle->Fill(sortedconstituents[1].pt());
         if (jetid == 1 && sortedconstituents.size() > 1) SubLeadJetSubLeadParticle->Fill(sortedconstituents[1].pt());

         if (jetid == 0 && sortedconstituents.size() > 1) LeadJetSubLeadParticleW->Fill(sortedconstituents[1].pt(), w);
         if (jetid == 1 && sortedconstituents.size() > 1) SubLeadJetSubLeadParticleW->Fill(sortedconstituents[1].pt(), w);
      }

      for (int i = 0; i < mParticles_; i++){
         TVector3 p(mParticles_mPx[i], mParticles_mPy[i], mParticles_mPz[i]);

         if (mParticles_mId[i] != 22) continue;
         if (mParticles_mStatus[i] <= 0) continue;
         if (abs(p.Eta()) > 1.) continue;
         // double phi1 = abs(selectedJets[0].phi() - standardPhi(p.Phi()));
         // if (selectedJets.size() > 1) phi2 = abs(selectedJets[1].phi() - standardPhi(p.Phi()));
         // LeadingJetPhotonCorr->Fill(selectedJets[0].phi() - standardPhi(p.Phi()));
         // if (selectedJets.size() > 1) SubLeadingJetPhotonCorr->Fill(selectedJets[1].phi() - standardPhi(p.Phi()));
         LeadingJetPhotonCorr->Fill(dPhi(selectedJets[0].phi(), standardPhi(p.Phi())));
         if (selectedJets.size() > 1) SubLeadingJetPhotonCorr->Fill(dPhi(selectedJets[1].phi(), standardPhi(p.Phi())));
      }
   }

   // TCanvas *c = new TCanvas("c", "c", 1000, 1000);
   // gPad->SetLogy();
   // LeadingParticlePt->Draw("EP");
   // SubLeadingParticlePt->Draw("EP SAME");
   // SubLeadingParticlePt->SetLineColor(kBlack);
   // SubLeadingParticlePt->SetMarkerColor(kBlack);

   TFile *f = new TFile("Output_NoQuenching_pthat15.root", "RECREATE");
   f->cd();

   LeadingParticlePt->Write();
   SubLeadingParticlePt->Write();
   LeadingJetPt->Write();
   SubLeadingJetPt->Write();
   LeadJetLeadParticle->Write();
   SubLeadJetLeadParticle->Write();
   LeadJetSubLeadParticle->Write();
   SubLeadJetSubLeadParticle->Write();

   LeadingParticlePtW->Write();
   SubLeadingParticlePtW->Write();
   LeadingJetPtW->Write();
   SubLeadingJetPtW->Write();
   LeadJetLeadParticleW->Write();
   SubLeadJetLeadParticleW->Write();
   LeadJetSubLeadParticleW->Write();
   SubLeadJetSubLeadParticleW->Write();

   LeadingJetPhotonCorr->Write();
   SubLeadingJetPhotonCorr->Write();
   f->Close();

}
