#include <TH1D.h>
#include <TCanvas.h>
#include "TTree.h"
#include <TFile.h>
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include <iostream>
#include <TH2D.h>
#include <stdio.h>
#include <TObject.h>
#include <vector>
#include <TCut.h>
#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentZVector>+;
#endif
#include "TROOT.h"
void Efficiency() {
  gROOT->ProcessLine("#include <vector>");
  TFile* file = TFile::Open("Total.root");
  TTree* Tree = (TTree*)file->Get("TreeMaker2/PreSelection");
  Tree->SetBranchStatus("*",0);
  Tree->SetBranchStatus("RunNum",1);
  Tree->SetBranchStatus("BTags",1);
  Tree->SetBranchStatus("DeltaPhi1",1);
  Tree->SetBranchStatus("DeltaPhi2",1);
  Tree->SetBranchStatus("DeltaPhi3",1);
  Tree->SetBranchStatus("DeltaPhi4",1);
  Tree->SetBranchStatus("HT",1);
  Tree->SetBranchStatus("MHT",1);
  Tree->SetBranchStatus("NJets",1);
  Tree->SetBranchStatus("isoMuonTracks",1);
  Tree->SetBranchStatus("isoElectronTracks",1);
  Tree->SetBranchStatus("isoPionTracks",1);
  Tree->SetBranchStatus("Muons",1);
  Tree->SetBranchStatus("Electrons",1);
  Tree->SetBranchStatus("Jets",1);
  Tree->SetBranchStatus("Weight",1);
  Tree->SetBranchStatus("JetID",1);
  Tree->SetBranchStatus("PFCaloMETRatio",1);
  Tree->SetBranchStatus("Jets_bDiscriminatorCSV",1);  
  Tree->SetBranchStatus("Jets_ID",1);
//Make histogram to store the events that pass cuts
  TH1D* skim = new TH1D("skim", "Passing Skim Events", 12, -.5, 11.5);
  TH1D* skimNo = new TH1D("skimNo", "Passing Skim Events, No Order", 12, -.5, 11.5);
//Declare variables
  UInt_t NRun=0;
  Int_t BTags=0, NJets=0;
  Double_t dPhi1=0, dPhi2=0, dPhi3=0, dPhi4=0, HT=0, MHT=0, weight=0, PFCaloMETRatio=0;
  Int_t isoMuon=0, isoElectron=0, isoPion=0;
  Bool_t JetID=0, noMuonJet=0;
  std::vector<TLorentzVector> *Electrons(0); TBranch *bElectrons = nullptr;
  std::vector<TLorentzVector> *Muons(0); TBranch *bMuons = nullptr;
  std::vector<TLorentzVector> *Jets(0); TBranch *bJets = nullptr;
  std::vector<double> *bDisc(0); TBranch *bbDisc = nullptr;
  std::vector<Bool_t> *JetsID(0); TBranch *bJetsID = nullptr;
//Set Branch Addresses
  Tree->SetBranchAddress("RunNum",&NRun);
  Tree->SetBranchAddress("NJets", &NJets);
  Tree->SetBranchAddress("BTags", &BTags);
  Tree->SetBranchAddress("DeltaPhi1", &dPhi1);
  Tree->SetBranchAddress("DeltaPhi2", &dPhi2);
  Tree->SetBranchAddress("DeltaPhi3", &dPhi3);
  Tree->SetBranchAddress("DeltaPhi4", &dPhi4);
  Tree->SetBranchAddress("HT", &HT);
  Tree->SetBranchAddress("MHT", &MHT);
  Tree->SetBranchAddress("Weight",&weight);
  Tree->SetBranchAddress("isoMuonTracks", &isoMuon);
  Tree->SetBranchAddress("isoElectronTracks", &isoElectron);
  Tree->SetBranchAddress("isoPionTracks", &isoPion);
  Tree->SetBranchAddress("JetID", &JetID);
  Tree->SetBranchAddress("Electrons", &Electrons, &bElectrons);
  Tree->SetBranchAddress("Muons", &Muons, &bMuons);
  Tree->SetBranchAddress("Jets", &Jets, &bJets);
  Tree->SetBranchAddress("PFCaloMETRatio", &PFCaloMETRatio);
  Tree->SetBranchAddress("Jets_bDiscriminatorCSV", &bDisc, &bbDisc);
  Tree->SetBranchAddress("Jets_ID", &JetsID, &bJetsID);
  int jetcheck = 0;
  int BTagcheck = 0;
//Loop through events and check cuts
  for(int i = 0; i < Tree->GetEntries(); i++) {
        Tree->GetEntry(i);
        noMuonJet = true;
        skim->Fill(11,weight);
        skimNo->Fill(11,weight);
	jetcheck = 0;
	BTagcheck = 0;
	for(UInt_t j = 0; j < Jets->size(); j++) {
		if(abs(Jets->at(j).Eta())<2.4 && Jets->at(j).Pt() > 20 && JetsID->at(j)) {
			jetcheck++;
			if(bDisc->at(j) > .898) {
				BTagcheck++; }
			}
		}
        if(NRun <= 274443) {skimNo->Fill(0.0,weight);}
        if(jetcheck==4||jetcheck==5) {skimNo->Fill(1,weight);}
        if(HT>200) {skimNo->Fill(2,weight);}
        if(MHT>200) {skimNo->Fill(3,weight);}
        if(Muons->size()==0) {skimNo->Fill(4,weight);}
        if(Electrons->size()==0) {skimNo->Fill(5,weight);}
        if(isoMuon==0) {skimNo->Fill(6,weight);}
        if(isoElectron==0) {skimNo->Fill(7,weight);}
        if(isoPion==0) {skimNo->Fill(8,weight);}
        if(dPhi1 > .5 && dPhi2 > .5 && dPhi3 > .3 && dPhi4 > .3) {skimNo->Fill(9,weight);}
        if(JetID && PFCaloMETRatio < 5 && noMuonJet) {skimNo->Fill(10,weight);}

        if(NRun <= 274443) {
         skim->Fill(0.0,weight);
         if(jetcheck==4||jetcheck==5) {
          skim->Fill(1,weight);
          if(HT>200) {
           skim->Fill(2,weight);
           if(MHT>200) {
            skim->Fill(3,weight);
            if(Muons->size()==0) {
             skim->Fill(4,weight);
             if(Electrons->size()==0) {
              skim->Fill(5,weight);
              if(isoMuon==0) {
               skim->Fill(6,weight);
               if(isoElectron==0) {
                skim->Fill(7,weight);
                if(isoPion==0) {
                 skim->Fill(8,weight);
                 if(dPhi1 > .5 && dPhi2 > .5 && dPhi3 > .3 && dPhi4 > .3) {
                  skim->Fill(9,weight);
                  if(JetID && PFCaloMETRatio < 5 && noMuonJet) {
                   skim->Fill(10,weight);
                        } } } } } } } } } } }

        }
  cout << "Efficiency Table" << endl;
  cout << "Selector| Cut Details | Efficiency 1 | Efficiency 2" << endl;
  cout << "   0    |RunNum<=27443|  " << skim->GetBinContent(1)/skim->GetBinContent(12) << "  |  " << skimNo->GetBinContent(1)/skimNo->GetBinContent(12) << endl;
  cout << "   1A   |NJets == 4||5|  " << skim->GetBinContent(2)/skim->GetBinContent(1) << "  |  " << skimNo->GetBinContent(2)/skimNo->GetBinContent(12) << endl;
  cout << "   2A   |HT > 200     |  " << skim->GetBinContent(3)/skim->GetBinContent(2) << "  |  " << skimNo->GetBinContent(3)/skimNo->GetBinContent(12) << endl;
  cout << "   3A   |MHT > 200    |  " << skim->GetBinContent(4)/skim->GetBinContent(3) << "  |  " << skimNo->GetBinContent(4)/skimNo->GetBinContent(12) << endl;
  cout << "   4A   |MuonSize == 0|  " << skim->GetBinContent(5)/skim->GetBinContent(4) << "  |  " << skimNo->GetBinContent(5)/skimNo->GetBinContent(12) << endl;
  cout << "   5A   |ElecSize == 0|  " << skim->GetBinContent(6)/skim->GetBinContent(5) << "  |  " << skimNo->GetBinContent(6)/skimNo->GetBinContent(12) << endl;
  cout << "   7    |IsoMuon == 0 |  " << skim->GetBinContent(7)/skim->GetBinContent(6) << "  |  " << skimNo->GetBinContent(7)/skimNo->GetBinContent(12) << endl;
  cout << "   8    |IsoElec == 0 |  " << skim->GetBinContent(8)/skim->GetBinContent(7) << "  |  " << skimNo->GetBinContent(8)/skimNo->GetBinContent(12) << endl;  
  cout << "   9    |IsoPion == 0 |  " << skim->GetBinContent(9)/skim->GetBinContent(8) << "  |  " << skimNo->GetBinContent(9)/skimNo->GetBinContent(12) << endl;
  cout << "   10A  |Delta Phi Cut|  " << skim->GetBinContent(10)/skim->GetBinContent(9) << "  |  " << skimNo->GetBinContent(10)/skimNo->GetBinContent(12) << endl;
  cout << "   11B  |Clean Cuts   |  " << skim->GetBinContent(11)/skim->GetBinContent(10) << "  |  " << skimNo->GetBinContent(11)/skimNo->GetBinContent(12) << endl;

  TFile out("Skim.root", "RECREATE");
  skim->Write();
  skimNo->Write();
  out.Close();
  return; }





