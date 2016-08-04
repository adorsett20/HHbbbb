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
#include <math.h>
#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentZVector>+;
#endif
#include "TROOT.h"
void TreeMaker(TString filename) {
  gROOT->ProcessLine("#include <vector>");
  TFile* file = TFile::Open(Form("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/%s", filename.Data()));
  TTree* Tree = (TTree*)file->Get("tree");
  TFile out(Form("/eos/uscms/store/user/adorsett/HHbbbb/BackTrees/%s",filename.Data()), "RECREATE");
  out.cd();
  Tree->SetBranchStatus("*",0);
  Tree->SetBranchStatus("RunNum",1);
  Tree->SetBranchStatus("BTags",1);
  Tree->SetBranchStatus("DeltaPhi1",1);
  Tree->SetBranchStatus("DeltaPhi2",1);
  Tree->SetBranchStatus("DeltaPhi3",1);
  Tree->SetBranchStatus("DeltaPhi4",1);
  Tree->SetBranchStatus("Electrons",1);
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
// Declare Variables
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
// Set Branch Addresses
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
  std::vector<double> CSV(0);
  std::vector<UInt_t> jets(0);
  std::vector<double>::iterator itC;
  std::vector<UInt_t>::iterator itJ;
  double DeltaMjj=0, Rmax=0, Mjj=0, r1=0,r2=0;
  double CSV1=0, CSV2=0, CSV3=0, CSV4=0;
  int config = 0;
  double HTnew=0,MHTnew=0,DeltaPhi1=0,DeltaPhi2=0,DeltaPhi3=0,DeltaPhi4=0;
  int NJetsnew=0, BTagsnew=0, MuonSize=0, ElectronSize=0;
  TLorentzVector pair1, pair2;
// Initialize output tree
  TTree *tree = new TTree("tree", "HHbbbb stats");
  tree->Branch("NJets20", &jetcheck);
  tree->Branch("DeltaMjj", &DeltaMjj);
  tree->Branch("Rmax", &Rmax);
  tree->Branch("AverageMjj", &Mjj);
  tree->Branch("HiggsCan1", &pair1);
  tree->Branch("HiggsCan2", &pair2);
  tree->Branch("CSV1", &CSV1);
  tree->Branch("CSV2", &CSV2);
  tree->Branch("CSV3", &CSV3);
  tree->Branch("CSV4", &CSV4); 
  tree->Branch("HT", &HTnew);
  tree->Branch("MHT", &MHTnew);
  tree->Branch("NJets", &NJetsnew);
  tree->Branch("BTags", &BTagsnew);
  tree->Branch("DeltaPhi1", &DeltaPhi1); 
  tree->Branch("DeltaPhi2", &DeltaPhi2);
  tree->Branch("DeltaPhi3", &DeltaPhi3);
  tree->Branch("DeltaPhi4", &DeltaPhi4);
  tree->Branch("RunNum", &NRun);
  tree->Branch("MuonSize", &MuonSize);
  tree->Branch("ElectronSize", &ElectronSize);
  tree->Branch("isoMuonTracks", &isoMuon);
  tree->Branch("isoElectronTracks", &isoElectron);
  tree->Branch("isoPionTracks", &isoPion);
  tree->Branch("Weight", &weight);
  tree->Branch("PFCaloMETRatio", &PFCaloMETRatio);
  tree->Branch("Jets", &Jets);
  tree->Branch("Muons", &Muons);
  tree->Branch("Electrons", &Electrons);
// Loop over Tree entries
  for(int i = 0; i < Tree->GetEntries(); i++) {
	Tree->GetEntry(i);
	jetcheck = 0;
	CSV.clear();
	jets.clear();
	for(int k = 0; k < 5; k++) {
		 CSV.push_back(0);
		 jets.push_back(0);  }
	for(UInt_t j = 0; j < Jets->size(); j++) {
                if(abs(Jets->at(j).Eta())<2.4 && Jets->at(j).Pt() > 20 && JetsID->at(j)) {
			itC = CSV.begin();
			itJ = jets.begin();
			if(bDisc->at(j) > CSV.at(0)) {
				CSV.insert(itC,bDisc->at(j));
				jets.insert(itJ,j); }
			else if(bDisc->at(j) > CSV.at(1)) {
				CSV.insert(itC+1,bDisc->at(j));
				jets.insert(itJ+1,j); }
			else if(bDisc->at(j) > CSV.at(2)) {
                                CSV.insert(itC+2,bDisc->at(j));
                                jets.insert(itJ+2,j); }
			else if(bDisc->at(j) > CSV.at(3)) {
                                CSV.insert(itC+3,bDisc->at(j));
                                jets.insert(itJ+3,j); } 
                        jetcheck++;  }
          	}        	
	CSV1 = CSV.at(0);
	CSV2 = CSV.at(1);
	CSV3 = CSV.at(2);
	CSV4 = CSV.at(3);
	if(CSV1 > 1) {CSV1 = 1;}
	if(CSV2 > 1) {CSV2 = 1;}
	if(CSV3 > 1) {CSV3 = 1;}
	if(CSV4 > 1) {CSV4 = 1;}
	config = 0;
	pair1 = Jets->at(jets.at(0)) + Jets->at(jets.at(1));
	pair2 = Jets->at(jets.at(2)) + Jets->at(jets.at(3));
	DeltaMjj = abs(pair1.M() - pair2.M());	
	pair1 = Jets->at(jets.at(0)) + Jets->at(jets.at(2));
        pair2 = Jets->at(jets.at(1)) + Jets->at(jets.at(3));
	if(abs(pair1.M() - pair2.M()) < DeltaMjj) {
		DeltaMjj = abs(pair1.M() - pair2.M());
		config = 1; }
	pair1 = Jets->at(jets.at(0)) + Jets->at(jets.at(3));
        pair2 = Jets->at(jets.at(1)) + Jets->at(jets.at(2));	
	if(abs(pair1.M() - pair2.M()) < DeltaMjj) {
		DeltaMjj = abs(pair1.M() - pair2.M());
		config = 2; }
	if(config == 0) {
		r1 = sqrt(pow(Jets->at(jets.at(0)).Eta() - Jets->at(jets.at(1)).Eta(),2) + pow(Jets->at(jets.at(0)).Phi() - Jets->at(jets.at(1)).Phi(),2));
		r2 = sqrt(pow(Jets->at(jets.at(2)).Eta() - Jets->at(jets.at(3)).Eta(),2) + pow(Jets->at(jets.at(2)).Phi() - Jets->at(jets.at(3)).Phi(),2));
		if(r1>r2) {Rmax = r1;}
		else{Rmax = r2;} 
		pair1 = Jets->at(jets.at(0)) + Jets->at(jets.at(1));
        	pair2 = Jets->at(jets.at(2)) + Jets->at(jets.at(3));
		Mjj = (pair1.M() + pair2.M())/2; }
	if(config == 1) {
                r1 = sqrt(pow(Jets->at(jets.at(0)).Eta() - Jets->at(jets.at(2)).Eta(),2) + pow(Jets->at(jets.at(0)).Phi() - Jets->at(jets.at(2)).Phi(),2));
                r2 = sqrt(pow(Jets->at(jets.at(1)).Eta() - Jets->at(jets.at(3)).Eta(),2) + pow(Jets->at(jets.at(1)).Phi() - Jets->at(jets.at(3)).Phi(),2));
                if(r1>r2) {Rmax = r1;}
                else{Rmax = r2;}
                pair1 = Jets->at(jets.at(0)) + Jets->at(jets.at(2));
                pair2 = Jets->at(jets.at(1)) + Jets->at(jets.at(3));
                Mjj = (pair1.M() + pair2.M())/2; }
	if(config == 2) {
                r1 = sqrt(pow(Jets->at(jets.at(0)).Eta() - Jets->at(jets.at(3)).Eta(),2) + pow(Jets->at(jets.at(0)).Phi() - Jets->at(jets.at(3)).Phi(),2));
                r2 = sqrt(pow(Jets->at(jets.at(1)).Eta() - Jets->at(jets.at(2)).Eta(),2) + pow(Jets->at(jets.at(1)).Phi() - Jets->at(jets.at(2)).Phi(),2));
                if(r1>r2) {Rmax = r1;}
                else{Rmax = r2;}
                pair1 = Jets->at(jets.at(0)) + Jets->at(jets.at(3));
                pair2 = Jets->at(jets.at(1)) + Jets->at(jets.at(2));
                Mjj = (pair1.M() + pair2.M())/2; }
	HTnew = HT;
	MHTnew = MHT;
	NJetsnew = NJets;
	BTagsnew = BTags;
	DeltaPhi1 = dPhi1;
	DeltaPhi2 = dPhi2;
	DeltaPhi3 = dPhi3;
	DeltaPhi4 = dPhi4;
	MuonSize = Muons->size();
	ElectronSize = Electrons->size();
	tree->Fill();
	}

  tree->Write();
  out.Close();

  return; }
