#include <TH1D.h>
#include <TCanvas.h>
#include "TTree.h"
#include <TFile.h>
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include <iostream>
#include <TVector.h>
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

double deltaR(double eta1, double phi1, double eta2, double phi2){
  return hypot(TVector2::Phi_mpi_pi(phi2-phi1), eta2-eta1);
}

void TreeMaker(TString filename) {
  gROOT->ProcessLine("#include <vector>");
//  TFile* file = TFile::Open(Form("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/%s", filename.Data()));
  TFile* file = TFile::Open("Total.root");
  TTree* Tree = (TTree*)file->Get("TreeMaker2/PreSelection");
//  TFile out(Form("/eos/uscms/store/user/adorsett/HHbbbb/BackTrees/%s",filename.Data()), "RECREATE");
  TFile out("Signal.root","RECREATE");
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
  Tree->SetBranchStatus("MET",1);
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
  Tree->SetBranchStatus("madHT",1);
  Tree->SetBranchStatus("GenEls",1);
  Tree->SetBranchStatus("GenMus",1);
  Tree->SetBranchStatus("GenTaus",1);
  Tree->SetBranchStatus("HBHENoiseFilter",1);
  Tree->SetBranchStatus("HBHEIsoNoiseFilter",1);
  Tree->SetBranchStatus("eeBadScFilter",1);
  Tree->SetBranchStatus("EcalDeadCellTriggerPrimitiveFilter",1);
  Tree->SetBranchStatus("NVtx",1);
  Tree->SetBranchStatus("LumiBlockNum",1);
  Tree->SetBranchStatus("TriggerPass", 1);
// Declare Variables
  UInt_t NRun=0;
  Int_t BTags=0, NJets=0;
  Double_t dPhi1=0, dPhi2=0, dPhi3=0, dPhi4=0, HT=0, MHT=0, weight=0, PFCaloMETRatio=0;
  Double_t madHT=0, MET=0;
  Int_t NVtx=0;
  UInt_t Lumi=0;
  Int_t isoMuon=0, isoElectron=0, isoPion=0;
  Bool_t JetID=0, noMuonJet=0;
  Int_t eeBad=0,EcalDead=0;
  Int_t HBHE=0,HBHEiso=0;
  std::vector<int> *Trigger(0); TBranch *bTrigger = nullptr;
  std::vector<TLorentzVector> *Electrons(0); TBranch *bElectrons = nullptr;
  std::vector<TLorentzVector> *Muons(0); TBranch *bMuons = nullptr;
  std::vector<TLorentzVector> *Jets(0); TBranch *bJets = nullptr;
  std::vector<TLorentzVector> *GenEls(0); TBranch *bGenEls = nullptr;
  std::vector<TLorentzVector> *GenMus(0); TBranch *bGenMus = nullptr;
  std::vector<TLorentzVector> *GenTaus(0); TBranch *bGenTaus = nullptr;
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
  Tree->SetBranchAddress("MET", &MET);
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
  Tree->SetBranchAddress("madHT", &madHT);
  Tree->SetBranchAddress("GenEls", &GenEls, &bGenEls);
  Tree->SetBranchAddress("GenMus", &GenMus, &bGenMus);
  Tree->SetBranchAddress("GenTaus", &GenTaus, &bGenTaus);
  Tree->SetBranchAddress("HBHENoiseFilter", &HBHE);
  Tree->SetBranchAddress("HBHEIsoNoiseFilter", &HBHEiso);
  Tree->SetBranchAddress("eeBadScFilter", &eeBad);
  Tree->SetBranchAddress("EcalDeadCellTriggerPrimitiveFilter", &EcalDead);
  Tree->SetBranchAddress("NVtx", &NVtx);
  Tree->SetBranchAddress("LumiBlockNum", &Lumi);
  Tree->SetBranchAddress("TriggerPass", &Trigger, &bTrigger);
  int jetcheck = 0, jet20=0;
  int BTagcheck = 0;
  std::vector<double> CSV(0);
  std::vector<UInt_t> jets(0);
  std::vector<double>::iterator itC;
  std::vector<UInt_t>::iterator itJ;
  double DeltaMjj=0, Rmax=0, Mjj=0, r1=0,r2=0;
  double CSV1=0, CSV2=0, CSV3=0, CSV4=0;
  int config = 0;
  int MuonSize=0, ElectronSize=0;
  Int_t Bloose=0, Bmedium=0, Btight=0;
  TLorentzVector pair1, pair2, empty(0.,0.,0.,0.);
// Initialize output tree
  TTree *tree = new TTree("tree", "HHbbbb stats");
  tree->Branch("NJets20", &jet20);
  tree->Branch("DeltaMjj", &DeltaMjj);
  tree->Branch("Rmax", &Rmax);
  tree->Branch("AverageMjj", &Mjj);
  tree->Branch("HiggsCan1", &pair1);
  tree->Branch("HiggsCan2", &pair2);
  tree->Branch("CSV1", &CSV1);
  tree->Branch("CSV2", &CSV2);
  tree->Branch("CSV3", &CSV3);
  tree->Branch("CSV4", &CSV4); 
  tree->Branch("HT", &HT);
  tree->Branch("MHT", &MHT);
  tree->Branch("MET", &MET);
  tree->Branch("NJets", &NJets);
  tree->Branch("BTags", &BTags);
  tree->Branch("DeltaPhi1", &dPhi1); 
  tree->Branch("DeltaPhi2", &dPhi2);
  tree->Branch("DeltaPhi3", &dPhi3);
  tree->Branch("DeltaPhi4", &dPhi4);
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
  tree->Branch("madHT", &madHT);
  tree->Branch("GenEls", &GenEls);
  tree->Branch("GenMus", &GenMus);
  tree->Branch("GenTaus",&GenTaus);
  tree->Branch("HBHENoiseFilter", &HBHE);
  tree->Branch("HBHEIsoNoiseFilter", &HBHEiso);
  tree->Branch("eeBadScFilter", &eeBad);
  tree->Branch("EcalDeadCellTriggerPrimitiveFilter", &EcalDead);
  tree->Branch("NVtx", &NVtx);
  tree->Branch("LumiBlockNum", &Lumi);
  tree->Branch("Jets_ID", &JetsID);
  tree->Branch("Jets_bDiscriminatorCSV", &bDisc);
  tree->Branch("Bloose", &Bloose);
  tree->Branch("Bmedium", &Bmedium);
  tree->Branch("Btight", &Btight);
// Loop over Tree entries
  for(int i = 0; i < Tree->GetEntries(); i++) {
	Tree->GetEntry(i);
	jetcheck = 0;
	jet20 = 0;
	Bloose = 0;
	Bmedium = 0;
	Btight = 0;
	CSV.clear();
	jets.clear();
	DeltaMjj = -1;
	Rmax = -1;
	Mjj = -1;
	pair1 = empty;
	pair2 = empty;
	if(Trigger->at(40)||Trigger->at(41)||Trigger->at(44)||Trigger->at(45)){
	for(int k = 0; k < 5; k++) {
		 CSV.push_back(0);
		 jets.push_back(0);  }
	for(UInt_t j = 0; j < Jets->size(); j++) {
                if(abs(Jets->at(j).Eta())<2.4 && Jets->at(j).Pt() > 30 && JetsID->at(j)) {
			if(Jets->at(j).Pt() > 20) jet20++;
			if(bDisc->at(j) > .970) Btight++;
			if(bDisc->at(j) > .890) Bmedium++;
			if(bDisc->at(j) > .605) Bloose++;
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
	if(NJets >= 4) {
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
		r1 = deltaR(Jets->at(jets.at(0)).Eta(), Jets->at(jets.at(0)).Phi(), Jets->at(jets.at(1)).Eta(), Jets->at(jets.at(1)).Phi());
		r2 = deltaR(Jets->at(jets.at(2)).Eta(), Jets->at(jets.at(2)).Phi(), Jets->at(jets.at(3)).Eta(), Jets->at(jets.at(3)).Phi());
		if(r1>r2) {Rmax = r1;}
		else{Rmax = r2;} 
		pair1 = Jets->at(jets.at(0)) + Jets->at(jets.at(1));
        	pair2 = Jets->at(jets.at(2)) + Jets->at(jets.at(3));
		Mjj = (pair1.M() + pair2.M())/2; }
	if(config == 1) {
		r1 = deltaR(Jets->at(jets.at(0)).Eta(), Jets->at(jets.at(0)).Phi(), Jets->at(jets.at(2)).Eta(), Jets->at(jets.at(2)).Phi());
                r2 = deltaR(Jets->at(jets.at(1)).Eta(), Jets->at(jets.at(1)).Phi(), Jets->at(jets.at(3)).Eta(), Jets->at(jets.at(3)).Phi());
                if(r1>r2) {Rmax = r1;}
                else{Rmax = r2;}
                pair1 = Jets->at(jets.at(0)) + Jets->at(jets.at(2));
                pair2 = Jets->at(jets.at(1)) + Jets->at(jets.at(3));
                Mjj = (pair1.M() + pair2.M())/2; }
	if(config == 2) {
		r1 = deltaR(Jets->at(jets.at(0)).Eta(), Jets->at(jets.at(0)).Phi(), Jets->at(jets.at(3)).Eta(), Jets->at(jets.at(3)).Phi());
                r2 = deltaR(Jets->at(jets.at(1)).Eta(), Jets->at(jets.at(1)).Phi(), Jets->at(jets.at(2)).Eta(), Jets->at(jets.at(2)).Phi());
                if(r1>r2) {Rmax = r1;}
                else{Rmax = r2;}
                pair1 = Jets->at(jets.at(0)) + Jets->at(jets.at(3));
                pair2 = Jets->at(jets.at(1)) + Jets->at(jets.at(2));
                Mjj = (pair1.M() + pair2.M())/2; }
	}
	MuonSize = Muons->size();
	ElectronSize = Electrons->size();
	tree->Fill();
	} }

  tree->Write();
  out.Close();

  return; }




void BigTree(){

  TreeMaker("WORDS");

/*
  TreeMaker("tree_ZZTo2Q2Nu.root");
  TreeMaker("tree_ZZTo2L2Q.root");
  TreeMaker("tree_ZJetsToNuNu_HT-100to200.root");
  TreeMaker("tree_ZJetsToNuNu_HT-200to400.root");
  TreeMaker("tree_ZJetsToNuNu_HT-400to600.root");
  TreeMaker("tree_ZJetsToNuNu_HT-600to800.root");
  TreeMaker("tree_ZJetsToNuNu_HT-800to1200.root");
  TreeMaker("tree_ZJetsToNuNu_HT-1200to2500.root");
  TreeMaker("tree_ZJetsToNuNu_HT-2500toInf.root");
  TreeMaker("tree_WZZ.root");
  TreeMaker("tree_WZTo1L3Nu.root");
  TreeMaker("tree_WZTo1L1Nu2Q.root");
  TreeMaker("tree_WWZ.root");
  TreeMaker("tree_WWTo2L2Nu.root");
  TreeMaker("tree_WWTo1L1Nu2Q.root");
  TreeMaker("tree_WJetsToLNu_HT-100to200.root");
  TreeMaker("tree_WJetsToLNu_HT-200to400.root");
  TreeMaker("tree_WJetsToLNu_HT-400to600.root");
  TreeMaker("tree_WJetsToLNu_HT-600to800.root");
  TreeMaker("tree_WJetsToLNu_HT-800to1200.root");
  TreeMaker("tree_WJetsToLNu_HT-1200to2500.root");
  TreeMaker("tree_WJetsToLNu_HT-2500toInf.root");
  TreeMaker("tree_TTZToQQ.root");
  TreeMaker("tree_TTZToLLNuNu.root");
  TreeMaker("tree_TTWJetsToQQ.root");
  TreeMaker("tree_TTWJetsToLNu.root");
  TreeMaker("tree_TTTT.root");
  TreeMaker("tree_TTJets_SingleLeptFromT.root");
  TreeMaker("tree_TTJets_SingleLeptFromTbar.root");
  TreeMaker("tree_TTJets.root");
  TreeMaker("tree_TTJets_HT-600to800.root");
  TreeMaker("tree_TTJets_HT-800to1200.root");
  TreeMaker("tree_TTJets_HT-1200to2500.root");
  TreeMaker("tree_TTJets_HT-2500toInf.root");
  TreeMaker("tree_TTJets_DiLept.root");
  TreeMaker("tree_ST_tW_top.root");
  TreeMaker("tree_ST_tW_antitop.root");
  TreeMaker("tree_ST_t-channel_top.root");
  TreeMaker("tree_ST_t-channel_antitop.root");
  TreeMaker("tree_ST_s-channel.root");
  TreeMaker("tree_QCD_HT-200to300.root");
  TreeMaker("tree_QCD_HT-300to500.root");
  TreeMaker("tree_QCD_HT-500to700.root");
  TreeMaker("tree_QCD_HT-700to1000.root");
  TreeMaker("tree_QCD_HT-1000to1500.root");
  TreeMaker("tree_QCD_HT-1500to2000.root");
  TreeMaker("tree_QCD_HT-2000toInf.root");
  TreeMaker("tree_DYJetsToLL_M-50_HT-600toInf.root");
  TreeMaker("tree_DYJetsToLL_M-50_HT-400to600.root");
  TreeMaker("tree_DYJetsToLL_M-50_HT-200to400.root");
  TreeMaker("tree_DYJetsToLL_M-50_HT-100to200.root");
*/
  return; }



