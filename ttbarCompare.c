#include <iostream>
#include <vector>
#include <math.h>
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TCut.h"
#include "THStack.h"
#include "TLine.h"
#include "TFileCollection.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "jack_style.h"
#include <string>
using namespace std;

TString plotdir = "StackedPlots/";

void Plot(TString var, TCut OtherCuts, TCut trig, TCut weight, TH1D* temp, TString skimdir="tree_signal", TString data_skimdir="tree_signalUnblind",
	  TString options="plotSig:plotLog", TString title="default") {
//Boolians for plot drawing options
  bool plotLog = options.Contains("plotLog") && (!options.Contains("!plotLog"));
  bool entries = options.Contains("entries") && (!options.Contains("!entries"));
//TChain setup
  TChain* ttbar = new TChain("tree");
  ttbar->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTJets_HT-2500toInf.root");
  ttbar->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTJets_HT-1200to2500.root");
  ttbar->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTJets_HT-600to800.root");
  ttbar->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTJets_HT-800to1200.root");
  TChain* ttbar_lep = new TChain("tree");
  ttbar_lep->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTJets_SingleLeptFromT.root");
  ttbar_lep->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTJets_SingleLeptFromTbar.root");
  ttbar_lep->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTJets_DiLept.root");
//Make output file
  TFile* outfile = new TFile(Form("hists/%s.root", title.Data()), "recreate"); 
//Make Histograms from input temp
  TH1D* httbar2 = (TH1D*)temp->Clone("ttbar2");
  TH1D* httbar3 = (TH1D*)temp->Clone("ttbar3");
  TH1D* httbar4 = (TH1D*)temp->Clone("ttbar4");
//Make Historgrams with SIG/SBD cuts
  TH1D* sig2B = (TH1D*)temp->Clone("sig2B");
  TH1D* sbd2B = (TH1D*)temp->Clone("sbd2B");
  TH1D* sig3B = (TH1D*)temp->Clone("sig3B");
  TH1D* sbd3B = (TH1D*)temp->Clone("sbd3B");
  TH1D* sig4B = (TH1D*)temp->Clone("sig4B");
  TH1D* sbd4B = (TH1D*)temp->Clone("sbd4B");
//Additional Cuts
  TCut selection(OtherCuts);
  TCut tempo("@GenEls.size()+@GenMus.size() == 0 && MHT > 200 && CSV2 > .898 && NJets <= 5 && DeltaPhi1>0.5&&DeltaPhi2>0.5&&DeltaPhi3>0.3&&DeltaPhi4>0.3 && isoMuonTracks+isoElectronTracks+isoPionTracks==0 && NJets >=4 && Rmax < 2.2");
//Sig/Back Cuts
  TCut SIG("DeltaMjj < 40 && AverageMjj > 100 && AverageMjj < 140");
  TCut SBD("DeltaMjj > 40 || AverageMjj < 100 || AverageMjj > 140");
//4b Cuts
  TCut backB4("CSV2 > .935 && CSV3 > .8 && CSV4 > .460");
//3b Cuts
  TCut backB3("CSV2 > .935 && CSV3 > .8 && CSV4 < .460");
//2b Cuts
  TCut backB2("CSV2 > .935 && CSV3 < .8");
//Fill histograms
  ttbar->Project("ttbar2", var, (selection+backB2)*weight);
  ttbar->Project("ttbar3", var, (selection+backB3)*weight);
  ttbar->Project("ttbar4", var, (selection+backB4)*weight);
  ttbar_lep->Project("+ttbar2", var, ((selection+backB2)+"madHT<600")*weight);
  ttbar_lep->Project("+ttbar3", var, ((selection+backB3)+"madHT<600")*weight);
  ttbar_lep->Project("+ttbar4", var, ((selection+backB4)+"madHT<600")*weight);
//Fil Signal/Background histograms
  ttbar->Project("sig2B", var, (SIG+selection+backB2)*weight);
  ttbar->Project("sig3B", var, (SIG+selection+backB3)*weight);
  ttbar->Project("sig4B", var, (SIG+selection+backB4)*weight);
  ttbar_lep->Project("+sig2B", var, ((SIG+selection+backB2)+"madHT<600")*weight);
  ttbar_lep->Project("+sig3B", var, ((SIG+selection+backB3)+"madHT<600")*weight);
  ttbar_lep->Project("+sig4B", var, ((SIG+selection+backB4)+"madHT<600")*weight);
  ttbar->Project("sbd2B", var, (SBD+selection+backB2)*weight);
  ttbar->Project("sbd3B", var, (SBD+selection+backB3)*weight);
  ttbar->Project("sbd4B", var, (SBD+selection+backB4)*weight);
  ttbar_lep->Project("+sbd2B", var, ((SBD+selection+backB2)+"madHT<600")*weight);
  ttbar_lep->Project("+sbd3B", var, ((SBD+selection+backB3)+"madHT<600")*weight);
  ttbar_lep->Project("+sbd4B", var, ((SBD+selection+backB4)+"madHT<600")*weight);
  cout << "SIG/SBD ratios" << endl;
  cout << "  2B: " << sig2B->Integral()/sbd2B->Integral() << endl;
  cout << "  3B: " << sig3B->Integral()/sbd3B->Integral() << endl;
  cout << "  4B: " << sig4B->Integral()/sbd4B->Integral() << endl;
  cout << "Kappa factors" << endl;
  cout << "  3B/2B: " << (sig3B->Integral()/sbd3B->Integral())/(sig2B->Integral()/sbd2B->Integral()) << endl;
  cout << "  4B/2B: " << (sig4B->Integral()/sbd4B->Integral())/(sig2B->Integral()/sbd2B->Integral()) << endl;
//Fix overflow by adding to last bin in hist
  bool addOverflow(true);
  Double_t e_overflow(0.), i_overflow(0.);
  int nbinsx=temp->GetNbinsX();
  if (addOverflow) {
    i_overflow=httbar2->IntegralAndError(nbinsx,nbinsx+1,e_overflow);
    httbar2->SetBinContent(nbinsx, i_overflow);
    httbar2->SetBinError(nbinsx, e_overflow);
    i_overflow=httbar3->IntegralAndError(nbinsx,nbinsx+1,e_overflow);
    httbar3->SetBinContent(nbinsx, i_overflow);
    httbar3->SetBinError(nbinsx, e_overflow);
    i_overflow=httbar4->IntegralAndError(nbinsx,nbinsx+1,e_overflow);
    httbar4->SetBinContent(nbinsx, i_overflow);
    httbar4->SetBinError(nbinsx, e_overflow);
  }
//Use set_style function on histograms
  set_style_lite(httbar2, "ttbar2");
  httbar2->SetLineColor(kGreen);
  set_style_lite(httbar3, "ttbar3");
  httbar3->SetLineColor(kAzure);
  set_style_lite(httbar4, "ttbar4");
  httbar4->SetLineColor(kRed);
//Setup histogram labels
  TH1D* hists[5] = {httbar2, httbar3, httbar4};
  TString hlabels[5] = {"t#bar{t} N_{b} = 2", "t#bar{t} N_{b} = 3", "t#bar{t} N_{b} = 4"};
//Write histograms to output file
  outfile->cd();
  httbar2->Write();
  httbar3->Write();
  httbar3->Write();
//Setup Legends
  TLegend * leg2 = new TLegend(0.72, 0.74, 0.94, 0.92);
  set_style(leg2,0.035);
  leg2->AddEntry(hists[0], hlabels[0].Data(), "l");
  leg2->AddEntry(hists[1], hlabels[1].Data(), "l");
  leg2->AddEntry(hists[2], hlabels[2].Data(), "l");
  double ymax = httbar4->GetMaximum();
  if (httbar4->GetMaximum()>ymax) ymax=httbar4->GetMaximum();
  if(plotLog) {
    httbar4->SetMaximum(1000*ymax);
    httbar4->SetMinimum(0.08);
    if (var.Contains("num_")) httbar4->SetMaximum(1000*ymax);
    if (var.Contains("eta")) httbar4->SetMaximum(1000*ymax);
    if (var.Contains("dphi")) httbar4->SetMaximum(1000*ymax);
  }
  else {
    httbar4->SetMinimum(0);
    httbar4->SetMaximum(1.75*ymax);
    if (title.Contains("cdtt")) httbar4->SetMaximum(1.75*ymax);
    if (title.Contains("/")) httbar4->SetMaximum(1.75*ymax);
    if (var.Contains("eta")) httbar4->SetMaximum(3*ymax);
  }
//Normalize Histograms
  httbar2->Scale(100./httbar2->Integral());
  httbar3->Scale(100./httbar3->Integral());
  httbar4->Scale(100./httbar4->Integral());
  httbar4->GetYaxis()->SetTitle("% events");
  ymax = httbar4->GetMaximum();
  httbar4->SetMaximum(1.8*ymax);
  cout << "Histogram Averages" << endl;
  cout << "  2B:" << httbar2->GetMean() << endl;
  cout << "  3B:" << httbar3->GetMean() << endl;
  cout << "  4B:" << httbar4->GetMean() << endl;
//Canvas and Pad Setup
  TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
 // c1->SetLogy(plotLog);
  httbar4->SetMarkerColor(kRed);
  httbar4->SetMarkerSize(.3);
  httbar4->Draw("e1");
  httbar4->GetXaxis()->SetLabelSize(0.03);
  httbar4->GetXaxis()->SetTitleSize(0.025);
  httbar4->GetYaxis()->SetLabelSize(0.04);
//Draw Histograms
  httbar3->SetMarkerColor(kAzure);
  httbar3->SetMarkerSize(.3);
  httbar3->Draw("e1same");
  httbar2->Draw("histsame");
//Draw Legends
  leg2->Draw();
  TLatex * latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextFont(62);
  latex->SetTextSize(0.032);
  latex->DrawLatex(0.19, 0.89, "CMS Simulation");
  latex->SetTextSize(0.022);
  latex->DrawLatex(0.19, 0.86, "#sqrt{s} = 13 TeV, L = 30 fb^{-1}");
  c1->cd();
  if (title.EqualTo("default")) title=plotdir+var;
  gPad->Print(plotdir+title+".pdf");   
//Clean up histograms
  delete ttbar;
  delete leg2;
  delete latex;
//  delete pad1;
  delete c1;

  outfile->Close();
  return; }

void ttbarCompare() { 
  TH1::SetDefaultSumw2(1);
  if (gSystem->AccessPathName(plotdir))
    gSystem->mkdir(plotdir);
//Style Setup
  cout << "Setting tdr style."  << endl;
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  setTDRStyle(tdrStyle);
  tdrStyle->cd();
//Put in Cuts
  TCut baseline("MHT>250&&HT>250");
  TCut ZL("@Muons.size()+@Electrons.size()==0 && isoMuonTracks+isoElectronTracks+isoPionTracks==0");
  TCut dPhi("DeltaPhi1>0.5&&DeltaPhi2>0.5&&DeltaPhi3>0.3&&DeltaPhi4>0.3");

  TCut antiPhi("DeltaPhi1>2&&DeltaPhi2>2&&DeltaPhi3>1.8&&DeltaPhi4>1.8");

  TCut JetFix("Sum$(abs(Jets.Eta())<2.4&&Jets.Pt()>20&&Jets_ID) == 4 || Sum$(abs(Jets.Eta())<2.4&&Jets.Pt()>20&&Jets_ID) == 5");
  TCut Bcheck("CSV2 > .898 && CSV3 > .679 && CSV4 > .245");
 // TCut Bcheck("CSV2 > .898 && CSV4 < .244");
  //TCut Bcheck("CSV4 > .898");
  TCut bad_lumi_filter("!(RunNum==275936 && LumiBlockNum==430)");

  TCut skim("MHT > 250 && HT > 250 && DeltaPhi1>0.5&&DeltaPhi2>0.5&&DeltaPhi3>0.3&&DeltaPhi4>0.3 && @Muons.size()+@Electrons.size()==0 && isoMuonTracks+isoElectronTracks+isoPionTracks==0");

  TCut trk("@Muons.size()+@Electrons.size()==0 && isoMuonTracks+isoElectronTracks+isoPionTracks==0");
 
  TCut NJet("NJets == 4 || NJets == 5");
  TCut DeltaM("DeltaMjj < 40");
  TCut AverageM("AverageMjj > 100 && AverageMjj < 140");
  TCut Rmax("Rmax < 2.5");

  TCut Smu("@Muons.size()==1&&@Electrons.size()==0&&selectedIDIsoMuons_MTW[0]<100");
  TCut Sel("@Electrons.size()==1&&@Muons.size()==0&&selectedIDIsoElectrons_MTW[0]<100");

  TCut filters("HBHENoiseFilter && HBHEIsoNoiseFilter && eeBadScFilter && EcalDeadCellTriggerPrimitiveFilter && NVtx > 0");
  TCut v2_filters("globalTightHalo2016Filter==1&&BadChargedCandidateFilter&&BadPFMuonFilter");
  TCut mu_jet_veto = "Sum$(Jets.Pt()>200&&Jets_muonEnergyFraction>0.5&&(TVector2::Phi_mpi_pi(Jets.Phi()-METPhi)>(TMath::Pi()-0.4)))==0";
  TCut qcd_cleaning = "MET/CaloMET<5"; 
//Setup Histograms for kinematic distributions
  Double_t mht_bins[8] = {300, 350, 400, 500, 600, 750, 1000, 1500};
  TH1D* hMHT = new TH1D("hMHT", ";E_{T}^{miss} [GeV]", 25, 250, 700);
  TH1D* hHT = new TH1D("hHT", ";H_{T} [GeV]", 28, 0, 1400);
  TH1D* hNJets = new TH1D("hNJets", ";N_{jet} (p_{T} > 30 GeV)", 10, 3.5, 13.5);
  TH1D* hBTags = new TH1D("hBTags", ";N_{b-jet} (p_{T} > 30 GeV)", 5, -0.5, 4.5);
  TH1D* res = new TH1D("res", ";H_{T}^{miss}/#sqrt{H_{T}}", 18, 0,36);
  TH1D* AK8 = new TH1D("AK8", ";M_(Jet)", 30,0,300);
  TH1D* Rma = new TH1D("Rma", ";R_{Max}", 20,0,4);
  TH1D* Mav = new TH1D("Mav", ";<M_{jj}> [GeV]", 12, 0, 240);
  TH1D* Dmass = new TH1D("Dmass", ";#Deltam_{jj} [GeV]", 16, 0, 160);
  TH1D* Dphi2 = new TH1D("Dphi2", ";#Delta#phi_{2}", 32, 0, 3.2);
  TH1D* Dphi3 = new TH1D("Dphi3", ";#Delta#phi_{3}", 32, 0, 3.2);
  TH1D* Dphi4 = new TH1D("Dphi4", ";#Delta#phi_{4}", 32, 0, 3.2);
  TH1D* Jet1PT = new TH1D("Jet1PT", ";Jet 1 p_{T} [GeV]", 30, 0, 600);
  TH1D* Jet2PT = new TH1D("Jet2PT", ";Jet 2 p_{T} [GeV]", 35, 0, 350);
  TH1D* Jet3PT = new TH1D("Jet3PT", ";Jet 3 p_{T} [GeV]", 25, 0, 250);
  TH1D* BL = new TH1D("BL", ";N_{b-jet}^{L}", 7, -.5, 6.5);
  TH1D* BM = new TH1D("BM", ";N_{b-jet}^{M}", 7, -.5, 6.5);
  TH1D* BT = new TH1D("BT", ";N_{b-jet}^{T}", 7, -.5, 6.5);
//Make Plots

  Plot("DeltaMjj","MET > 250"+NJet+Rmax, bad_lumi_filter, "Weight*30000", Dmass,
           "tree_signal", "tree_signalUnblind",
           "plotLog",
           "ttbar-DeltaMjj-Nul");

  Plot("AverageMjj","MET > 250"+NJet+Rmax, bad_lumi_filter, "Weight*30000", Mav,
           "tree_signal", "tree_signalUnblind",
           "plotLog",
           "ttbar-AverageMjj-Nul");

/* 
  Plot("Rmax","MET > 250"+Rmax+AverageM+DeltaM+NJet, bad_lumi_filter, "Weight*30000", Rma,
           "tree_signal", "tree_signalUnblind",
           "plotData:plotLog",
           "New-DeltaPhi3-4b");

  Plot("Rmax","MET > 250"+Rmax+DeltaM+AverageM+NJet, bad_lumi_filter, "Weight*30000", Rma,
           "tree_signal", "tree_signalUnblind",
           "plotData:plotLog",
           "New-DeltaPhi4-4b");

  Plot("DeltaMjj",baseline+dPhi+ZL+NJ20+Rmax+AverageM+Bcheck, bad_lumi_filter, "Weight*30000", Dmass,
           "tree_signal", "tree_signalUnblind",
           "plotData:plotLog",
           "New-DeltaM-4b");

  Plot("AverageMjj",baseline+ZL+dPhi+NJ20+Rmax+DeltaM+Bcheck, bad_lumi_filter, "Weight*30000", Mav,
           "tree_signal", "tree_signalUnblind",
           "plotData:plotLog",
           "New-Maverage-4b");

  Plot("Rmax",baseline+ZL+dPhi+NJ20+DeltaM+AverageM+Bcheck, bad_lumi_filter, "Weight*30000", Rma,
           "tree_signal", "tree_signalUnblind",
           "plotData:plotLog",
           "New-Rmax-4b");

  Plot("HT","MHT>250"+ZL+dPhi+NJ20+Rmax+AverageM+DeltaM+Bcheck, bad_lumi_filter, "Weight*30000", hHT,
           "tree_signal", "tree_signalUnblind",
           "plotData:plotLog",
           "New-HT-4b");

  Plot("NJets20",baseline+dPhi+ZL+Rmax+AverageM+DeltaM+Bcheck, bad_lumi_filter, "Weight*30000", hNJets,
           "tree_signal", "tree_signalUnblind",
           "plotData:plotLog",
           "New-NJets-4b");
 */
  return; }
  





