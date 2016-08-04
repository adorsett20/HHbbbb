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
	  TString options="plotSig:plotLog:plotData", TString title="default") {
//Boolians for plot drawing options
  bool plotLog = options.Contains("plotLog") && (!options.Contains("!plotLog"));
  bool plotData = options.Contains("plotData") && (!options.Contains("!plotData"));
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
  TChain* ttbar_had = new TChain("tree");
  ttbar_had->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTJets.root");
  TChain* qcd = new TChain("tree");
  qcd->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_QCD_HT-700to1000.root");
  qcd->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_QCD_HT-200to300.root");
  qcd->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_QCD_HT-300to500.root");
  qcd->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_QCD_HT-500to700.root");
  qcd->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_QCD_HT-2000toInf.root");
  qcd->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_QCD_HT-1500to2000.root");
  qcd->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_QCD_HT-1000to1500.root");
  TChain* znn = new TChain("tree");
  znn->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ZJetsToNuNu_HT-400to600.root");
  znn->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ZJetsToNuNu_HT-600to800.root");
  znn->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ZJetsToNuNu_HT-200to400.root");
  znn->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ZJetsToNuNu_HT-100to200.root");
  znn->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ZJetsToNuNu_HT-800to1200.root");
  znn->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ZJetsToNuNu_HT-1200to2500.root");
  znn->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ZJetsToNuNu_HT-2500toInf.root");
  TChain* wjets = new TChain("tree");
  wjets->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WJetsToLNu_HT-1200to2500.root");
  wjets->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WJetsToLNu_HT-2500toInf.root");
  wjets->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WJetsToLNu_HT-400to600.root");
  wjets->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WJetsToLNu_HT-800to1200.root");
  wjets->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WJetsToLNu_HT-600to800.root");
  wjets->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WJetsToLNu_HT-200to400.root");
  wjets->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WJetsToLNu_HT-100to200.root");
  TChain* other = new TChain("tree");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WWTo2L2Nu.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTWJetsToLNu.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTZToLLNuNu.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTTT.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTZToQQ.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_TTWJetsToQQ.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ZZZ.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WZZ.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WWZ.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ST_s-channel.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_DYJetsToLL_M-50_HT-200to400.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WZTo1L3Nu.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_DYJetsToLL_M-50_HT-600toInf.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_DYJetsToLL_M-50_HT-400to600.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_DYJetsToLL_M-50_HT-100to200.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WWTo1L1Nu2Q.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ZZTo2L2Q.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_WZTo1L1Nu2Q.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ZZTo2Q2Nu.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ST_tW_top.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ST_tW_antitop.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ST_t-channel_antitop.root");
  other->Add("root://cmseos.fnal.gov//store/user/adorsett/HHbbbb/BackTrees/tree_ST_t-channel_top.root");
  TChain* data = new TChain("TreeMaker2/PreSelection");
//  TChain* sig1 = new TChain("TreeMaker2/PreSelection");
  data->Add("Signal.root");
//Make output file
  TFile* outfile = new TFile(Form("hists/%s.root", title.Data()), "recreate"); 
//Make Histograms from input temp
  TH1D * httbar = (TH1D*)temp->Clone("ttbar");
  TH1D * hqcd = (TH1D*)temp->Clone("qcd");
  TH1D * hznn = (TH1D*)temp->Clone("znn");
  TH1D * hwjets = (TH1D*)temp->Clone("wjets");
  TH1D * hother = (TH1D*)temp->Clone("other");
  TH1D * hdata = (TH1D*)temp->Clone("data");
  TH1D * hmc = (TH1D*)temp->Clone("mc"); //Used to keep track of stack values
//Additional Cuts
  TCut selection(OtherCuts);
//Fill histograms
  ttbar->Project("ttbar", var, selection*weight);
  ttbar_lep->Project("+ttbar", var, (selection+"madHT<600")*weight);
  qcd->Project("qcd", var, (selection)*weight);
  ttbar_had->Project("+qcd", var, (selection+"madHT<600&&@GenEls.size()+@GenMus.size()+@GenTaus.size()==0")*weight);
  znn->Project("znn", var, selection*weight);
  wjets->Project("wjets", var, selection*weight);
  other->Project("other", var, selection*weight);
  data->Project("data",var,selection*"30000.*0.00887325/28405.");
//  hdata->Scale(30000./7628.574);
//Fix overflow by adding to last bin in hist
  bool addOverflow(true);
  Double_t e_overflow(0.), i_overflow(0.);
  int nbinsx=temp->GetNbinsX();
  if (addOverflow) {
    i_overflow=httbar->IntegralAndError(nbinsx,nbinsx+1,e_overflow);
    httbar->SetBinContent(nbinsx, i_overflow);
    httbar->SetBinError(nbinsx, e_overflow);
    i_overflow=hqcd->IntegralAndError(nbinsx,nbinsx+1,e_overflow);
    hqcd->SetBinContent(nbinsx, i_overflow);
    hqcd->SetBinError(nbinsx, e_overflow);
    i_overflow=hznn->IntegralAndError(nbinsx,nbinsx+1,e_overflow);
    hznn->SetBinContent(nbinsx, i_overflow);
    hznn->SetBinError(nbinsx, e_overflow);
    i_overflow=hwjets->IntegralAndError(nbinsx,nbinsx+1,e_overflow);
    hwjets->SetBinContent(nbinsx, i_overflow);
    hwjets->SetBinError(nbinsx, e_overflow);
    i_overflow=hother->IntegralAndError(nbinsx,nbinsx+1,e_overflow);
    hother->SetBinContent(nbinsx, i_overflow);
    hother->SetBinError(nbinsx, e_overflow);
    i_overflow=hdata->IntegralAndError(nbinsx,nbinsx+1,e_overflow);
    hdata->SetBinContent(nbinsx, i_overflow);
    hdata->SetBinError(nbinsx, e_overflow);
  }
//Output number of events present in data
  printf("Selected %3.0f events in data.\n", hdata->Integral());
//Use set_style function on histograms
  set_style_lite(httbar, "ttbar");
  httbar->SetFillColor(kAzure+7);
  set_style(hqcd, "qcd");
  set_style_lite(hznn, "znn");
  set_style_lite(hwjets, "wjets");
  set_style_lite(hother, "single_top");
  set_style(hdata, "data");
  hdata->SetFillColor(0);
 // hdata->SetMarkerColor(kMagenta+1);
  hdata->SetLineColor(kMagenta+1);
//Create MC histogram with backgrounds
  hmc->Add(httbar);
  hmc->Add(hqcd);
  hmc->Add(hznn);
  hmc->Add(hwjets);
  hmc->Add(hother);
//Print out total events in MC background
  printf("Selected %3.2f events in MC.\n", hmc->Integral());
//Write histograms to output file
  outfile->cd();
  hdata->Write();
  httbar->Write();
  hqcd->Write();
  hznn->Write();
  hwjets->Write();
  hother->Write();
  hmc->Write();
//Make MC histogram invsible
  hmc->SetMarkerSize(0);
  hmc->SetLineWidth(0);
//Use MC hist to label Axes
  TString ytitle = Form("Events");
  if (entries) ytitle = Form("Entries");
  hmc->GetXaxis()->SetTitle(httbar->GetXaxis()->GetTitle());
  hmc->GetYaxis()->SetTitle(ytitle);
//Make HStack for backgrounds
  THStack * hs = new THStack("hs", "");
//Stack histograms in order from highest to lowest contribution
  TH1D* hists[5] = {httbar, hqcd, hznn, hwjets, hother};
  TString hlabels[5] = {"t#bar{t}", "QCD", "Z+jets", "W+jets", "Other"};
  vector<int> hist_ids;
  bool hist_done[5] = {false, false, false, false, false};
  for (int ih=0; ih<5; ih++) {
    int jmin=-1;
    for (int jh=0; jh<5; jh++) {
      if (hist_done[jh]) continue;
      if (jmin<0) jmin = jh;
      if (hists[jh]->Integral()<hists[jmin]->Integral()) jmin=jh;
    }
    hist_done[jmin]=true;
    hist_ids.push_back(jmin);
    hs->Add(hists[jmin]);
  }
//Setup additional histograms
  TH1D * staterr = (TH1D *) hmc->Clone("staterr");
  staterr->SetFillColor(1);
  staterr->SetLineColor(1);
  staterr->SetMarkerSize(0);
  staterr->SetLineWidth(0);
  staterr->SetFillStyle(3354);
  TH1D * ratio = (TH1D *) hdata->Clone("ratio");
  ratio->Divide(hdata, hmc);
  ratio->SetStats(0);
  ratio->SetTitle(hmc->GetTitle());
  ratio->GetYaxis()->SetTitle("S/B");
  ratio->SetMaximum(.7);
  ratio->SetMinimum(0.);
  ratio->GetXaxis()->SetLabelSize(0.15);
  ratio->GetXaxis()->SetLabelOffset(0.03);
  ratio->GetXaxis()->SetTitleSize(0.14);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetLabelSize(0.10);
  ratio->GetYaxis()->SetTitleSize(0.12);
  ratio->GetYaxis()->SetTitleOffset(0.6);
  ratio->GetYaxis()->SetNdivisions(505);
//Setup Legends
  TLegend * leg1 = new TLegend(0.53, 0.6, 0.77, 0.92);
  set_style(leg1,0.035);
  if (plotData) leg1->AddEntry(hdata, "Signal", "pel");
  leg1->AddEntry(hists[hist_ids[4]], hlabels[hist_ids[4]].Data(), "f");
  leg1->AddEntry(hists[hist_ids[3]], hlabels[hist_ids[3]].Data(), "f");
  TLegend * leg2 = new TLegend(0.72, 0.6, 0.94, 0.92);
  set_style(leg2,0.035);
  leg2->AddEntry(hists[hist_ids[2]], hlabels[hist_ids[2]].Data(), "f");
  leg2->AddEntry(hists[hist_ids[1]], hlabels[hist_ids[1]].Data(), "f");
  leg2->AddEntry(hists[hist_ids[0]], hlabels[hist_ids[0]].Data(), "f");
  double ymax = hs->GetMaximum();
  if (hdata->GetMaximum()>ymax) ymax=hdata->GetMaximum();
  if(plotLog) {
    hmc->SetMaximum(1000*ymax);
    hmc->SetMinimum(0.08);
    if (var.Contains("num_")) hmc->SetMaximum(1000*ymax);
    if (var.Contains("eta")) hmc->SetMaximum(1000*ymax);
    if (var.Contains("dphi")) hmc->SetMaximum(1000*ymax);
  }
  else {
    hmc->SetMinimum(0);
    hmc->SetMaximum(1.75*ymax);
    if (title.Contains("cdtt")) hmc->SetMaximum(1.75*ymax);
    if (title.Contains("/")) hmc->SetMaximum(1.75*ymax);
    if (var.Contains("eta")) hmc->SetMaximum(3*ymax);
  }
//Canvas and Pad Setup
  TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
  TPad * pad1 = new TPad("pad1", "top pad" , 0.0, 0.31, 1.0, 1.0);
  TPad * pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.3);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad1->cd();
  pad1->SetLogy(plotLog);
  hmc->Draw("hist");
  hmc->GetXaxis()->SetLabelSize(0.03);
  hmc->GetXaxis()->SetTitleSize(0.025);
  hmc->GetYaxis()->SetLabelSize(0.04);
  hmc->SetMarkerSize(0);
  hmc->SetLineWidth(0);
//Draw Histograms
  hs->Draw("hist,same");
  staterr->Draw("e2 same");
  if (plotData) hdata->Draw("hist,same");
  hmc->GetXaxis()->SetLabelSize(0);
//Draw Legends
  leg1->Draw();
  leg2->Draw();
  TLatex * latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextFont(62);
  latex->SetTextSize(0.052);
  latex->DrawLatex(0.19, 0.89, "CMS Simulation");
  latex->SetTextSize(0.04);
  latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 13 TeV, L = 30 fb^{-1}");
//Setup Second Pad
  pad2->cd();
  pad2->SetGridy(0);
  ratio->Draw("e1");
  TLine* ratiounity = new TLine(hmc->GetBinLowEdge(1),1,hmc->GetBinLowEdge(hmc->GetNbinsX()+1),1);
  ratiounity->SetLineStyle(1);
  ratiounity->Draw("same");
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetLabelSize(0.12);
  hmc->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetLabelSize(0.12);
  TPaveText * pave = new TPaveText(0.18, 0.86, 0.4, 0.96, "brNDC");
  pave->SetLineColor(0);
  pave->SetFillColor(kWhite);
  pave->SetShadowColor(0);
  pave->SetBorderSize(1);

  pad1->cd();
  gPad->RedrawAxis();
  gPad->Modified();
  gPad->Update();
  pad2->cd();
  gPad->RedrawAxis();
  gPad->Modified();
  gPad->Update();
  c1->cd();

  c1->cd();
  if (title.EqualTo("default")) title=plotdir+var;
  gPad->Print(plotdir+title+".pdf");   
//Clean up histograms
  delete ttbar;
  delete qcd;
  delete znn;
  delete wjets;
  delete other;
  delete data;
  delete staterr;
  delete ratio;
  delete leg1;
  delete leg2;
  delete latex;
  delete pave;
  delete hs;
  delete pad1;
  delete pad2;
  delete c1;

  outfile->Close();
  return; }

void HHstackPlus() { 
  TH1::SetDefaultSumw2(1);
  if (gSystem->AccessPathName(plotdir))
    gSystem->mkdir(plotdir);
//Style Setup
  cout << "Setting tdr style."  << endl;
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  setTDRStyle(tdrStyle);
  tdrStyle->cd();
//Put in Cuts
  TCut baseline("MHT>250&&HT>300&&BTags>1");
  TCut ZL("@Muons.size()+@Electrons.size()==0 && isoMuonTracks+isoElectronTracks+isoPionTracks==0");
  TCut dPhi("DeltaPhi1>0.5&&DeltaPhi2>0.5&&DeltaPhi3>0.3&&DeltaPhi4>0.3");

  TCut antiPhi("DeltaPhi1>2&&DeltaPhi2>2&&DeltaPhi3>1.8&&DeltaPhi4>1.8");

  TCut JetFix("Sum$(abs(Jets.Eta())<2.4&&Jets.Pt()>20&&Jets_ID) == 4 || Sum$(abs(Jets.Eta())<2.4&&Jets.Pt()>20&&Jets_ID) == 5");
  TCut Bcheck("Sum$(abs(Jets.Eta())<2.4&&Jets.Pt()>20&&Jets_ID&&Jets_bDiscriminatorCSV>0.898)>=3");

  TCut bad_lumi_filter("!(RunNum==275936 && LumiBlockNum==430)");

  TCut Smu("@Muons.size()==1&&@Electrons.size()==0&&selectedIDIsoMuons_MTW[0]<100");
  TCut Sel("@Electrons.size()==1&&@Muons.size()==0&&selectedIDIsoElectrons_MTW[0]<100");

  TCut filters("HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && NVtx > 0");
  TCut v2_filters("globalTightHalo2016Filter==1&&BadChargedCandidateFilter&&BadPFMuonFilter");
  TCut mu_jet_veto = "Sum$(Jets.Pt()>200&&Jets_muonEnergyFraction>0.5&&(TVector2::Phi_mpi_pi(Jets.Phi()-METPhi)>(TMath::Pi()-0.4)))==0";
  TCut qcd_cleaning = "MET/CaloMET<5"; 
//Setup Histograms for kinematic distributions
  Double_t mht_bins[8] = {300, 350, 400, 500, 600, 750, 1000, 1500};
  TH1D* hMHT = new TH1D("hMHT", ";H_{T}^{miss} [GeV]", 20, 250, 1250);
  TH1D* hHT = new TH1D("hHT", ";H_{T} [GeV]", 23, 300, 2600);
  TH1D* hNJets = new TH1D("hNJets", ";N_{jet} (p_{T} > 50 GeV)", 9, -0.5, 8.5);
  TH1D* hBTags = new TH1D("hBTags", ";N_{b-jet} (p_{T} > 30 GeV)", 5, -0.5, 4.5);
  TH1D* res = new TH1D("res", ";H_{T}^{miss}/#sqrt{H_{T}}", 18, 0,36);
  TH1D* AK8 = new TH1D("AK8", ";M_(Jet)", 30,0,300);
//Make Plots

  Plot("MHT/sqrt(HT)",ZL+baseline+filters+dPhi+Bcheck+JetFix, bad_lumi_filter, "Weight*30000", res,
           "tree_signal", "tree_signalUnblind",
           "plotData:plotLog",
           "New-MHTrootHT");
/*
  Plot("HT",ZL+baseline+filters+antiPhi+"MHT>500", htmht_trigger+v2_filters+bad_lumi_filter, "Weight*7628.574", hHT,
           "tree_signal","tree_signalUnblind",
           "plotData",
           "Signal-AntiPhi-HT");
  Plot("NJets",ZL+baseline+filters+antiPhi+"MHT>500", htmht_trigger+v2_filters+bad_lumi_filter, "Weight*7628.574", hNJets,
           "tree_signal","tree_signalUnblind",
           "plotData",
           "Signal-AntiPhi-NJets");
  Plot("BTags",ZL+baseline+filters+antiPhi+"MHT>500", htmht_trigger+v2_filters+bad_lumi_filter, "Weight*7628.574", hBTags,
           "tree_signal","tree_signalUnblind",
           "plotData",
           "Signal-AntiPhi-BTags");
*/
/*
  Plot("HT", ZL+baseline+filters+"(NJets>5&&(abs(TVector2::Phi_mpi_pi(Jets[5].Phi()-MHTPhi))<.3||abs(TVector2::Phi_mpi_pi(Jets[6].Phi()-MHTPhi)) <.3))||(NJets==5&&abs(TVector2::Phi_mpi_pi(Jets[5].Phi()-MHTPhi))<.3)"+"BTags==0", htmht_trigger+v2_filters+bad_lumi_filter, "Weight*7628.574", hHT,
	   "tree_signal", "tree_signalUnblind",
	   "plotData",
	   "Signal-LowPhi56-HT");  
  Plot("MHT", ZL+baseline+filters+"(NJets>5&&(abs(TVector2::Phi_mpi_pi(Jets[5].Phi()-MHTPhi))<.3||abs(TVector2::Phi_mpi_pi(Jets[6].Phi()-MHTPhi)) <.3))||(NJets==5&&abs(TVector2::Phi_mpi_pi(Jets[5].Phi()-MHTPhi))<.3)"+"BTags==0", htmht_trigger+v2_filters+bad_lumi_filter, "Weight*7628.574", hMHT,
           "tree_signal", "tree_signalUnblind",
           "plotData",
           "Signal-LowPhi56-MHT");
  Plot("NJets", ZL+baseline+filters+"(NJets>5&&(abs(TVector2::Phi_mpi_pi(Jets[5].Phi()-MHTPhi))<.3||abs(TVector2::Phi_mpi_pi(Jets[6].Phi()-MHTPhi)) <.3))||(NJets==5&&abs(TVector2::Phi_mpi_pi(Jets[5].Phi()-MHTPhi))<.3)"+"BTags==0", htmht_trigger+v2_filters+bad_lumi_filter, "Weight*7628.574", hNJets,
           "tree_signal", "tree_signalUnblind",
           "plotData",
           "Signal-LowPhi56-NJets");
  Plot("BTags", ZL+baseline+filters+"(NJets>5&&(abs(TVector2::Phi_mpi_pi(Jets[5].Phi()-MHTPhi))<.3||abs(TVector2::Phi_mpi_pi(Jets[6].Phi()-MHTPhi)) <.3))||(NJets==5&&abs(TVector2::Phi_mpi_pi(Jets[5].Phi()-MHTPhi))<.3)"+"BTags==0", htmht_trigger+v2_filters+bad_lumi_filter, "Weight*7628.574", hBTags,
           "tree_signal", "tree_signalUnblind",
           "plotData",
           "Signal-LowPhi56-BTags");  
*/

  return; }
  





