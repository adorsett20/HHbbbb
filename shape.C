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
#include "TH1.h"
#include "TH2.h"
#include "TFileCollection.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TText.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "jack_style.h"
#include "THStack.h" 

void shape() {
  TH1::SetDefaultSumw2(); //Activates storage of the sum of squares of errors
  
  //Put cuts here
  TCut filters("HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && NVtx > 0");
//  TCut ZL("MHT>500&&HT>500&&NJets>=3 && isoMuonTracks==0 && isoElectronTracks==0 && isoPionTracks==0 &&@Electrons.size()==0 && @Muons.size()==0 && DeltaPhi1>1.5 && DeltaPhi2>1.5 && DeltaPhi3>1.8 && DeltaPhi4>1.8");
  TCut ZL("MHT>200&&HT>200 && isoMuonTracks==0 && isoElectronTracks==0 && isoPionTracks==0 &&@Electrons.size()==0 && @Muons.size()==0 && DeltaPhi1>.5 && DeltaPhi2>.5 && DeltaPhi3>.3 && DeltaPhi4>.3");
  ZL+=filters;
  TCut box1("MHT>200&&MHT<500&&HT>500&&HT<800"), box2("MHT>200&&MHT<500&&HT>800&&HT<1200"), box3("MHT>200&&MHT<500&&HT>1200"),
    box4("MHT>500&&MHT<750&&HT>500&&HT<1200"), box5("MHT>500&&MHT<750&&HT>1200"), box6("MHT>750&&HT>800");
  TCut weight = "Weight";
//Setup TChains for background processes
  TChain* ttbar = new TChain("tree");
  ttbar->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_TTJets_HT-2500toInf.root");
  ttbar->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_TTJets_HT-1200to2500.root");
  ttbar->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_TTJets_HT-600to800.root");
  ttbar->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_TTJets_HT-800to1200.root");
  TChain * ttbar_lep = new TChain("tree");
  ttbar_lep->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_TTJets_SingleLeptFromT.root");
  ttbar_lep->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_TTJets_SingleLeptFromTbar.root");
  ttbar_lep->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_TTJets_DiLept.root");
  TChain * ttbar_had = new TChain("tree");
  ttbar_had->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_TTJets.root");
  TChain * qcd = new TChain("tree");
  qcd->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_QCD_HT-700to1000.root");
  qcd->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_QCD_HT-200to300.root");
  qcd->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_QCD_HT-300to500.root");
  qcd->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_QCD_HT-500to700.root");
  qcd->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_QCD_HT-2000toInf.root");
  qcd->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_QCD_HT-1500to2000.root");
  qcd->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_QCD_HT-1000to1500.root");
  TChain * znn = new TChain("tree");
  znn->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_ZJetsToNuNu_HT-400to600.root");
  znn->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_ZJetsToNuNu_HT-600to800.root");
  znn->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_ZJetsToNuNu_HT-200to400.root");
  znn->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_ZJetsToNuNu_HT-100to200.root");
  znn->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_ZJetsToNuNu_HT-800to1200.root");
  znn->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_ZJetsToNuNu_HT-1200to2500.root");
  znn->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_ZJetsToNuNu_HT-2500toInf.root");
  TChain * wjets = new TChain("tree");
  wjets->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_WJetsToLNu_HT-1200to2500.root");
  wjets->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_WJetsToLNu_HT-2500toInf.root");
  wjets->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_WJetsToLNu_HT-400to600.root");
  wjets->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_WJetsToLNu_HT-800to1200.root");
  wjets->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_WJetsToLNu_HT-600to800.root");
  wjets->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_WJetsToLNu_HT-200to400.root");
  wjets->Add("root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV9/tree_signal/tree_WJetsToLNu_HT-100to200.root");
//Setup TChains for signal processes
  TChain* sig1 = new TChain("TreeMaker2/PreSelection");
  sig1->Add("Total.root");
//Initialize Histograms
  TH1D* hMHT = new TH1D("hMHT", ";H_{T}^{miss} [GeV]; Fraction of events", 16, 200, 1000);
  TH1D* hBTags = new TH1D("hBTags",";NJets; Fraction of events",8, -.5,7.5);
  TH1D* ttbarh = (TH1D*)hBTags->Clone("ttbar");
  TH1D* qcdh = (TH1D*)hBTags->Clone("qcd");
  TH1D* znnh = (TH1D*)hBTags->Clone("znn");
  TH1D* wjetsh = (TH1D*)hBTags->Clone("wjets");
  TH1D* sig1h = (TH1D*)hBTags->Clone("sig1");
//Project tree variables into histograms
  ttbar->Project("ttbar", "Sum$(abs(Jets.Eta())<2.4&&Jets.Pt()>20&&Jets_ID)", ZL*weight);
  ttbar_lep->Project("+ttbar", "Sum$(abs(Jets.Eta())<2.4&&Jets.Pt()>20&&Jets_ID)", (ZL+"madHT<600")*weight);
  qcd->Project("qcd", "Sum$(abs(Jets.Eta())<2.4&&Jets.Pt()>20&&Jets_ID)", ZL*weight);
  ttbar_had->Project("+qcd","Sum$(abs(Jets.Eta())<2.4&&Jets.Pt()>20&&Jets_ID)",(ZL+"madHT<600"+"@GenEls.size()+@GenMus.size()+@GenTaus.size()==0")*weight);
  znn->Project("znn","Sum$(abs(Jets.Eta())<2.4&&Jets.Pt()>20&&Jets_ID)",ZL*weight);
  wjets->Project("wjets","Sum$(abs(Jets.Eta())<2.4&&Jets.Pt()>20&&Jets_ID)",ZL*weight);
  sig1->Project("sig1","Sum$(abs(Jets.Eta())<2.4&&Jets.Pt()>20&&Jets_ID)",ZL*weight);
//Fix overflow by adding to the last bin in histogram
//Also Normalizes Histograms
  int nbinsx=hMHT->GetNbinsX();
  Double_t e_overflow(0.);
  ttbarh->SetBinContent(nbinsx, ttbarh->IntegralAndError(nbinsx, nbinsx+1, e_overflow));
  ttbarh->SetBinError(nbinsx, e_overflow);
  ttbarh->Scale(1/ttbarh->Integral(1, nbinsx));
  qcdh->SetBinContent(nbinsx, qcdh->IntegralAndError(nbinsx, nbinsx+1, e_overflow));
  qcdh->SetBinError(nbinsx, e_overflow);
  qcdh->Scale(1/qcdh->Integral(1, nbinsx));
  znnh->SetBinContent(nbinsx, znnh->IntegralAndError(nbinsx, nbinsx+1, e_overflow));
  znnh->SetBinError(nbinsx, e_overflow);
  znnh->Scale(1/znnh->Integral(1, nbinsx));
  wjetsh->SetBinContent(nbinsx, wjetsh->IntegralAndError(nbinsx, nbinsx+1, e_overflow));
  wjetsh->SetBinError(nbinsx, e_overflow);
  wjetsh->Scale(1/wjetsh->Integral(1, nbinsx));
  sig1h->SetBinContent(nbinsx, sig1h->IntegralAndError(nbinsx, nbinsx+1, e_overflow));
  sig1h->SetBinError(nbinsx, e_overflow);
  sig1h->Scale(1/sig1h->Integral(1, nbinsx));

//Setup Canvas and Plot
  int W = 800;
  int H = 800;

  // 
  int H_ref = 800;
  int W_ref = 800; 
// references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.1*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  TCanvas* canv = new TCanvas("canvName","canvName",50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);

  ttbarh->SetStats(0);
  qcdh->SetStats(0);
  znnh->SetStats(0);
  wjetsh->SetStats(0);
  sig1h->SetStats(0);
 
  ttbarh->SetLineWidth(3);
  qcdh->SetLineWidth(3);
  znnh->SetLineWidth(3);
  wjetsh->SetLineWidth(3);
  sig1h->SetLineWidth(4);
  ttbarh->SetLineColor(2);
  qcdh->SetLineColor(3); 
  znnh->SetLineColor(4);
  wjetsh->SetLineColor(5);

  sig1h->SetLineStyle(7);
  sig1h->SetLineColor(kBlack);

//  ttbarh->SetFillColor(2);
//  qcdh->SetFillColor(3);
 // znnh->SetFillColor(4);
//  wjetsh->SetFillColor(5);
//  sig1h->SetFillColor(1);

  ttbarh->GetYaxis()->SetLabelSize(0.035*1.);
  ttbarh->GetYaxis()->SetTitleSize(0.04*1.);
  ttbarh->GetYaxis()->SetTitleOffset(1.35);
  ttbarh->GetYaxis()->SetTitleFont(42);

  ttbarh->GetXaxis()->SetLabelSize(0.035*1.);
  ttbarh->GetXaxis()->SetTitleSize(0.04*1.);
  ttbarh->GetXaxis()->SetTitleOffset(1.15);
  ttbarh->GetXaxis()->SetTitleFont(42);

  ttbarh->SetMaximum(wjetsh->GetMaximum()*1.5);

  znnh->GetYaxis()->SetRangeUser(0,1.0);
 // canv->SetLogy();
  THStack* stack = new THStack("stack","");
  stack->Add(ttbarh);
  stack->Add(qcdh);
  stack->Add(znnh);
  stack->Add(wjetsh);
  stack->Add(sig1h);

  znnh->Draw("hist");
  qcdh->Draw("hist,same");
  ttbarh->Draw("hist,same");
  wjetsh->Draw("hist,same");
  sig1h->Draw("hist,same");
//  stack->Draw("hist");

//  stack->GetXaxis()->SetTitle("H_{T}^{miss} [GeV]");
 // stack->GetYaxis()->SetTitle("Fraction of Events");

  TLegend * leg1;
  TLegend * leg2;
  TLegend * leg3;
  TLegend * leg4;
  TLegend * legsig;

  leg1 = new TLegend(0.17, 0.856, 0.58, 0.92);
  leg2 = new TLegend(0.33, 0.856, 0.73, 0.92);
  leg3 = new TLegend(0.52, 0.856, 0.94, 0.92);
  leg4 = new TLegend(0.735, 0.856, 1.12, 0.92);
  set_style(leg1,0.03);
  set_style(leg2,0.03);
  set_style(leg3,0.03);
  set_style(leg4,0.03);
  leg1->AddEntry(ttbarh, "t#bar{t}", "l");
  leg2->AddEntry(qcdh, "QCD", "l");
  leg3->AddEntry(znnh, "Z+jets", "l");
  leg4->AddEntry(wjetsh, "W+jets", "l");

  legsig = new TLegend(0.17, 0.77, 0.65, 0.87);
  set_style(legsig,0.03);
  legsig->SetMargin(0.2);
  legsig->AddEntry(sig1h, "pp #rightarrow #tilde{#chi}_{1}^{0}#tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow H #tilde{G} (m_{#tilde{G}} = 1 GeV, m_{#tilde{#chi}_{1}^{0}} = 400 GeV)", "l");

  leg1->Draw();
  leg2->Draw();
  leg3->Draw();
  leg4->Draw();
  legsig->Draw();

  TLatex * latex = new TLatex();
  latex->SetNDC();

  float t=canv->GetTopMargin();
  float r=canv->GetRightMargin();
  TString cmsText     = "CMS";
  float cmsTextFont   = 61;
  float cmsTextSize      = 0.55;
  latex->SetTextFont(cmsTextFont);
  latex->SetTextSize(cmsTextSize*t);
  latex->DrawLatex(0.12, 0.95, cmsText);

  TString extraText   = "            Supplementary (Simulation)";
  latex->SetTextFont(52);
  float extraOverCmsTextSize  = 0.76;
  float extraTextSize = extraOverCmsTextSize*cmsTextSize;
  latex->SetTextSize(extraTextSize*t);
  float relExtraDY = 1.2;
  latex->DrawLatex(0.12, 0.95, extraText);

  latex->SetTextFont(42);
  latex->SetTextAlign(31);
  float lumiTextSize=0.6;
  float lumiTextOffset=0.2;

  latex->SetTextSize(lumiTextSize*t);
  latex->DrawLatex(1-r,1-t+lumiTextOffset*canv->GetTopMargin(), "(13 TeV)");

  canv->Print("NJets-sig20.pdf");
  canv->Print("NJets-sig20.png");

/*  delete canv;
  delete ttbarh;
  delete qcdh ;
  delete znnh ;
  delete wjetsh;
  delete sig1h;
  delete ttbar;
  delete qcd ;
  delete znn ;
  delete wjets;
  delete sig1;
*/
  return; }
