//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 24 21:50:43 2023 by ROOT version 6.24/08
// from TTree outTree/outTree
// found on file: user.rlindley.34485398._000003.tree.root
//////////////////////////////////////////////////////////
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "vector"
#include <TMath.h>
#include "algorithm"
#include <string>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

int analysis(){

  TH1F *histo_1around;
  TFile *f_1around = new TFile("road_barcodeeff_nopileup_1around.root");
  f_1around->GetObject("eff_d0",histo_1around);

  TH1F *histo_2around_avg;
  TFile *f_2around_avg = new TFile("road_barcodeeff_nopileup_2around_avg.root");
  f_2around_avg->GetObject("eff_d0",histo_2around_avg);
  
  TH1F *histo_2around;
  TFile *f_2around = new TFile("road_barcodeeff_nopileup_2around.root");
  f_2around->GetObject("eff_d0",histo_2around);

  TCanvas *c = new TCanvas();
  histo_1around->SetMinimum(0.87);
  histo_1around->SetMaximum(1.0);
  histo_1around->SetLineColor(kBlack);
  histo_1around->SetLineStyle(2);
  histo_2around->SetLineColor(kBlue);
  histo_2around->SetLineStyle(2);
  histo_2around_avg->SetLineColor(kRed);
  histo_2around_avg->SetLineStyle(2);
  histo_1around->Draw("Hist");
  histo_2around->Draw("Same");
  histo_2around_avg->Draw("Same");

  auto leg = new TLegend(.55, .6, .85, .89, "Road Efficiency");
  leg->AddEntry(histo_1around,"8 bins");
  leg->AddEntry(histo_2around,"24 bins");
  leg->AddEntry(histo_2around_avg,"24 bins + avg");
  leg->Draw("HistSame");
  c->Print("road_barcodeeff_comparison_d0.pdf"); 

  return 0;
}
