//Chad Harrington 7/28/2015, root -l -b -q l1fastjet.c

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"

using namespace std;

void l1fastjet(){

  int Rlabel = 4;
  TString pf_type = "chs";
  bool nPU_derived = false;

  ifstream scale_file( Form("./plots/scalefactor/DataMcSF_L1FastJet_AK%iPF", Rlabel) + pf_type + ".txt" );

  ifstream mc_file( Form("Summer15_25nsV2_MC_L1FastJet_AK%iPF", Rlabel) + pf_type + ".txt" );
  string scale_line, mc_line;

  //read first line
  getline(scale_file, scale_line);
  getline(mc_file, mc_line);

  double scale_eta1[100], scale_eta2[100], scale_p0[100], scale_p1[100], scale_p2[100];
  double mc_eta1[100], mc_eta2[100], mc_p0[100], mc_p1[100];

  for (int i=0; getline(scale_file,scale_line); i++){

    string str;
    int delim_pos;

    while (scale_line.at(0) == ' ')  //eta values have space in front of them
      scale_line.erase(0, 1);

    //loop over strings in data line
    for (int col_num=0; (delim_pos = scale_line.find(' ')) != -1; col_num++){

      str = scale_line.substr(0, delim_pos);
      scale_line.erase(0, delim_pos + 1);

      while (scale_line.at(0) == ' ')  //get rid of white space between columns
        scale_line.erase(0, 1);

      if (col_num == 0) scale_eta1[i] = stod(str);
      else if (col_num == 1) scale_eta2[i] = stod(str);
      else if (col_num == 5) scale_p0[i] = stod(str);
      else if (col_num == 6) scale_p1[i] = stod(str);
    }
    //last column after loop
    scale_p2[i] = stod(scale_line);
  }
  scale_file.close();

  int lines;
  for (int i=0; getline(mc_file,mc_line); i++){

    lines = i;
    string str;
    int delim_pos;

    while (mc_line.at(0) == ' ')  //eta values have space in front of them
      mc_line.erase(0, 1);

    //loop over strings in mc line
    for (int col_num=0; (delim_pos = mc_line.find(' ')) != -1; col_num++){

      str = mc_line.substr(0, delim_pos);
      mc_line.erase(0, delim_pos + 1);

      while (!mc_line.empty() && mc_line.at(0) == ' ')  //get rid of white space between columns
        mc_line.erase(0, 1);

      if (col_num == 0) mc_eta1[i] = stod(str);
      else if (col_num == 1) mc_eta2[i] = stod(str);
      else if (col_num == 9) mc_p0[i] = stod(str);
      else if (col_num == 10) mc_p1[i] = stod(str);
    }
  }
  const int mc_lines = lines + 1;

  TFile* data_root = TFile::Open( Form("RunII_DataD_2509_R%i.root", Rlabel) );
  TFile* mc_root = TFile::Open( Form("RunII_MCD_2509_R%i.root", Rlabel) );

  TString hname;

  if (nPU_derived) hname = "nPU";
  else hname = "nPV";

  TH1F* h_nPU = (TH1F*) data_root->Get(hname);
  double mean_nPU = h_nPU->GetMean();

  int mean_bin = h_nPU->FindBin(mean_nPU);
  int low_bin, high_bin;
  double total_area = h_nPU->Integral();

  for (low_bin = mean_bin; (h_nPU->Integral(low_bin, mean_bin) / total_area) < 0.34; low_bin--){}
  for (high_bin = mean_bin; (h_nPU->Integral(mean_bin, high_bin) / total_area) < 0.34; high_bin++){}

  double nPU_low = h_nPU->GetBinCenter(low_bin);
  double nPU_high = h_nPU->GetBinCenter(high_bin);

  if (nPU_derived) hname = "p_rho_nPU";
  else hname = "p_rho_nPV";

  TProfile* data_rho_nPU = (TProfile*) data_root->Get(hname);
  TProfile* mc_rho_nPU = (TProfile*) mc_root->Get(hname);

  double rho_nominal = data_rho_nPU->GetBinContent( data_rho_nPU->FindBin( mean_nPU ) );
  double rho_low = data_rho_nPU->GetBinContent( data_rho_nPU->FindBin( nPU_low ) );
  double rho_high = data_rho_nPU->GetBinContent( data_rho_nPU->FindBin( nPU_high ) );

  double mc_rho_mean = mc_rho_nPU->GetBinContent( mc_rho_nPU->FindBin( mean_nPU ) );
  double mc_rho_low = mc_rho_nPU->GetBinContent( mc_rho_nPU->FindBin( nPU_low ) );
  double mc_rho_high = mc_rho_nPU->GetBinContent( mc_rho_nPU->FindBin( nPU_high ) );

  double eta[mc_lines], sf[mc_lines], sf_low[mc_lines], sf_high[mc_lines], new_p0[mc_lines], new_p1[mc_lines];

  int j = 0; //the scale index

  for (int i=0; i<mc_lines; i++){

    double scale_eta = (scale_eta1[j]+scale_eta2[j])/2;
    eta[i] = (mc_eta1[i]+mc_eta2[i])/2;

    //match eta bins
    if (j == 0){
      double scale_etaUP = (scale_eta1[j+1]+scale_eta2[j+1])/2;

      if ( abs(scale_etaUP - eta[i]) < abs(scale_eta - eta[i]) ) j++;
    }
    else if (j == 99){
      double scale_etaDOWN = (scale_eta1[j-1]+scale_eta2[j-1])/2;

      if ( abs(scale_etaDOWN - eta[i]) < abs(scale_eta - eta[i]) ) j--;
    }
    else{
      double scale_etaUP = (scale_eta1[j+1]+scale_eta2[j+1])/2;
      double scale_etaDOWN = (scale_eta1[j-1]+scale_eta2[j-1])/2;

      double diff = abs(scale_eta - eta[i]); double diffUP = abs(scale_etaUP - eta[i]); double diffDOWN = abs(scale_etaDOWN - eta[i]);

      if (diffUP < diff && diffUP < diffDOWN) j++;
      else if (diffDOWN < diffUP && diffDOWN < diff) j--;
    }

    sf_low[i] = scale_p0[j] + scale_p1[j]*rho_low + scale_p2[j]*rho_low*rho_low;
    sf_high[i] = scale_p0[j] + scale_p1[j]*rho_high + scale_p2[j]*rho_high*rho_high;

    sf[i] = scale_p0[j] + scale_p1[j]*rho_nominal + scale_p2[j]*rho_nominal*rho_nominal;
    new_p0[i] = sf[i] * mc_p0[i];
    new_p1[i] = sf[i] * mc_p1[i];

    j++;
  }

  ofstream writeFile( Form("Summer15_25nsV0_DATA_L1FastJet_AK%iPF", Rlabel) + pf_type + ".txt" );
  mc_file.clear();
  mc_file.seekg(0, mc_file.beg);

  //write first line
  getline(mc_file,mc_line);
  writeFile << mc_line << endl;

  for (int i=0; getline(mc_file,mc_line); i++){

    string str;
    int delim_pos;

    while (mc_line.at(0) == ' ')  //eta values have space in front of them
      mc_line.erase(0, 1);

    //loop over strings in mc line
    for (int col_num=0; (delim_pos = mc_line.find(' ')) != -1; col_num++){

      str = mc_line.substr(0, delim_pos);
      mc_line.erase(0, delim_pos + 1);

      while (!mc_line.empty() && mc_line.at(0) == ' ')  //get rid of white space between columns
        mc_line.erase(0, 1);

      if (col_num == 9) str = to_string(new_p0[i]);
      else if (col_num == 10) str = to_string(new_p1[i]);

      writeFile << str << setw(15);
    }
    writeFile << mc_line << endl;
  }
  mc_file.close();
  writeFile.close();

//Style//

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(0.9);

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  //tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  tdrStyle->SetPaperSize(20.,20.);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->SetOptStat(0);

  tdrStyle->cd();

//End Style//

  TGraphErrors* graph = new TGraphErrors(mc_lines, eta, sf, 0, 0);
  TGraphErrors* g_high = new TGraphErrors(mc_lines, eta, sf_high, 0, 0);
  TGraphErrors* g_low = new TGraphErrors(mc_lines, eta, sf_low, 0, 0);  

  TCanvas* c = new TCanvas("c", "c", 600, 600);
  TH1F* h = new TH1F("h", "h", 100, -5, 5);
  float topY = 1.2;

  h->GetXaxis()->SetTitle("#eta");
  h->GetYaxis()->SetTitle("Scale Factor");
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetRangeUser(0.7, topY);
  h->Draw();

  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(1);
  graph->Draw("psame");

  g_high->SetMarkerStyle(2);
  g_high->SetMarkerColor(kRed);
  g_high->Draw("psame");

  g_low->SetMarkerStyle(5);
  g_low->SetMarkerColor(kBlue);
  g_low->Draw("psame");

  TLatex text;
  text.SetTextSize(0.04);

  if (pf_type.EqualTo("all"))
    text.DrawLatex(-4.9, topY+0.003, "Run II, PF");
  else
    text.DrawLatex(-4.9, topY+0.003, "Run II, PF " + pf_type);

  //text.DrawLatex(n2-4, topY-spacing, "R = 0.4" );

  string var_string;
  if (nPU_derived) var_string = "#mu";
  else var_string = "N_{PV}";

  text.DrawLatex(-1.5, 0.88, Form( "#bf{%s = %4.1f^{ #color[2]{%4.1f}}_{ #color[4]{%4.1f}}}" , var_string.c_str(), mean_nPU, nPU_high, nPU_low ) );
  text.DrawLatex(-1.5, 0.84, Form( "#bf{#rho_{Data} = %4.1f^{ #color[2]{%4.1f}}_{ #color[4]{%4.1f}} GeV}" , rho_nominal, rho_high, rho_low ) );
  text.DrawLatex(-1.5, 0.8, Form( "#bf{#rho_{MC} = %4.1f^{ #color[2]{%4.1f}}_{ #color[4]{%4.1f}} GeV}" , mc_rho_mean, mc_rho_high, mc_rho_low ) );

  text.SetTextSize(0.035);
  text.SetTextColor(1);
  text.SetTextFont(42);
  text.DrawLatex(3, topY+0.003, "#sqrt{s} = 13 TeV");

  c->Print("l1fastjet_scalefactor_eta_PF" + pf_type + ".pdf");
}
