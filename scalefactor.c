//Chad Harrington 7/28/2015, root -l -b -q scalefactor.c

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"

using namespace std;

void scalefactor(){

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

  //Start Program//

  int Rlabel = 4;
  TString pf_type = "chs";
  bool nPU_derived = false;

  ifstream data_file("./plots/indirectRho/" + pf_type + Form("/RunII_V0_Data_RC_AK%iPF", Rlabel) + pf_type + ".txt");
  ifstream mc_file("./plots/indirectRho/" + pf_type + Form("/RunII_V0_MC_RC_AK%iPF", Rlabel) + pf_type + ".txt");
  string data_line, mc_line;

  //read first line
  getline(data_file, data_line);
  getline(mc_file, mc_line);

  double eta1[100], eta2[100], data_p0[100], data_p1[100], data_p2[100];
  double mc_p0[100], mc_p1[100], mc_p2[100];

  for (int i=0; getline(data_file,data_line); i++){

    string str;
    int delim_pos;

    while (data_line.at(0) == ' ')  //positive values of eta have space in front of them
      data_line.erase(0, 1);

    //loop over strings in data line
    for (int col_num=0; (delim_pos = data_line.find(' ')) != -1; col_num++){

      str = data_line.substr(0, delim_pos);
      data_line.erase(0, delim_pos + 1);

      while (data_line.at(0) == ' ')  //get rid of white space between columns
        data_line.erase(0, 1);

      if (col_num == 0) eta1[i] = stod(str);
      else if (col_num == 1) eta2[i] = stod(str);
      else if (col_num == 9) data_p0[i] = stod(str);
      else if (col_num == 10) data_p1[i] = stod(str);
    }
    //last column after loop
    data_p2[i] = stod(data_line);

    //mc line
    getline(mc_file, mc_line);

    while (mc_line.at(0) == ' ')  //positive values of eta have space in front of them
      mc_line.erase(0, 1);

    //loop over strings in mc line
    for (int col_num=0; (delim_pos = mc_line.find(' ')) != -1; col_num++){

      str = mc_line.substr(0, delim_pos);
      mc_line.erase(0, delim_pos + 1);

      while (mc_line.at(0) == ' ')  //get rid of white space between columns
        mc_line.erase(0, 1);

      if (col_num == 9) mc_p0[i] = stod(str);
      else if (col_num == 10) mc_p1[i] = stod(str);
    }
    //last column after loop
    mc_p2[i] = stod(mc_line);
  }
  data_file.close();
  mc_file.close();

  TFile* data_root = TFile::Open( Form("RunII_DataD_2509_R%i.root", Rlabel) );
  TFile* mc_root = TFile::Open( Form("RunII_MCD_2509_R%i.root", Rlabel) );

  TString hname;
  if (nPU_derived) hname = "p_rho_nPU";
  else hname = "p_rho_nPV";

  TProfile* data_rho_nPU = (TProfile*) data_root->Get(hname);

  for (int i=1; i<=100; i++)
    cout << i << "\t" << data_rho_nPU->GetBinCenter(i) << "\t" << data_rho_nPU->GetBinContent(i) << endl;

  int low_bin = 1, high_bin = 29;
  cout << low_bin << "\t" << high_bin << endl;

  TProfile* mc_rho_nPU = (TProfile*) mc_root->Get(hname);

  ofstream writeFile(Form("./plots/scalefactor/DataMcSF_L1FastJet_AK%iPF", Rlabel) + pf_type + ".txt");
  writeFile << "{1   JetEta   1   Rho   [0]+[1]*x+[2]*pow(x,2)   Data/MC   L1FastJet}" << endl;

  TF1* fit = new TF1("fit", "1++x++x*x");
  fit->SetLineColor(1);
  fit->SetLineWidth(2);

  float rho_start = 3;
  float rho_end = 12;

  int size;
  if (nPU_derived) size = (rho_end-rho_start+0.5)*2;
  else size = rho_end-rho_start+1;
  const int rho_size = size;

  for (int i=0; i<100; i++){

    double rho[rho_size], scale_factor[rho_size], sf_error[rho_size];

    for (int j=0; j<rho_size; j++){

      if (nPU_derived) rho[j] = rho_start+j/2.0;
      else rho[j] = rho_start+j;

      double data_offset = data_p0[i] + data_p1[i]*rho[j] + data_p2[i]*rho[j]*rho[j];

      int data_mu_bin = 0;
      data_rho_nPU->GetBinWithContent(rho[j], data_mu_bin, low_bin, high_bin, 1);
      double data_mu = data_rho_nPU->GetBinCenter( data_mu_bin );

      double mc_rho = mc_rho_nPU->GetBinContent( mc_rho_nPU->FindBin( data_mu ) );
      double mc_offset = mc_p0[i] + mc_p1[i]*mc_rho + mc_p2[i]*mc_rho*mc_rho;

      scale_factor[j] = data_offset / mc_offset;
      sf_error[j] = 0.05*scale_factor[j];
    }

    TGraphErrors* graph = new TGraphErrors(rho_size, rho, scale_factor, 0, sf_error);
    graph->Fit(fit);

    writeFile << setw(4) << eta1[i] << setw(6) << eta2[i] << setw(4) << 5 << setw(4) << 0 << setw(6) << 200
              << setw(15) << fit->GetParameter(0) << setw(15) << fit->GetParameter(1) << setw(15) << fit->GetParameter(2) << endl;

    TCanvas* c = new TCanvas("c", "c", 600, 600);
    TH1F* h = new TH1F("h", "h", rho_size-1, rho_start, rho_end);
    int topY = 2;

    h->GetXaxis()->SetTitle("#rho_{Data} (GeV)");
    h->GetYaxis()->SetTitle("Scale Factor");
    h->GetYaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetRangeUser(0, topY);
    h->Draw();

    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(1);
    graph->Draw("psame");

    float spacing = 0.1;
    float legOff = rho_start+0.05*(rho_end-rho_start);

    TLatex text;
    text.SetTextSize(0.04);

    if (pf_type.EqualTo("all")) 
      text.DrawLatex(rho_start+0.025*(rho_end-rho_start), topY+0.01, Form("Run II, PF %4.1f #leq #eta #leq %4.1f", eta1[i], eta2[i]) );
    else
      text.DrawLatex(rho_start+0.025*(rho_end-rho_start), topY+0.01, "Run II, PF " + pf_type + Form(" %4.1f #leq #eta #leq %4.1f", eta1[i], eta2[i]) );

    //text.DrawLatex(n2-4, topY-spacing, "R = 0.4" );

    text.DrawLatex(legOff, topY-0.1-spacing, Form("#chi^{2}/ndof = %4.2f/%i", fit->GetChisquare(), fit->GetNDF() ) );
    text.DrawLatex(legOff, topY-0.1-2*spacing, Form("p0 = %4.3f #pm %4.3f", fit->GetParameter(0), fit->GetParError(0) ) );
    text.DrawLatex(legOff, topY-0.1-3*spacing, Form("p1 = %4.3f #pm %4.3f", fit->GetParameter(1), fit->GetParError(1) ) );
    text.DrawLatex(legOff, topY-0.1-4*spacing, Form("p2 = %4.4f #pm %4.4f", fit->GetParameter(2), fit->GetParError(2) ) );

    text.SetTextSize(0.035);
    text.SetTextColor(1);
    text.SetTextFont(42);
    text.DrawLatex(rho_start+0.8*(rho_end-rho_start), topY+0.01, "#sqrt{s} = 13 TeV");

    c->Print("./plots/scalefactor/scalefactor_PF" + pf_type + Form("_eta%4.2f.pdf", eta1[i]+0.05) );
    delete h;
    delete c;
  }
  writeFile.close();
}
