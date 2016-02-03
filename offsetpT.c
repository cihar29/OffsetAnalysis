//Chad Harrington 1/20/2015
//EXECUTE as root -l -b -q offsetpT.c

#include "TFile.h"
#include <vector>
#include "TProfile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

void setStyle();

void offsetpT(int n1=1, int n2=25, float topX=15, float topY=30){

  setStyle();
        
  const double R = 0.8;
  int Rlabel = R*10;

  TFile* mcFile = TFile::Open( Form("MC_reweight_R%i.root", Rlabel) );
  TFile* dataFile = TFile::Open( Form("DataD_R%i.root", Rlabel) );

  TString var_type = "nPV"; // nPU, nPV
  bool isIndirect = true;   // indirectRho
  bool rhoCentral = false;

  //int n1 = 1, n2 = 25;
                
  const int nPoints = n2-n1;

  int pf_choice;
  cout << "\n1) All\n2) Charged Hadron\n3) Electron\n4) Muon\n5) Photon\n6) Neutral Hadron\n7) HF Tower Hadron\n"
       << "8) HF Tower EM Particle\n9) Charged Hadron Subtraction\n\n### Enter PF Particle type: ";
  cin >> pf_choice;
  cout << endl;

  TString pf_type = "all";    //default choice
  if (pf_choice == 2) pf_type = "h";
  else if (pf_choice == 3) pf_type = "e";
  else if (pf_choice == 4) pf_type = "mu";
  else if (pf_choice == 5) pf_type = "gamma";
  else if (pf_choice == 6) pf_type = "h0";
  else if (pf_choice == 7) pf_type = "h_HF";
  else if (pf_choice == 8) pf_type = "egamma_HF";
  else if (pf_choice == 9) pf_type = "chs";

  vector<TProfile*> mcProfiles;
  vector<TProfile*> dataProfiles;
	
  for (int n=0; n<nPoints; n++){
    TString hname = "p_cone_eta_totalpT_" + var_type + Form("%i_PF", n1+n) + pf_type;

    mcProfiles.push_back( (TProfile*) mcFile->FindObjectAny(hname) );
    dataProfiles.push_back( (TProfile*) dataFile->FindObjectAny(hname) );
  }
	
  TProfile* NPUtoRho_mc;
  TProfile* NPUtoRho_data;

  if (isIndirect){

    TString hname;
    if (var_type.EqualTo("nPU")) {
      if (rhoCentral) hname = "p_rhoCentral0_nPU";
      else hname = "p_rho_nPU";
    }
    else {
      if (rhoCentral) hname = "p_rhoCentral0_nPV";
      else hname = "p_rho_nPV";
    }

    NPUtoRho_mc = (TProfile*) mcFile->Get(hname);
    NPUtoRho_data = (TProfile*) dataFile->Get(hname);

    if (NPUtoRho_data->GetXaxis()->GetBinWidth(1) == 0.5){
      NPUtoRho_mc->Rebin();
      NPUtoRho_data->Rebin();
    }
    var_type = "indirectRho";
  }

  ofstream writeMC("./plots/" + var_type + "/" + pf_type + Form("/Fall15_25nsV1_MC_L1RC_AK%iPF", Rlabel) + pf_type + ".txt");
  ofstream writeData("./plots/" + var_type + "/" + pf_type + Form("/Fall15_25nsV1_DATA_L1RC_AK%iPF", Rlabel) + pf_type + ".txt");

  TString header;
  if ( var_type.EqualTo("nPV") || var_type.EqualTo("nPU") )
    header = "{1\tJetEta\t3\tJetPt\tJetA\t" + var_type + "\t\tmax(0.0001,1-y*([0]+[1]*(z-1)+[2]*pow(z-1,2))/x)\tCorrection L1Offset}";
  else
    header = "{1         JetEta              3          JetPt           JetA            Rho               max(0.0001,1-y*([0]+[1]*(z-1.519)+[2]*pow(z-1.519,2))/x)     Correction      L1FastJet}";
	
  writeMC << header << endl;
  writeData << header << endl;
	
  TF1* f_mc = new TF1("f_mc", "1++x++x*x");
  TF1* f_data = new TF1("f_data", "1++x++x*x");
  f_mc->SetLineColor(2);
  f_mc->SetLineWidth(2);
  f_data->SetLineColor(1);
  f_data->SetLineWidth(2);

  const int ETA_BINS = 82;
  float etabins[ETA_BINS+1] =
  {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65,
   -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957,
   -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0,
   0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479,
   1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013,
   4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

  for (int i=0; i<ETA_BINS; i++) {
    vector<double> mc_x, data_x, mc_y, data_y, mc_error, data_error;
		
    for (int n=0; n<nPoints; n++){
//      if (n == 5 || n == 6) continue;

      double mcX = n1+n+0.5;
      double dataX = n1+n+0.5;

      if (isIndirect){
        mcX = NPUtoRho_mc->GetBinContent( NPUtoRho_mc->FindBin(mcX) );
        dataX = NPUtoRho_data->GetBinContent( NPUtoRho_data->FindBin(dataX) );
      }
      mc_x.push_back( mcX );
      data_x.push_back( dataX );

      mc_y.push_back( mcProfiles[n]->GetBinContent(i+1) );
      data_y.push_back( dataProfiles[n]->GetBinContent(i+1) );

      mc_error.push_back( 0.02 * mc_y.back() );
      data_error.push_back( 0.01 * data_y.back() );
    }
    TGraphErrors* mcGraph = new TGraphErrors(mc_x.size(), &mc_x[0], &mc_y[0], 0, &mc_error[0]);
    TGraphErrors* dataGraph = new TGraphErrors(data_x.size(), &data_x[0], &data_y[0], 0, &data_error[0]);

    mcGraph->Fit(f_mc, "Q");
    dataGraph->Fit(f_data, "Q");

    const double PI = 3.14159265359;
    double area = PI*R*R;

    writeMC << etabins[i] << setw(8) << etabins[i+1] << setw(8) << 9 << setw(6) << 1 << setw(6) << 3500 << setw(6)
            << 0 << setw(6) << 10 << setw(6) << 0 << setw(6) << 200
            << setw(15) << f_mc->GetParameter(0)/area << setw(15) << f_mc->GetParameter(1)/area
            << setw(15) << f_mc->GetParameter(2)/area << endl;
    writeData << etabins[i] << setw(8) << etabins[i+1] << setw(8) << 9 << setw(6) << 1 << setw(6) << 3500 << setw(6)
              << 0 << setw(6) << 10 << setw(6) << 0 << setw(6) << 200
              << setw(15) << f_data->GetParameter(0)/area << setw(15) << f_data->GetParameter(1)/area
              << setw(15) << f_data->GetParameter(2)/area << endl;

    TString xTitle = "";
    if ( var_type.EqualTo("nPV") ) xTitle = "N_{PV}";
    else if ( var_type.EqualTo("nPU") ) xTitle = "#mu";
    else if ( var_type.EqualTo("indirectRho") ) xTitle = "<#rho> (GeV)";

    TCanvas* c = new TCanvas("c", "c", 600, 600);
    //float topX = 15;
    TH1F* h = new TH1F("h", "h", 100, 0, topX);
    //float topY = 30;

    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle("Offset p_{T} (GeV)");
    h->GetYaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetRangeUser(0, topY);

    h->Draw();
    dataGraph->SetMarkerStyle(20);
    dataGraph->SetMarkerColor(1);
    dataGraph->Draw("Psame");
    mcGraph->SetMarkerStyle(24);
    mcGraph->SetMarkerColor(2);
    mcGraph->SetLineColor(2);
    mcGraph->Draw("Psame");

    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.04);

    if (pf_type.EqualTo("all"))
      text.DrawLatex(0.17, 0.96, Form("AK%i PF %4.3f #leq #eta #leq %4.3f", Rlabel, etabins[i], etabins[i+1]) );
    else
      text.DrawLatex(0.17, 0.96, Form("AK%i PF%s %4.3f #leq #eta #leq %4.3f", Rlabel, pf_type.Data(), etabins[i], etabins[i+1]) );

    text.DrawLatex(0.2, 0.88, "Data");
    text.DrawLatex(0.2, 0.84, Form("#chi^{2}/ndof = %4.2f/%i", f_data->GetChisquare(), f_data->GetNDF() ) );
    text.DrawLatex(0.2, 0.8, Form("p0 = %4.3f #pm %4.3f", f_data->GetParameter(0), f_data->GetParError(0) ) );
    text.DrawLatex(0.2, 0.76, Form("p1 = %4.3f #pm %4.3f", f_data->GetParameter(1), f_data->GetParError(1) ) );
    text.DrawLatex(0.2, 0.72, Form("p2 = %4.4f #pm %4.4f", f_data->GetParameter(2), f_data->GetParError(2) ) );
    text.SetTextColor(2);
    text.DrawLatex(0.2, 0.64, "MC");
    text.DrawLatex(0.2, 0.6, Form("#chi^{2}/ndof = %4.2f/%i", f_mc->GetChisquare(), f_mc->GetNDF() ) );
    text.DrawLatex(0.2, 0.56, Form("p0 = %4.3f #pm %4.3f", f_mc->GetParameter(0), f_mc->GetParError(0) ) );
    text.DrawLatex(0.2, 0.52, Form("p1 = %4.3f #pm %4.3f", f_mc->GetParameter(1), f_mc->GetParError(1) ) );
    text.DrawLatex(0.2, 0.48, Form("p2 = %4.4f #pm %4.4f", f_mc->GetParameter(2), f_mc->GetParError(2) ) );

    text.SetTextSize(0.035);
    text.SetTextFont(42);
    text.SetTextColor(1);
    text.DrawLatex( 0.8, 0.96, "#sqrt{s} = 13 TeV" );

    cout << f_data->GetChisquare() / f_data->GetNDF() << "\t" << f_mc->GetChisquare() / f_mc->GetNDF() << endl;

    c->Print("./plots/" + var_type + "/" + pf_type + "/" + pf_type + "_pT_" + var_type + Form("_eta%4.3f", etabins[i]) + ".pdf");
    delete h;
    delete c;
  }
  writeMC.close();
  writeData.close();
}

void setStyle(){

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
}
