//Chad Harrington 12/7/2015
//EXECUTE as root -l -b -q offsetpT_stack.c

#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"

#include <vector>
using namespace std;

void setStyle();

const int ETA_BINS = 82;
double etabins[ETA_BINS+1] =
  {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65,
   -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957,
   -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0,
   0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479,
   1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013,
   4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

TH1F* getHistogram(TProfile* prof, TString name){
  TH1F* histo = new TH1F(name, name, ETA_BINS, etabins);
  int size = prof->GetNbinsX();
  for (int i=1; i<=size; i++)
    histo->SetBinContent( i, prof->GetBinContent(i) );
  return histo;
}

void offsetpT_stack(){

  TFile* mcFile = TFile::Open("MC_reweight_R4.root");
  TFile* dataFile = TFile::Open("DataD_R4.root");

  int var_choice;
  cout << "\n1) nPV\n2) nPU\n\n### Enter Variable type (number): ";
  cin >> var_choice;
  cout << endl;

  TString var_type;
  int nPoints;
  if (var_choice == 1) { var_type = "nPV"; nPoints = 50; }
  else { var_type = "nPU"; nPoints = 50; }

  int n1;
  cout << "### Enter starting number of " << var_type << " ( integer from [0:" << nPoints << ") ): ";
  cin >> n1;
  cout << endl;
  int n2;
  cout << "### Enter ending number of " << var_type << " ( integer from [" << n1 << ":" << nPoints << ") ): ";
  cin >> n2;
  cout << endl;

  if (n1 < 0 || n1 > n2)
    n1 = 10;
  if (n2 >= nPoints || n2 < n1)
    n2 = 35;

  int simple_weight = -1;
  if (n1 != n2){
    cout << "1) Simple Weighting\n2) Regular Weighting\n\n### Enter weighting option: ";
    cin >> simple_weight;
    cout << endl;
  }

  enum Id{
    gamma=0, egamma_HF, h0, h_HF, chs, h, numId, all, hf_dep
  };
  TString ids[] = {"gamma", "egamma_HF", "h0", "h_HF", "chs", "h"};

  int pf_choice;
  cout << "1) All\n2) Photons\n3) EM deposits\n4) Neutral Hadrons\n5) Hadronic Deposits\n"
       << "6) Unassociated Charged Hadrons\n7) Associated Charged Hadrons\n8) HF Deposits\n\n### Enter PF Particle type: ";
  cin >> pf_choice;
  cout << endl;

  if (pf_choice == 1) pf_choice = all;
  else if (gamma <= pf_choice-2 && pf_choice-2 < numId) pf_choice -= 2;
  else if (pf_choice == numId+2) pf_choice = hf_dep;
  else pf_choice = all;

  char plot_ratio = 'n';
  cout << "### Include ratio plot (y/n)? ";
  cin >> plot_ratio;
  cout << endl;

  vector<TH1F*> v_MC (numId);
  vector<TH1F*> v_Data (numId);

  TString hname;
  for (int i=0; i<numId; i++){

    hname = "p_cone_eta_totalpT_" + var_type + Form("%i_PF", n1) + ids[i];

    v_MC[i] = getHistogram( (TProfile*) mcFile->FindObjectAny(hname), ids[i]+"MC" );
    v_Data[i] = getHistogram( (TProfile*) dataFile->FindObjectAny(hname), ids[i]+"Data" );
  }

  //Weights//

  if (n1 == n2){
    for (int i=0; i<numId; i++){
      v_MC[i]->Scale( 1.0 / n1 );
      v_Data[i]->Scale( 1.0 / n1 );
    }
  }
  else if (simple_weight == 1){
    for (int i=0; i<numId; i++){

      v_MC[i]->Scale(-1);
      v_Data[i]->Scale(-1);

      hname = "p_cone_eta_totalpT_" + var_type + Form("%i_PF", n2) + ids[i];
      v_MC[i]->Add( getHistogram( (TProfile*) mcFile->FindObjectAny(hname), ids[i]+"MC2" ) );
      v_Data[i]->Add( getHistogram( (TProfile*) dataFile->FindObjectAny(hname), ids[i]+"Data2" ) );

      v_MC[i]->Scale( 1.0 / (n2-n1) );
      v_Data[i]->Scale( 1.0 / (n2-n1) );
    }
  }
  else{

    TH1F* h_weight_MC = (TH1F*) mcFile->Get(var_type);
    TH1F* h_weight_Data = (TH1F*) dataFile->Get(var_type);

    if (h_weight_Data->GetXaxis()->GetBinWidth(1) == 0.5){
      h_weight_MC->Rebin();
      h_weight_Data->Rebin();
    }

    h_weight_MC->Scale( 1 / h_weight_MC->Integral() );
    h_weight_Data->Scale( 1 / h_weight_Data->Integral() );

    for (int i=0; i<numId; i++){

      double total_weight_MC = 0;
      double total_weight_Data = 0;

      TH1F* sum_MC = new TH1F("sum_MC", "sum_MC", ETA_BINS, etabins);
      TH1F* sum_Data = new TH1F("sum_Data", "sum_Data", ETA_BINS, etabins);

      for (int j=n1+1; j<=n2; j++){

        hname = "p_cone_eta_totalpT_" + var_type + Form("%i_PF", j) + ids[i];

        TH1F* temp_MC = getHistogram( (TProfile*) mcFile->FindObjectAny(hname), "temp_MC");
        TH1F* temp_Data = getHistogram( (TProfile*) dataFile->FindObjectAny(hname), "temp_Data");

        temp_MC->Add(v_MC[i], -1);
        temp_Data->Add(v_Data[i], -1);

        double weight_MC = h_weight_MC->GetBinContent( h_weight_MC->FindBin(j) );
        double weight_Data = h_weight_Data->GetBinContent( h_weight_Data->FindBin(j) );

        temp_MC->Scale( weight_MC / (j-n1) );
        temp_Data->Scale( weight_Data / (j-n1) );

        sum_MC->Add(temp_MC);
        sum_Data->Add(temp_Data);

        delete temp_MC;
        delete temp_Data;

        total_weight_MC += weight_MC;
        total_weight_Data += weight_Data;
      }
      v_MC[i] = (TH1F*) sum_MC->Clone("v_MC[i]");
      v_Data[i] = (TH1F*) sum_Data->Clone("v_Data[i]");

      delete sum_MC;
      delete sum_Data;

      v_MC[i]->Scale( 1 / total_weight_MC );
      v_Data[i]->Scale( 1 /total_weight_Data );
    }
  }

  //Make PF all and PF chs for ratio plot before alterations//

  TH1F* all_Data = (TH1F*) v_Data[h]->Clone("all_Data");
  TH1F* all_MC = (TH1F*) v_MC[h]->Clone("all_MC");

  for (int i=0; i<chs; i++){  //don't include chs
    all_Data->Add(v_Data[i]);
    all_MC->Add(v_MC[i]);
  }

  TH1F* chs_Data = (TH1F*) v_Data[chs]->Clone("chs_Data");
  TH1F* chs_MC = (TH1F*) v_MC[chs]->Clone("chs_MC");

  //Unassociated Charged Hadrons//

  for (int i=0; i<chs; i++){
    v_MC[chs]->Add(v_MC[i], -1);
    v_Data[chs]->Add(v_Data[i], -1);
  }

  //Associated Charged Hadrons//

  v_MC[h]->Add(v_MC[chs], -1);
  v_Data[h]->Add(v_Data[chs], -1);

  //Draw Markers for EM Deposits and Hadronic Deposits in two separate regions//

  TH1F* EM_clone = (TH1F*) v_Data[egamma_HF]->Clone("EM_clone");
  TH1F* Had_clone = (TH1F*) v_Data[h_HF]->Clone("Had_clone");

  //Stacks//

  THStack* mcStack = new THStack();
  THStack* dataStack = new THStack();
  THStack* cloneStack = new THStack();

  for (int i=0; i<numId; i++){
    mcStack->Add(v_MC[i]);
    dataStack->Add(v_Data[i]);
  }

  cloneStack->Add(v_Data[gamma]);
  cloneStack->Add(EM_clone);
  cloneStack->Add(v_Data[h0]);
  cloneStack->Add(Had_clone);

  setStyle();

  TString yTitle = "";
  if ( var_type.EqualTo("nPV") ) yTitle = "<Offset p_{T}> / <N_{PV}> (GeV)";
  else yTitle = "<Offset p_{T}> / <#mu> (GeV)";

  TCanvas* c = new TCanvas("c", "c", 600, 600);
  gStyle->SetOptStat(0);

  float b_scale = 0.3;
  float t_scale = 1 - b_scale;

  TPad* top = new TPad("top", "top", 0, b_scale, 1, 1);
  TPad* bottom = new TPad("bottom", "bottom", 0, 0, 1, b_scale);

  if ( plot_ratio == 'y' ){
    top->SetTopMargin(0.05);
    top->SetBottomMargin(0.05);
    top->Draw();
    bottom->SetBottomMargin(0.35);
    bottom->SetGridy();
    bottom->Draw();
    top->cd();
  }

  TH1F* h1 = new TH1F("h1", "h1", ETA_BINS, etabins);
  h1->GetYaxis()->SetTitle(yTitle);

  if ( plot_ratio == 'y' ) {
    h1->GetXaxis()->SetTickLength(0.03/t_scale);
    h1->GetXaxis()->SetLabelSize(0);
    h1->GetYaxis()->SetTitleSize(0.06/t_scale);
    h1->GetYaxis()->SetTitleOffset(0.8);
    h1->GetYaxis()->SetLabelSize(0.05/t_scale);
  }
  else
    h1->GetXaxis()->SetTitle("#eta");

  h1->Draw();

  v_MC[gamma]    ->SetMarkerStyle(kMultiply);
  v_MC[egamma_HF]->SetMarkerStyle(kOpenStar);
  v_MC[h0]       ->SetMarkerStyle(kOpenDiamond);
  v_MC[h_HF]     ->SetMarkerStyle(kOpenTriangleUp);
  v_MC[chs]      ->SetMarkerStyle(kOpenCircle);
  v_MC[h]        ->SetMarkerStyle(kOpenCircle);

  v_Data[gamma]    ->SetMarkerStyle(kMultiply);
  v_Data[egamma_HF]->SetMarkerStyle(kOpenStar);
  EM_clone         ->SetMarkerStyle(kOpenStar);
  v_Data[h0]       ->SetMarkerStyle(kOpenDiamond);
  v_Data[h_HF]     ->SetMarkerStyle(kOpenTriangleUp);
  Had_clone        ->SetMarkerStyle(kOpenTriangleUp);
  v_Data[chs]      ->SetMarkerStyle(kOpenCircle);
  v_Data[h]        ->SetMarkerStyle(kOpenCircle);

  v_MC[gamma]    ->SetFillColor(kBlue);
  v_MC[egamma_HF]->SetFillColor(kViolet+2);
  v_MC[h0]       ->SetFillColor(kGreen);
  v_MC[h_HF]     ->SetFillColor(kPink+6);
  v_MC[chs]      ->SetFillColor(kRed-9);
  v_MC[h]        ->SetFillColor(kRed);

  v_MC[gamma]    ->SetLineColor(kBlack);
  v_MC[egamma_HF]->SetLineColor(kBlack);
  v_MC[h0]       ->SetLineColor(kBlack);
  v_MC[h_HF]     ->SetLineColor(kBlack);
  v_MC[chs]      ->SetLineColor(kBlack);
  v_MC[h]        ->SetLineColor(kBlack);

  //CMS and Lumi Text//

  TLatex text;
  text.SetNDC();

  //text.SetTextSize(0.045);
  //text.SetTextFont(42);
  //text.DrawLatex(2, 1.41, "19.7 fb^{-1} (8 TeV)");

  if (pf_choice == all){
    float topY = 0.7;

    v_Data[gamma]    ->SetAxisRange(-3.4,3.3);
    v_Data[egamma_HF]->SetAxisRange(-5,-2.6);
    EM_clone         ->SetAxisRange(2.6,5);
    v_Data[h_HF]     ->SetAxisRange(-5,-2.6);
    Had_clone        ->SetAxisRange(2.6,5);
    v_Data[h0]       ->SetAxisRange(-3.5,3.5);
    v_Data[chs]      ->SetAxisRange(-3,3);
    v_Data[h]        ->SetAxisRange(-3,3);

    h1->GetYaxis()->SetTitleOffset(1.1);
    h1->GetYaxis()->SetRangeUser(0, topY);
    mcStack->Draw("same");
    dataStack->Draw("samep");
    cloneStack->Draw("samep");

    TLegend* leg = new TLegend(.37,.6,.8,.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->SetHeader("#bf{Markers: Data, Histograms: MC}");
    leg->AddEntry(v_MC[gamma],"Photons","PF");
    leg->AddEntry(v_MC[egamma_HF],"EM Deposits","PF");
    leg->AddEntry(v_MC[h0],"Neutral Hadrons","PF");
    leg->AddEntry(v_MC[h_HF],"Hadronic Deposits","PF");
    leg->AddEntry(v_MC[chs],"Unassoc. Charged Hadrons","PF");
    leg->AddEntry(v_MC[h],"Assoc. Charged Hadrons","PF");
    leg->Draw();

    text.SetTextSize(0.065);
    text.SetTextFont(61);
    text.DrawLatex(0.2, 0.85, "CMS");

    //text.SetTextSize(0.055);
    //text.SetTextFont(52);
    //text.DrawLatex(-4.5, 0.73, "Preliminary");

    text.SetTextSize(0.045);
    text.SetTextFont(42);
    text.DrawLatex(0.6, 0.96, "Run 2015D - 2.1 fb^{-1} (13 TeV)");
    //text.DrawLatex(-4.5, 0.5, "R = 0.4");

    gPad->RedrawAxis();

    if (plot_ratio == 'y'){
      h1->GetYaxis()->SetTitleOffset(0.8);

      bottom->cd();
      TH1F* h2 = new TH1F("h2", "h2", ETA_BINS, etabins);

      chs_Data->Divide(chs_MC);
      all_Data->Divide(all_MC);

      h2->GetXaxis()->SetLabelSize(0.05/b_scale);
      h2->GetXaxis()->SetTickLength(0.03/b_scale);
      h2->GetXaxis()->SetTitleSize(0.06/b_scale);
      h2->GetXaxis()->SetTitleOffset(0.75);
      h2->GetXaxis()->SetTitle("#eta");
      h2->GetYaxis()->SetRangeUser(0.8, 1.2);
      h2->GetYaxis()->SetNdivisions(5, 3, 0);
      h2->GetYaxis()->SetLabelSize(0.05/b_scale);
      h2->GetYaxis()->SetTitle("Data/MC");
      h2->GetYaxis()->SetTitleSize(0.05/b_scale);
      h2->GetYaxis()->SetTitleOffset(0.4);

      chs_Data->SetMarkerStyle(24);
      chs_Data->SetMarkerColor(2);
      all_Data->SetMarkerStyle(24);
      h2->Draw();
      chs_Data->Draw("sameP");
      all_Data->Draw("sameP");

      TLegend* leg = new TLegend(.55,.5,.65,.65);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.06);
      leg->SetTextFont(42);
      leg->AddEntry(all_Data,"PF","P");
      leg->AddEntry(chs_Data,"PF chs","P");
      leg->Draw();
    }
  }
  else {
    if (plot_ratio == 'y') h1->GetYaxis()->SetTitleOffset(0.8);
    else h1->GetYaxis()->SetTitleOffset(1.1);

    h1->GetYaxis()->SetRangeUser(0, 0.5);
    TH1F* hist_MC;
    TH1F* hist_Data;

    if (pf_choice != hf_dep){
      hist_MC = v_MC[pf_choice];
      hist_Data = v_Data[pf_choice];
    }
    else{
      //HF = EM Deposits + Hadronic Deposits
      hist_MC = (TH1F*) v_MC[egamma_HF]->Clone("hist_MC");
      hist_MC->Add( v_MC[h_HF] );
      hist_Data = (TH1F*) v_Data[egamma_HF]->Clone("hist_Data");
      hist_Data->Add( v_Data[h_HF] );
    }

    TString title;
    if (pf_choice == gamma) { title = "Photons"; hist_Data->SetAxisRange(-3.2, 3.2); }
    else if (pf_choice == egamma_HF) title = "EM Deposits";
    else if (pf_choice == h0) { title = "Neutral Hadrons"; hist_Data->SetAxisRange(-3.2, 3.2); }
    else if (pf_choice == h_HF) title = "Hadronic Deposits";
    else if (pf_choice == chs) { title = "Unassoc. Charged Hadrons"; hist_Data->SetAxisRange(-2.9, 2.9); }
    else if (pf_choice == h) { title = "Assoc. Charged Hadrons"; hist_Data->SetAxisRange(-2.9, 2.9); }
    else title = "HF Deposits";

    hist_MC->Draw("same");
    hist_Data->Draw("sameP");

    TLegend* leg = new TLegend(.7,.7,.9,.8);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.05);
    leg->SetTextFont(42);
    leg->AddEntry(hist_Data,"Data","P");
    leg->AddEntry(hist_MC,"MC","F");
    leg->Draw();

    text.SetTextSize(0.035/t_scale);
    text.SetTextFont(61);
    text.DrawLatex(0.17, 0.96, title);

    text.SetTextSize(0.065);
    text.SetTextFont(61);
    text.DrawLatex(0.25, 0.8, "CMS");

    //text.SetTextSize(0.055);
    //text.SetTextFont(52);
    //text.DrawLatex(-4.5, 0.73, "Preliminary");

    text.SetTextSize(0.045);
    text.SetTextFont(42);
    if (plot_ratio == 'y') text.DrawLatex(0.6, 0.96, "Run 2015D - 2.1 fb^{-1} (13 TeV)");
    else text.DrawLatex(2.5, 0.502, "#sqrt{s} = 13 TeV");
    //text.DrawLatex(-4, 0.4, "R = 0.4");

    //text.SetTextSize(0.04/t_scale);
    //text.SetTextFont(61);
    //text.DrawLatex(-4.5, 0.72, "CMS");
      
    //text.SetTextSize(0.035/t_scale);
    //text.SetTextFont(52);
    //text.DrawLatex(-4.5, 0.67, "Preliminary");

    gPad->RedrawAxis();

    if (plot_ratio == 'y'){

      bottom->cd();
      TH1F* h2 = new TH1F("h2", "h2", ETA_BINS, etabins);

      TH1F* ratio_MC = (TH1F*) hist_MC->Clone("ratio_MC");
      TH1F* ratio_Data = (TH1F*) hist_Data->Clone("ratio_Data");
      ratio_Data->Divide(ratio_MC);

      if (pf_choice == gamma) ratio_Data->SetAxisRange(-3.2, 3.2);
      //else if (pf_choice == egamma_HF) //em deposits
      else if (pf_choice == h0) ratio_Data->SetAxisRange(-3.2, 3.2);
      //else if (pf_choice == h_HF) //hadronic deposits
      else if (pf_choice == chs) ratio_Data->SetAxisRange(-2.9, 2.9);
      else if (pf_choice == h) ratio_Data->SetAxisRange(-2.9, 2.9);
      else if (pf_choice == hf_dep) ratio_Data->SetAxisRange(-3, 3);
      //else

      h2->GetXaxis()->SetLabelSize(0.05/b_scale);
      h2->GetXaxis()->SetTickLength(0.03/b_scale);
      h2->GetXaxis()->SetTitleSize(0.06/b_scale);
      h2->GetXaxis()->SetTitleOffset(0.75);
      h2->GetXaxis()->SetTitle("#eta");
      h2->GetYaxis()->SetRangeUser(0, 2);
      h2->GetYaxis()->SetNdivisions(5, 3, 0);
      h2->GetYaxis()->SetLabelSize(0.05/b_scale);
      h2->GetYaxis()->SetTitle("Data/MC");
      h2->GetYaxis()->SetTitleSize(0.05/b_scale);
      h2->GetYaxis()->SetTitleOffset(0.4);

      ratio_Data->SetMarkerStyle(24);
      h2->Draw();
      ratio_Data->Draw("sameP");
    }
  }

  c->Print("./plots/" + var_type + "/stack_" + var_type + Form("%i-%i_%i", n1, n2, pf_choice) + ".pdf");
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

  tdrStyle->cd();
}
