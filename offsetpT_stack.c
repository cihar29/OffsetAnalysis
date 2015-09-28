//Chad Harrington 2/20/2015
//EXECUTE as root -l -b -q offsetpT_stack.c

#include "TFile.h"
#include <vector>
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
using namespace std;

TH1F* getHistogram(TProfile* prof, TString name){
  TH1F* histo = new TH1F(name, name, 100, -5, 5);
  int size = prof->GetNbinsX();
  for (int i=1; i<=size; i++)
    histo->SetBinContent( i, prof->GetBinContent(i) );
  return histo;
}

void offsetpT_stack(){

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

  //By end of code we will have:
  //   pfID[] = { "gamma e mu", ..., ..., "egamma_HF", "h0", "h_HF", "chs", "h" }
  //                 0                        3         4      5       6     7   

  int pf_choice;
  cout << "1) All\n2) Photons\n3) EM deposits\n4) Neutral Hadrons\n5) Hadronic Deposits\n"
       << "6) Unassociated Charged Hadrons\n7) Associated Charged Hadrons\n8) HF Deposits\n\n### Enter PF Particle type: ";
  cin >> pf_choice;
  cout << endl;

  if (pf_choice == 1) pf_choice = -1;
  else if (pf_choice == 2) pf_choice = 0;
  else if (pf_choice != 3 && pf_choice != 4 && pf_choice != 5 &&
           pf_choice != 6 && pf_choice != 7 && pf_choice != 8) pf_choice = -1; //default option

  char plot_ratio = 'n';
  cout << "### Include ratio plot (y/n)? ";
  cin >> plot_ratio;
  cout << endl;

  TFile* mcFile = TFile::Open("RunII_MCD_2509_p865_R4.root");
  TFile* dataFile = TFile::Open("RunII_DataD_2509_p865_R4.root");

  vector<TH1F*> h_MC (8);
  vector<TH1F*> h_Data (8);

  //PF particle Histograms//

  TString pfID[] = { "gamma", "e", "mu", "egamma_HF", "h0", "h_HF", "chs", "h" };
  //ROOT doesn't like maps, so we're stuck with shoddy parallel arrays

  int size = h_MC.size();
  TString hname;

  for (int i=0; i<size; i++){

    hname = "p_cone_eta_totalpT_" + var_type + Form("%i_PF", n1) + pfID[i];

    h_MC[i] = getHistogram( (TProfile*) mcFile->FindObjectAny(hname), "h_" + pfID[i] + "_MC" );
    h_Data[i] = getHistogram( (TProfile*) dataFile->FindObjectAny(hname), "h_" + pfID[i] + "_Data" );
  }

  //Weights//

  vector<TH1F*> w_MC (size);
  vector<TH1F*> w_Data (size);

  if (n1 == n2){
    for (int j=0; j<size; j++){    //Loop over PF particles

      hname = Form("w_mc%i", j);
      w_MC[j] = new TH1F(hname, hname, 100, -5, 5);  //Easier to clone array object this way
      w_MC[j]->Add( h_MC[j] );

      hname = Form("w_data%i", j);
      w_Data[j] = new TH1F(hname, hname, 100, -5, 5);
      w_Data[j]->Add( h_Data[j] );

      w_MC[j]->Scale( 1.0 / n1 );
      w_Data[j]->Scale( 1.0 / n1 );
    }
  }

  //Simple Scaling//

  else if (simple_weight == 1){
    for (int j=0; j<size; j++){    //Loop over PF particles

      hname = Form("w_mc%i", j);
      w_MC[j] = new TH1F(hname, hname, 100, -5, 5);
      w_MC[j]->Add( h_MC[j] );

      hname = Form("w_data%i", j);
      w_Data[j] = new TH1F(hname, hname, 100, -5, 5);
      w_Data[j]->Add( h_Data[j] );

      w_MC[j]->Scale(-1);
      w_Data[j]->Scale(-1);

      hname = "p_cone_eta_totalpT_" + var_type + Form("%i_PF", n2) + pfID[j];
      w_MC[j]->Add( getHistogram( (TProfile*) mcFile->FindObjectAny(hname), Form("temp_MC%i_", n2) + pfID[j] ) );
      w_Data[j]->Add( getHistogram( (TProfile*) dataFile->FindObjectAny(hname), Form("temp_Data%i_", n2) + pfID[j] ) );

      w_MC[j]->Scale( 1.0 / (n2-n1) );
      w_Data[j]->Scale( 1.0 / (n2-n1) );
    }
  }

  //Regular Scaling//

  else{

    if ( var_type.EqualTo("nPV") ) hname = "nPV";
    else hname = "nPU";

    TH1F* h_weight_MC = (TH1F*) mcFile->Get(hname);
    TH1F* h_weight_Data = (TH1F*) dataFile->Get(hname);

    if (h_weight_Data->GetXaxis()->GetBinWidth(1) == 0.5){
      h_weight_MC->Rebin();
      h_weight_Data->Rebin();
    }

    h_weight_MC->Scale( 1 / h_weight_MC->Integral() );
    h_weight_Data->Scale( 1 / h_weight_Data->Integral() );

    for (int j=0; j<size; j++){    //Loop over PF particles

      double total_weight_MC = 0;
      double total_weight_Data = 0;

      hname = Form("w_mc%i", j);
      w_MC[j] = new TH1F(hname, hname, 100, -5, 5);
      
      hname = Form("w_data%i", j);
      w_Data[j] = new TH1F(hname, hname, 100, -5, 5);

      for (int i=n1+1; i<=n2; i++){

        hname = "p_cone_eta_totalpT_" + var_type + Form("%i_PF", i) + pfID[j];

        TH1F* temp_MC = getHistogram( (TProfile*) mcFile->FindObjectAny(hname), Form("temp_MC%i_", i) + pfID[j] );
        TH1F* temp_Data = getHistogram( (TProfile*) dataFile->FindObjectAny(hname), Form("temp_Data%i_", i) + pfID[j] );

        temp_MC->Add(h_MC[j], -1);
        temp_Data->Add(h_Data[j], -1);

        double weight_MC = h_weight_MC->GetBinContent( h_weight_MC->FindBin(i) );
        double weight_Data = h_weight_Data->GetBinContent( h_weight_Data->FindBin(i) );

        temp_MC->Scale( weight_MC / (i-n1) );
        temp_Data->Scale( weight_Data / (i-n1) );

        w_MC[j]->Add(temp_MC);
        w_Data[j]->Add(temp_Data);

        total_weight_MC += weight_MC;
        total_weight_Data += weight_Data;
      }
      w_MC[j]->Scale( 1 / total_weight_MC );
      w_Data[j]->Scale( 1 /total_weight_Data );
    }
  }

  //Make PF all and PF chs for ratio plot//

  TH1F* all_Data = (TH1F*) w_Data[7]->Clone("all_Data");
  TH1F* all_MC = (TH1F*) w_MC[7]->Clone("all_MC");

  for (int i=0; i<6; i++){
    all_Data->Add(w_Data[i]);
    all_MC->Add(w_MC[i]);
  }

  TH1F* chs_Data = (TH1F*) w_Data[6]->Clone("chs_Data");
  TH1F* chs_MC = (TH1F*) w_MC[6]->Clone("chs_MC");

  //Add together gamma e and mu//

  w_MC[0]->Add(w_MC[1]);
  w_MC[0]->Add(w_MC[2]);
  w_Data[0]->Add(w_Data[1]);
  w_Data[0]->Add(w_Data[2]);

  //Subtract Unassociated Charged Hadrons//

  for (int i=0; i<6; i++){
    w_MC[6]->Add(w_MC[i], -1);
    w_Data[6]->Add(w_Data[i], -1);
  }

  //Charged Hadrons//

  w_MC[7]->Add(w_MC[6], -1);
  w_Data[7]->Add(w_Data[6], -1);

  //Draw Markers for EM Deposits and Hadronic Deposits in two separate regions//

  TH1F* EM_clone = (TH1F*) w_Data[3]->Clone("EM_clone");
  TH1F* Had_clone = (TH1F*) w_Data[5]->Clone("Had_clone");

  //Stacks//

  THStack* mcStack = new THStack();
  THStack* dataStack = new THStack();
  THStack* cloneStack = new THStack();

  for (int i=0; i<size; i++){
    if (i != 1 && i != 2){
      mcStack->Add(w_MC[i]);
      dataStack->Add(w_Data[i]);
    }
  }

  cloneStack->Add(w_Data[0]);
  cloneStack->Add(EM_clone);
  cloneStack->Add(w_Data[4]);
  cloneStack->Add(Had_clone);

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

  TString yTitle = "";
  if ( var_type.EqualTo("nPV") ) yTitle = "<Offset p_{T}> / <N_{PV}> (GeV)";
  else yTitle = "<Offset p_{T}> / <#mu> (GeV)";

  TCanvas* c = new TCanvas("c", "Offset pT vs. eta", 600, 600);
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

  TH1* h = new TH1F("h", Form("PF Distribution %i #leq ", n1) + var_type + Form(" #leq %i", n2),100,-5,5);
  h->GetYaxis()->SetTitle(yTitle);

  if ( plot_ratio == 'y' ) {
    h->GetXaxis()->SetTickLength(0.03/t_scale);
    h->GetXaxis()->SetLabelSize(0);
    h->GetYaxis()->SetTitleSize(0.06/t_scale);
    h->GetYaxis()->SetTitleOffset(0.8);
    h->GetYaxis()->SetLabelSize(0.05/t_scale);
  }
  else
    h->GetXaxis()->SetTitle("#eta");

  h->Draw();

  //REMINDER: pfID[] = { "gamma e mu", ..., ..., "egamma_HF", "h0", "h_HF", "chs", "h" }
  //                        0                        3         4      5       6     7     

  w_MC[0]->SetMarkerStyle(kMultiply);        //For the legend
  w_MC[3]->SetMarkerStyle(kOpenStar);
  w_MC[4]->SetMarkerStyle(kOpenDiamond);
  w_MC[5]->SetMarkerStyle(kOpenTriangleUp);
  w_MC[6]->SetMarkerStyle(kOpenCircle);
  w_MC[7]->SetMarkerStyle(kOpenCircle);

  w_Data[0]->SetMarkerStyle(kMultiply);
  w_Data[3]->SetMarkerStyle(kOpenStar);
  EM_clone->SetMarkerStyle(kOpenStar);
  w_Data[4]->SetMarkerStyle(kOpenDiamond);
  w_Data[5]->SetMarkerStyle(kOpenTriangleUp);
  Had_clone->SetMarkerStyle(kOpenTriangleUp);
  w_Data[6]->SetMarkerStyle(kOpenCircle);
  w_Data[7]->SetMarkerStyle(kOpenCircle);

  w_MC[0]->SetFillColor(kBlue);
  w_MC[3]->SetFillColor(kViolet+2);
  w_MC[4]->SetFillColor(kGreen);
  w_MC[5]->SetFillColor(kPink+6);
  w_MC[6]->SetFillColor(kRed-9);
  w_MC[7]->SetFillColor(kRed);

  w_MC[0]->SetLineColor(kBlack);
  w_MC[3]->SetLineColor(kBlack);
  w_MC[4]->SetLineColor(kBlack);
  w_MC[5]->SetLineColor(kBlack);
  w_MC[6]->SetLineColor(kBlack);
  w_MC[7]->SetLineColor(kBlack);

  //CMS and Lumi Text//

  TLatex text;

  //text.SetTextSize(0.045);
  //text.SetTextFont(42);
  //text.DrawLatex(2, 1.41, "19.7 fb^{-1} (8 TeV)");

  if (pf_choice == -1){  //all PF particles

    w_Data[0]->SetAxisRange(-3.4,3.3);
    w_Data[3]->SetAxisRange(-5,-2.6);
    EM_clone->SetAxisRange(2.6,5);
    w_Data[5]->SetAxisRange(-5,-2.6);
    Had_clone->SetAxisRange(2.6,5);
    w_Data[4]->SetAxisRange(-3.5,3.5);
    w_Data[6]->SetAxisRange(-3,3);
    w_Data[7]->SetAxisRange(-3,3);

    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetRangeUser(0, 0.8);
    mcStack->Draw("same");
    dataStack->Draw("samep");
    cloneStack->Draw("samep");

    TLegend* leg = new TLegend(.4,.6,.8,.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->SetHeader("#bf{Markers: Data, Histograms: MC}");
    leg->AddEntry(w_MC[0],"Photons","PF");
    leg->AddEntry(w_MC[3],"EM Deposits","PF");
    leg->AddEntry(w_MC[4],"Neutral Hadrons","PF");
    leg->AddEntry(w_MC[5],"Hadronic Deposits","PF");
    leg->AddEntry(w_MC[6],"Unassoc. Charged Hadrons","PF");
    leg->AddEntry(w_MC[7],"Assoc. Charged Hadrons","PF");
    leg->Draw();

    //text.SetTextSize(0.06);
    //text.SetTextFont(61);
    //text.DrawLatex(-4.5, 0.8, "CMS");

    //text.SetTextSize(0.055);
    //text.SetTextFont(52);
    //text.DrawLatex(-4.5, 0.73, "Preliminary");

    text.DrawLatex(-4.8, 0.803, "Run II");

    text.SetTextSize(0.045);
    text.SetTextFont(42);
    if (plot_ratio == 'y') text.DrawLatex(3, 0.803, "#sqrt{s} = 13 TeV");
    else text.DrawLatex(2.5, 0.803, "#sqrt{s} = 13 TeV");
    //text.DrawLatex(-4.5, 0.5, "R = 0.4");

    gPad->RedrawAxis();

    if (plot_ratio == 'y'){
      h->GetYaxis()->SetTitleOffset(0.8);

      bottom->cd();
      TH1F* h2 = new TH1F("h2", "h2", 100, -5, 5);

      chs_Data->Divide(chs_MC);
      all_Data->Divide(all_MC);

      h2->GetXaxis()->SetLabelSize(0.05/b_scale);
      h2->GetXaxis()->SetTickLength(0.03/b_scale);
      h2->GetXaxis()->SetTitleSize(0.06/b_scale);
      h2->GetXaxis()->SetTitleOffset(0.75);
      h2->GetXaxis()->SetTitle("#eta");
      h2->GetYaxis()->SetRangeUser(0.8, 1.1);
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
    if (plot_ratio == 'y') h->GetYaxis()->SetTitleOffset(0.8);
    else h->GetYaxis()->SetTitleOffset(1.1);

    h->GetYaxis()->SetRangeUser(0, 0.5);
    TH1F* hist_MC;
    TH1F* hist_Data;

    if (pf_choice != 8){
      hist_MC = w_MC[pf_choice];
      hist_Data = w_Data[pf_choice];
    }
    else{
      //HF = EM Deposits + Hadronic Deposits
      hist_MC = (TH1F*) w_MC[3]->Clone("hist_MC");
      hist_MC->Add( w_MC[5] );
      hist_Data = (TH1F*) w_Data[3]->Clone("hist_Data");
      hist_Data->Add( w_Data[5] );
    }

    TString title;
    if (pf_choice == 0) { title = "Photons"; hist_Data->SetAxisRange(-3.5, 3.5); }
    else if (pf_choice == 3) title = "EM Deposits";
    else if (pf_choice == 4) { title = "Neutral Hadrons"; hist_Data->SetAxisRange(-3.5, 3.5); }
    else if (pf_choice == 5) title = "Hadronic Deposits";
    else if (pf_choice == 6) { title = "Unassoc. Charged Hadrons"; hist_Data->SetAxisRange(-3, 3); }
    else if (pf_choice == 7) { title = "Assoc. Charged Hadrons"; hist_Data->SetAxisRange(-3, 3); }
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

    text.SetTextSize(0.04/t_scale);
    text.SetTextFont(61);
    text.DrawLatex(-4, 0.4, title);

    text.DrawLatex(-4.8, 0.502, "Run II");

    text.SetTextSize(0.045);
    text.SetTextFont(42);
    if (plot_ratio == 'y') text.DrawLatex(3, 0.502, "#sqrt{s} = 13 TeV");
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
      TH1F* h2 = new TH1F("h2", "h2", 100, -5, 5);

      TH1F* ratio_MC = (TH1F*) hist_MC->Clone("ratio_MC");
      TH1F* ratio_Data = (TH1F*) hist_Data->Clone("ratio_Data");
      ratio_Data->Divide(ratio_MC);

      if (pf_choice == 0) ratio_Data->SetAxisRange(-3.5, 3.5);
      //else if (pf_choice == 3) //em deposits
      else if (pf_choice == 4) ratio_Data->SetAxisRange(-3.5, 3.5);
      //else if (pf_choice == 5) //hadronic deposits
      else if (pf_choice == 6) ratio_Data->SetAxisRange(-3, 3);
      else if (pf_choice == 7) ratio_Data->SetAxisRange(-3, 3);
      else if (pf_choice == 8) ratio_Data->SetAxisRange(-3, 3);
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

  c->Print("./plots/" + var_type + "/stack_pT_" + var_type + Form("%i-%i_%i", n1, n2, pf_choice) + ".pdf");
}
