// -*- C++ -*-
//
// Package:    offset/OffsetAnalysis
// Class:      OffsetAnalysis
// 
/**\class OffsetAnalysis OffsetAnalysis.cc offset/OffsetAnalysis/plugins/OffsetAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  charles harrington
//         Created:  Thur, 16 July 2015 16:54:04 GMT
//
//

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/Ref.h"
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include "parsePileUpJSON2.h"
#include <vector>
#include <cmath>

//root files
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom3.h>

using namespace std;

const int ETA_BINS = 82;
double etabins[ETA_BINS+1] =
  {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65,
   -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957,
   -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0,
   0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479,
   1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013,
   4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

class OffsetAnalysis : public edm::EDAnalyzer {
  public:
    explicit OffsetAnalysis(const edm::ParameterSet&);

  private:
    enum Flavor{
      all=0, h, e, mu, gamma, h0, h_HF, egamma_HF, chs, numFlavors, X //undefined
    };
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    void FillHist1D(const TString& histName, const Double_t& value, Double_t weight);
    void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, Double_t weight);
    void FillProfile(const TString& histName, const Double_t& value1, const Double_t& value2, Double_t weight);
    void FillProfile2D(const TString& histName, const Double_t& value1, const Double_t& value2, const Double_t& value3, Double_t weight);
    Flavor getFlavor(reco::PFCandidate::ParticleType id);
    double deltaPhi(double phi1, double phi2);

    //int counter;	  
    TRandom3* rand;
    TFile* root_file;
    TString ids[numFlavors];

    map<TString, TH1*> m_Histos1D;
    map<TString, TH2*> m_Histos2D;  
    map<TString, TProfile*> m_Profiles;
    map<TString, TProfile2D*> m_Profiles2D;

    TString RootFileName_;
    double coneDR_;
    int numSkip_, maxNPV_, maxNPU_, maxRho_;
    bool isMC_, reweight_;
    edm::EDGetTokenT< vector<reco::Vertex> > pvTag_;
    edm::EDGetTokenT< vector<PileupSummaryInfo> > puTag_;
    edm::EDGetTokenT< vector<reco::PFCandidate> > pfTag_;
    edm::EDGetTokenT<double> rhoTag_;
    edm::EDGetTokenT<double> rhoCentral0Tag_;
    edm::EDGetTokenT<double> rhoCentralChargedTag_;
};

OffsetAnalysis::OffsetAnalysis(const edm::ParameterSet& iConfig)
{
  numSkip_ = iConfig.getParameter<int> ("numSkip");
  RootFileName_ = iConfig.getParameter<string> ("RootFileName");
  coneDR_ = iConfig.getParameter<double> ("coneDR");
  maxNPV_ = iConfig.getParameter<int> ("maxNPV");
  maxNPU_ = iConfig.getParameter<int> ("maxNPU");
  maxRho_ = iConfig.getParameter<int> ("maxRho");
  isMC_ = iConfig.getParameter<bool> ("isMC");
  reweight_ = iConfig.getParameter<bool> ("reweight");
  pvTag_ = consumes< vector<reco::Vertex> >( iConfig.getParameter<edm::InputTag>("pvTag") );
  puTag_ = consumes< vector<PileupSummaryInfo> >( iConfig.getParameter<edm::InputTag>("puTag") );
  pfTag_ = consumes< vector<reco::PFCandidate> >( iConfig.getParameter<edm::InputTag>("pfTag") );
  rhoTag_ = consumes<double>( iConfig.getParameter<edm::InputTag>("rhoTag") );
  rhoCentral0Tag_ = consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"));
  rhoCentralChargedTag_ = consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralChargedPileUp"));
}

// ------------ method called once each job just before starting event loop  ------------
void OffsetAnalysis::beginJob()
{
  if (!isMC_)
    parsePileUpJSON2();

  ids[all] = "all"; ids[h] = "h"; ids[e] = "e"; ids[mu] = "mu"; ids[gamma] = "gamma";
  ids[h0] = "h0"; ids[h_HF] = "h_HF"; ids[egamma_HF] = "egamma_HF"; ids[chs] = "chs";

  //counter = -1;
  rand = new TRandom3;
  root_file = new TFile(RootFileName_,"RECREATE");
  TString hname;

  hname = "pv_all_z";
  m_Histos1D[hname] = new TH1F(hname,hname,100,-50,50);
  hname = "pv_z";
  m_Histos1D[hname] = new TH1F(hname,hname,100,-50,50);
  hname = "nPV_all";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,100);
  hname = "nPV0";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,100);
  hname = "nPV";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,100);
  hname = "nPU";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,50);
  hname = "rho";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,100); 

  hname = "2h_nPV_nPU";
  m_Histos2D[hname] = new TH2F(hname,hname,100,0,maxNPU_,maxNPV_,0,maxNPV_);
  hname = "p_nPV_nPU";
  m_Profiles[hname] = new TProfile(hname,hname,100,0,maxNPU_);
  hname = "2h_nPVall_nPU";
  m_Histos2D[hname] = new TH2F(hname,hname,100,0,maxNPU_,maxNPV_,0,maxNPV_);
  hname = "p_nPVall_nPU";
  m_Profiles[hname] = new TProfile(hname,hname,100,0,maxNPU_);
  hname = "2h_nPV0_nPU";
  m_Histos2D[hname] = new TH2F(hname,hname,100,0,maxNPU_,maxNPV_,0,maxNPV_);
  hname = "p_nPV0_nPU";
  m_Profiles[hname] = new TProfile(hname,hname,100,0,maxNPU_);

  hname = "2h_rho_nPU";
  m_Histos2D[hname] = new TH2F(hname,hname,100,0,maxNPU_,100,0,maxRho_);
  hname = "p_rho_nPU";
  m_Profiles[hname] = new TProfile(hname,hname,100,0,maxNPU_);

  hname = "2h_rho_nPV";
  m_Histos2D[hname] = new TH2F(hname,hname,maxNPV_,0,maxNPV_,100,0,maxRho_);
  hname = "p_rho_nPV";
  m_Profiles[hname] = new TProfile(hname,hname,maxNPV_,0,maxNPV_);

  hname = "rhoCentral0";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,100);
  hname = "2h_rhoCentral0_nPU";
  m_Histos2D[hname] = new TH2F(hname,hname,100,0,maxNPU_,100,0,maxRho_);
  hname = "p_rhoCentral0_nPU";
  m_Profiles[hname] = new TProfile(hname,hname,100,0,maxNPU_);
  hname = "2h_rhoCentral0_nPV";
  m_Histos2D[hname] = new TH2F(hname,hname,maxNPV_,0,maxNPV_,100,0,maxRho_);
  hname = "p_rhoCentral0_nPV";
  m_Profiles[hname] = new TProfile(hname,hname,maxNPV_,0,maxNPV_);

  hname = "rhoCentralCharged";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,100);
  hname = "2h_rhoCentralCharged_nPU";
  m_Histos2D[hname] = new TH2F(hname,hname,100,0,maxNPU_,100,0,maxRho_);
  hname = "p_rhoCentralCharged_nPU";
  m_Profiles[hname] = new TProfile(hname,hname,100,0,maxNPU_);
  hname = "2h_rhoCentralCharged_nPV";
  m_Histos2D[hname] = new TH2F(hname,hname,maxNPV_,0,maxNPV_,100,0,maxRho_);
  hname = "p_rhoCentralCharged_nPV";
  m_Profiles[hname] = new TProfile(hname,hname,maxNPV_,0,maxNPV_);

  for (int i_id=0; i_id<numFlavors; i_id++){
    for (int i_nPU=0; i_nPU < maxNPU_; i_nPU++){

      hname = Form("2h_cone_eta_totalpT_nPU%i_PF", i_nPU) + ids[i_id];
      m_Histos2D[hname] = new TH2F(hname,hname,ETA_BINS,etabins,1000,0,100);
      hname = Form("p_cone_eta_totalpT_nPU%i_PF", i_nPU) + ids[i_id];
      m_Profiles[hname] = new TProfile(hname, hname, ETA_BINS, etabins);

      hname = Form("2h_cone_eta_pT_nPU%i_PF", i_nPU) + ids[i_id];
      m_Histos2D[hname] = new TH2F(hname, hname,ETA_BINS,etabins,100,0,10);
      hname = Form("p_cone_eta_pT_nPU%i_PF", i_nPU) + ids[i_id];
      m_Profiles[hname] = new TProfile(hname, hname, ETA_BINS, etabins);
    }
    for (int i_nPV=0; i_nPV < maxNPV_; i_nPV++){

      hname = Form("2h_cone_eta_totalpT_nPV%i_PF", i_nPV) + ids[i_id];
      m_Histos2D[hname] = new TH2F(hname,hname,ETA_BINS,etabins,1000,0,100);
      hname = Form("p_cone_eta_totalpT_nPV%i_PF", i_nPV) + ids[i_id];
      m_Profiles[hname] = new TProfile(hname, hname, ETA_BINS, etabins);

      hname = Form("2h_cone_eta_pT_nPV%i_PF", i_nPV) + ids[i_id];
      m_Histos2D[hname] = new TH2F(hname, hname,ETA_BINS,etabins,100,0,10);
      hname = Form("p_cone_eta_pT_nPV%i_PF", i_nPV) + ids[i_id];
      m_Profiles[hname] = new TProfile(hname, hname, ETA_BINS, etabins);

      hname = Form("2p_cone_eta_nPU_totalpT_nPV%i_PF", i_nPV) + ids[i_id];
      m_Profiles2D[hname] = new TProfile2D(hname, hname,ETA_BINS,etabins, 50,0,50);
    }
  }

}

// ------------ method called for each event  ------------
void OffsetAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //counter++;
  //if (counter%numSkip_ != 0) return;

  double weight = 1;

//------------ Pileup ------------//
 
  double true_pileup;

  if (isMC_){
    edm::Handle< vector<PileupSummaryInfo> > pileups;
    iEvent.getByToken(puTag_, pileups);
    true_pileup = pileups->at(1).getTrueNumInteractions();

    if (reweight_){

      double weights[] =
      {0.76855, 1.00000, 0.36366, 0.42279, 0.15484, 0.17634, 0.29624, 0.27565, 0.10049, 0.09502, 
       0.10464, 0.01865, 0.01607, 0.01608, 0.00105, 0.00340, 0.00123, 0.00240, 0.00137, 0.03345, 
       0.01614, 0.02763, 0.02549, 0.03761, 0.01366, 0.00620, 0.01719, 0.03349, 0.09280, 0.04156, 
       0.00000, 0.00030, 0.00049, 0.00050, 0.03905, 0.06454, 0.13404, 0.11517, 0.09962, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000 };

      int index = m_Histos1D["nPU"]->FindBin(true_pileup);
      weight = weights[index-1];
    }
  }
  else
    true_pileup = getAvgPU( int(iEvent.id().run()), int(iEvent.getLuminosityBlock().luminosityBlock()) );

//------------ Primary Vertices ------------//

  edm::Handle< vector<reco::Vertex> > primaryVertices;
  iEvent.getByToken(pvTag_, primaryVertices);

  int nPV_all = primaryVertices->size();

  FillHist1D("nPV_all", nPV_all, weight);

  int nPV = 0, nPV0 = 0;
  vector<reco::Vertex>::const_iterator i_pv, endpv = primaryVertices->end();
  for (i_pv = primaryVertices->begin();  i_pv != endpv;  ++i_pv) {
	
    double z = i_pv->z();
    FillHist1D("pv_all_z", z, weight);

    if ( !i_pv->isFake() && z <= 24 && i_pv->position().rho() <= 2){
      nPV0++; //no cut on ndof

      if (i_pv->ndof() > 4){
        FillHist1D("pv_z", z, weight);
        nPV++;
      }
    }
  }
  FillHist1D("nPV", nPV, weight);
  FillHist1D("nPV0", nPV0, weight);

//------------ Rho ------------//

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoTag_, rhoHandle);
  double rho = *rhoHandle;

  FillHist1D("rho", rho, weight);

  edm::Handle<double> rhoCentral0Handle;
  iEvent.getByToken(rhoCentral0Tag_, rhoCentral0Handle);
  double rhoCentral0 = *rhoCentral0Handle;

  FillHist1D("rhoCentral0", rhoCentral0, weight);

  edm::Handle<double> rhoCentralChargedHandle;
  iEvent.getByToken(rhoCentralChargedTag_, rhoCentralChargedHandle);
  double rhoCentralCharged = *rhoCentralChargedHandle;

  FillHist1D("rhoCentralCharged", rhoCentralCharged, weight);

//------------ Histos ------------//

  int intPU = true_pileup + 0.5;
  FillHist1D("nPU", true_pileup, weight);

  FillHist2D("2h_nPV_nPU", true_pileup, nPV, weight);
  FillHist2D("2h_nPVall_nPU", true_pileup, nPV_all, weight);
  FillHist2D("2h_nPV0_nPU", true_pileup, nPV0, weight);
  FillHist2D("2h_rho_nPU", true_pileup, rho, weight);
  FillHist2D("2h_rho_nPV", nPV, rho, weight);

  FillProfile("p_nPV_nPU", true_pileup, nPV, weight);
  FillProfile("p_nPVall_nPU", true_pileup, nPV_all, weight);
  FillProfile("p_nPV0_nPU", true_pileup, nPV0, weight);
  FillProfile("p_rho_nPU", true_pileup, rho, weight);
  FillProfile("p_rho_nPV", nPV, rho, weight);

  FillHist2D("2h_rhoCentral0_nPU", true_pileup, rhoCentral0, weight);
  FillHist2D("2h_rhoCentral0_nPV", nPV, rhoCentral0, weight);
  FillProfile("p_rhoCentral0_nPU", true_pileup, rhoCentral0, weight);
  FillProfile("p_rhoCentral0_nPV", nPV, rhoCentral0, weight);

  FillHist2D("2h_rhoCentralCharged_nPU", true_pileup, rhoCentralCharged, weight);
  FillHist2D("2h_rhoCentralCharged_nPV", nPV, rhoCentralCharged, weight);
  FillProfile("p_rhoCentralCharged_nPU", true_pileup, rhoCentralCharged, weight);
  FillProfile("p_rhoCentralCharged_nPV", nPV, rhoCentralCharged, weight);

//------------ PF Candidates ------------//

  edm::Handle< vector<reco::PFCandidate> > pfCandidates;
  iEvent.getByToken(pfTag_, pfCandidates);

  double cone_phi = (rand->Uniform(-1.,1.))*M_PI;
  TString hname;

  for (int i=0; i<3; i++){

    cone_phi += 2.*M_PI*double(i)/3.;

    double px[numFlavors][ETA_BINS] = {};
    double py[numFlavors][ETA_BINS] = {};

    vector<reco::PFCandidate>::const_iterator i_pf, endpf = pfCandidates->end();
    for (i_pf = pfCandidates->begin();  i_pf != endpf;  ++i_pf) {

      Flavor pf_id = getFlavor( i_pf->particleId() );
      if (pf_id == X) continue;

      double pf_px  = i_pf->px();
      double pf_py  = i_pf->py();
      //double pf_pt = i_pf->pt();
      double pf_eta = i_pf->eta();
      double pf_phi = i_pf->phi();

      bool attached = false;  //if false, include in charged hadron subtraction
      reco::TrackRef pftrack(i_pf->trackRef());

      if (!pftrack.isNull() ) {
      
        vector<reco::Vertex>::const_iterator i_pv, endpv = primaryVertices->end();
        for (i_pv = primaryVertices->begin(); i_pv != endpv && !attached;  ++i_pv) {
        
          if ( !i_pv->isFake() && i_pv->ndof() >= 4 && fabs(i_pv->z()) < 24 ) { //Why >= 4 and not > 4?

            reco::Vertex::trackRef_iterator i_vtxTrk, endvtxTrk = i_pv->tracks_end();
            for(i_vtxTrk = i_pv->tracks_begin(); i_vtxTrk != endvtxTrk && !attached; ++i_vtxTrk) {
              
              reco::TrackRef vtxTrk(i_vtxTrk->castTo<reco::TrackRef>());
              if (vtxTrk == pftrack)
                attached = true;
            } 
          }
        }
      }

      // loop over cones centered at different eta values and add energy from the pfCandidate if it is within coneDR_ from the cone axis
      for (int i_eta=0; i_eta<ETA_BINS; i_eta++) {
        double cone_eta = 0.5*(etabins[i_eta]+etabins[i_eta+1]);
        double dPhi = deltaPhi(pf_phi, cone_phi);
        double dEta = pf_eta - cone_eta ;
        double dR   = sqrt(dPhi*dPhi + dEta*dEta);

        if (dR < coneDR_) { 

          px[all][i_eta] += pf_px;
          py[all][i_eta] += pf_py;

          px[pf_id][i_eta] += pf_px;
          py[pf_id][i_eta] += pf_py;

          if (!attached){
            px[chs][i_eta] += pf_px;
            py[chs][i_eta] += pf_py;
          }
/*
          hname = Form("2h_cone_eta_pT_nPU%i_PFall", intPU);
          FillHist2D(hname, cone_eta, pf_pt, weight);
          hname = Form("p_cone_eta_pT_nPU%i_PFall", intPU);
          FillProfile(hname, cone_eta, pf_pt, weight);

          hname = Form("2h_cone_eta_pT_nPV%i_PFall", nPV);
          FillHist2D(hname, cone_eta, pf_pt, weight);
          hname = Form("p_cone_eta_pT_nPV%i_PFall", nPV);
          FillProfile(hname, cone_eta, pf_pt, weight);

          hname = Form("2h_cone_eta_pT_nPU%i_PF", intPU) + ids[pf_id];
          FillHist2D(hname, cone_eta, pf_pt, weight);
          hname = Form("p_cone_eta_pT_nPU%i_PF", intPU) + ids[pf_id];
          FillProfile(hname, cone_eta, pf_pt, weight);

          hname = Form("2h_cone_eta_pT_nPV%i_PF", nPV) + ids[pf_id];
          FillHist2D(hname, cone_eta, pf_pt, weight);
          hname = Form("p_cone_eta_pT_nPV%i_PF", nPV) + ids[pf_id];
          FillProfile(hname, cone_eta, pf_pt, weight);

          if (!attached){
            hname = Form("2h_cone_eta_pT_nPU%i_PFchs", intPU);
            FillHist2D(hname, cone_eta, pf_pt, weight);
            hname = Form("p_cone_eta_pT_nPU%i_PFchs", intPU);
            FillProfile(hname, cone_eta, pf_pt, weight);

            hname = Form("2h_cone_eta_pT_nPV%i_PFchs", nPV);
            FillHist2D(hname, cone_eta, pf_pt, weight);
            hname = Form("p_cone_eta_pT_nPV%i_PFchs", nPV);
            FillProfile(hname, cone_eta, pf_pt, weight);
          }
*/
        } 
      } //end loop over cone        
    } //end loop over pf candidates

    // fill in histograms for total pT
    for (int i_id=0; i_id<numFlavors; i_id++){
      for (int i_eta=0; i_eta<ETA_BINS; i_eta++) {
        double cone_eta = 0.5*(etabins[i_eta]+etabins[i_eta+1]);

        double total_pt = sqrt( px[i_id][i_eta]*px[i_id][i_eta] + py[i_id][i_eta]*py[i_id][i_eta] );

        hname = Form("2h_cone_eta_totalpT_nPU%i_PF", intPU) + ids[i_id];
        FillHist2D( hname, cone_eta, total_pt, weight );
        hname = Form("p_cone_eta_totalpT_nPU%i_PF", intPU) + ids[i_id];
        FillProfile( hname, cone_eta, total_pt, weight );

        hname = Form("2h_cone_eta_totalpT_nPV%i_PF", nPV) + ids[i_id];
        FillHist2D( hname, cone_eta, total_pt, weight );
        hname = Form("p_cone_eta_totalpT_nPV%i_PF", nPV) + ids[i_id];
        FillProfile( hname, cone_eta, total_pt, weight );

        hname = Form("2p_cone_eta_nPU_totalpT_nPV%i_PF", nPV) + ids[i_id];
        FillProfile2D( hname, cone_eta, intPU, total_pt, weight );
      }
    }

  } //end loop over three cone_phi's

}

// ------------ method called once each job just after ending the event loop  ------------
void OffsetAnalysis::endJob() 
{
  if (root_file !=0) {

    TString hname, pf_id;
    map<TString, TDirectory*> directories;

    //nPU
    hname = "2h_cone_eta_totalpT_nPU";
    directories[hname] = root_file->mkdir(hname);
    hname = "p_cone_eta_totalpT_nPU";
    directories[hname] = root_file->mkdir(hname);
    hname = "2h_cone_eta_pT_nPU";
    directories[hname] = root_file->mkdir(hname);
    hname = "p_cone_eta_pT_nPU";
    directories[hname] = root_file->mkdir(hname);

    //nPV
    hname = "2h_cone_eta_totalpT_nPV";
    directories[hname] = root_file->mkdir(hname);
    hname = "p_cone_eta_totalpT_nPV";
    directories[hname] = root_file->mkdir(hname);
    hname = "2h_cone_eta_pT_nPV";
    directories[hname] = root_file->mkdir(hname);
    hname = "p_cone_eta_pT_nPV";
    directories[hname] = root_file->mkdir(hname);
    hname = "2p_cone_eta_nPU_totalpT_nPV";
    directories[hname] = root_file->mkdir(hname);

    for(map<TString, TDirectory*>::iterator i_dir = directories.begin(); i_dir != directories.end(); ++i_dir){
      for (int i_id=0; i_id<numFlavors; i_id++)
        i_dir->second->mkdir( ids[i_id] );
    }

    root_file->cd();
    for (map<TString, TH1*>::iterator hid = m_Histos1D.begin(); hid != m_Histos1D.end(); hid++)
      hid->second->Write();

    for (map<TString, TH2*>::iterator hid = m_Histos2D.begin(); hid != m_Histos2D.end(); hid++){
      root_file->cd();
      hname = hid->first;
      //NPU
      if ( hname.Contains("2h_cone_eta_pT_nPU") ){
        pf_id = hname( hname.Index("PF")+2, hname.Length() );              //make pf_id substring from hname
        root_file->cd(RootFileName_ + ":/2h_cone_eta_pT_nPU/" + pf_id);
      } else if ( hname.Contains("2h_cone_eta_totalpT_nPU") ){
        pf_id = hname( hname.Index("PF")+2, hname.Length() );
        root_file->cd(RootFileName_ + ":/2h_cone_eta_totalpT_nPU/" + pf_id);
      } //NPV
      else if ( hname.Contains("2h_cone_eta_pT_nPV") ){
        pf_id = hname( hname.Index("PF")+2, hname.Length() );
        root_file->cd(RootFileName_ + ":/2h_cone_eta_pT_nPV/" + pf_id);
      } else if ( hname.Contains("2h_cone_eta_totalpT_nPV") ){
        pf_id = hname( hname.Index("PF")+2, hname.Length() );
        root_file->cd(RootFileName_ + ":/2h_cone_eta_totalpT_nPV/" + pf_id);
      }

      hid->second->Write();
    }

    for (map<TString, TProfile*>::iterator hid = m_Profiles.begin(); hid != m_Profiles.end(); hid++){
      root_file->cd();
      hname = hid->first;
      //NPU
      if ( hname.Contains("p_cone_eta_pT_nPU") ){
        pf_id = hname( hname.Index("PF")+2, hname.Length() );
        root_file->cd(RootFileName_ + ":/p_cone_eta_pT_nPU/" + pf_id);
      } else if ( hname.Contains("p_cone_eta_totalpT_nPU") ){
        pf_id = hname( hname.Index("PF")+2, hname.Length() );
        root_file->cd(RootFileName_ + ":/p_cone_eta_totalpT_nPU/" + pf_id);
      } //NPV
      else if ( hname.Contains("p_cone_eta_pT_nPV") ){
        pf_id = hname( hname.Index("PF")+2, hname.Length() );
        root_file->cd(RootFileName_ + ":/p_cone_eta_pT_nPV/" + pf_id);
      } else if ( hname.Contains("p_cone_eta_totalpT_nPV") ){
        pf_id = hname( hname.Index("PF")+2, hname.Length() );
        root_file->cd(RootFileName_ + ":/p_cone_eta_totalpT_nPV/" + pf_id);
      }

      hid->second->Write();
    }

    for (map<TString, TProfile2D*>::iterator hid = m_Profiles2D.begin(); hid != m_Profiles2D.end(); hid++){
      root_file->cd();
      hname = hid->first;
      if ( hname.Contains("2p_cone_eta_nPU_totalpT_nPV") ){
        pf_id = hname( hname.Index("PF")+2, hname.Length() );
        root_file->cd(RootFileName_ + ":/2p_cone_eta_nPU_totalpT_nPV/" + pf_id);
      }

      hid->second->Write();
    }

    root_file->Write();
    delete root_file;
    root_file = 0;
  }
}

void OffsetAnalysis::FillHist1D(const TString& histName, const Double_t& value, Double_t weight) 
{
  map<TString, TH1*>::iterator hid=m_Histos1D.find(histName);
  if (hid==m_Histos1D.end())
    cout << "%FillHist1D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value, weight);
}

void OffsetAnalysis::FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, Double_t weight) 
{
  map<TString, TH2*>::iterator hid=m_Histos2D.find(histName);
  if (hid==m_Histos2D.end())
    cout << "%FillHist2D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, weight);
}

void OffsetAnalysis::FillProfile(const TString& histName, const Double_t& value1, const Double_t& value2, Double_t weight)
{
  map<TString, TProfile*>::iterator hid=m_Profiles.find(histName);
  if (hid==m_Profiles.end())
    cout << "%FillProfile -- Could not find profile with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, weight);
}

void OffsetAnalysis::FillProfile2D(const TString& histName, const Double_t& value1, const Double_t& value2, const Double_t& value3, Double_t weight)
{
  map<TString, TProfile2D*>::iterator hid=m_Profiles2D.find(histName);
  if (hid==m_Profiles2D.end())
    cout << "%FillProfile2D -- Could not find profile with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, value3, weight);
}

OffsetAnalysis::Flavor OffsetAnalysis::getFlavor(reco::PFCandidate::ParticleType id)
{
    if (id == reco::PFCandidate::h)
        return h;
    else if (id == reco::PFCandidate::e)
        return e;
    else if (id == reco::PFCandidate::mu)
        return mu;
    else if (id == reco::PFCandidate::gamma)
        return gamma;
    else if (id == reco::PFCandidate::h0)
        return h0;
    else if (id == reco::PFCandidate::h_HF)
        return h_HF;
    else if (id == reco::PFCandidate::egamma_HF)
        return egamma_HF;
    else
        return X;
}

double OffsetAnalysis::deltaPhi(double phi1, double phi2)
{
  double dPhi = fabs(phi1 - phi2);
  if (dPhi > M_PI) dPhi = 2.0*M_PI - dPhi;
  return dPhi;
}

//define this as a plug-in
DEFINE_FWK_MODULE(OffsetAnalysis);
