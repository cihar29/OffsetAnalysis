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
#include "parsePileUpJSON2.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/Ref.h"
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include <vector>
#include <map>
#include <utility>
#include <cmath>

//root files
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom3.h>

using namespace std;

class OffsetAnalysis : public edm::EDAnalyzer {
  public:
    explicit OffsetAnalysis(const edm::ParameterSet&);

  private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    void FillHist1D(const TString& histName, const Double_t& value);
    void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2);
    void FillProfile(const TString& histName, const Double_t& value1, const Double_t& value2);
    void FillProfile2D(const TString& histName, const Double_t& value1, const Double_t& value2, const Double_t& value3);
    double deltaPhi(double phi1, double phi2);
    TString getPFCandidateId(const int& id);
	  
    TRandom3* rand;
    TFile* root_file;
    map<TString, TH1*> m_Histos1D;
    map<TString, TH2*> m_Histos2D;  
    map<TString, TProfile*> m_Profiles;
    map<TString, TProfile2D*> m_Profiles2D;

    TString RootFileName_;
    double coneDR_;
    int maxNPV_, maxNPU_, maxRho_;
    bool isMC_, reweight_;
    edm::InputTag pvTag_, puTag_, pfTag_, rhoTag_;

    //map<int, map<int, pair<int, int> > > m_avgNPV;
};

OffsetAnalysis::OffsetAnalysis(const edm::ParameterSet& iConfig)
{
  RootFileName_ = iConfig.getParameter<string> ("RootFileName");
  coneDR_ = iConfig.getParameter<double> ("coneDR");
  maxNPV_ = iConfig.getParameter<int> ("maxNPV");
  maxNPU_ = iConfig.getParameter<int> ("maxNPU");
  maxRho_ = iConfig.getParameter<int> ("maxRho");
  isMC_ = iConfig.getParameter<bool> ("isMC");
  reweight_ = iConfig.getParameter<bool> ("reweight");
  pvTag_ = iConfig.getParameter<edm::InputTag>("pvTag");
  puTag_ = iConfig.getParameter<edm::InputTag>("puTag");
  pfTag_ = iConfig.getParameter<edm::InputTag>("pfTag");
  rhoTag_ = iConfig.getParameter<edm::InputTag>("rhoTag");
}

// ------------ method called once each job just before starting event loop  ------------
void OffsetAnalysis::beginJob()
{
  if (!isMC_)
    parsePileUpJSON2();

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

  //PF Histos//
  TString ids[] = {"all", "h", "e", "mu", "gamma", "h0", "h_HF", "egamma_HF", "chs"};

  for (int i_id=0; i_id < 9; i_id++){
    for (int i_nPU=0; i_nPU < maxNPU_; i_nPU++){

      hname = Form("2h_cone_eta_totalpT_nPU%i_PF", i_nPU) + ids[i_id];
      m_Histos2D[hname] = new TH2F(hname,hname,100,-5.,5.,1000,-0.5,99.5);
      hname = Form("p_cone_eta_totalpT_nPU%i_PF", i_nPU) + ids[i_id];
      m_Profiles[hname] = new TProfile(hname, hname, 100, -5., 5.);

      hname = Form("2h_cone_eta_pT_nPU%i_PF", i_nPU) + ids[i_id];
      m_Histos2D[hname] = new TH2F(hname, hname, 100,-5.,5.,100, 0.,10.);
      hname = Form("p_cone_eta_pT_nPU%i_PF", i_nPU) + ids[i_id];
      m_Profiles[hname] = new TProfile(hname, hname, 100, -5., 5.);
    }
    for (int i_nPV=0; i_nPV < maxNPV_; i_nPV++){

      hname = Form("2h_cone_eta_totalpT_nPV%i_PF", i_nPV) + ids[i_id];
      m_Histos2D[hname] = new TH2F(hname,hname,100,-5.,5.,1000,-0.5,99.5);
      hname = Form("p_cone_eta_totalpT_nPV%i_PF", i_nPV) + ids[i_id];
      m_Profiles[hname] = new TProfile(hname, hname, 100, -5., 5.);

      hname = Form("2h_cone_eta_pT_nPV%i_PF", i_nPV) + ids[i_id];
      m_Histos2D[hname] = new TH2F(hname, hname, 100,-5.,5.,100, 0.,10.);
      hname = Form("p_cone_eta_pT_nPV%i_PF", i_nPV) + ids[i_id];
      m_Profiles[hname] = new TProfile(hname, hname, 100, -5., 5.);

      hname = Form("2p_cone_eta_nPU_totalpT_nPV%i_PF", i_nPV) + ids[i_id];
      m_Profiles2D[hname] = new TProfile2D(hname, hname, 100, -5., 5., 50, -0.5, 49.5);
    }
  }

}

// ------------ method called for each event  ------------
void OffsetAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//------------ Pileup ------------//
 
  double true_pileup;

  if (isMC_){
    edm::Handle< edm::View<PileupSummaryInfo> > pileups;
    iEvent.getByLabel(puTag_, pileups);
    true_pileup = pileups->at(1).getTrueNumInteractions();

    if (reweight_){

      float weights[] =
      {0.00000, 304.47073, 109.69337, 42.31831, 98.72736, 11.58533, 3.49539, 4.43090, 5.55305, 0.69875, 
       0.35308, 0.37753, 0.03418, 0.00052, 0.00716, 0.00987, 0.00655, 0.00326, 0.00702, 0.48605, 
       0.38253, 0.30134, 0.71854, 1.00000, 0.80642, 0.53452, 0.32552, 0.23700, 0.02542, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000 };

      int index = m_Histos1D["nPU"]->FindBin(true_pileup);
      if ( rand->Rndm() >= weights[index-1] )
        return;
    }
  }
  else
    true_pileup = getAvgPU( int(iEvent.id().run()), int(iEvent.getLuminosityBlock().luminosityBlock()) );

//------------ Primary Vertices ------------//

  edm::Handle< edm::View<reco::Vertex> > primaryVertices;
  iEvent.getByLabel(pvTag_, primaryVertices);

  int nPV_all = primaryVertices->size();
  if (!isMC_ && nPV_all <= 1) return;

  FillHist1D("nPV_all", nPV_all);

  int nPV = 0; int nPV0 = 0;
  edm::View<reco::Vertex>::const_iterator i_pv, endpv = primaryVertices->end();
  for (i_pv = primaryVertices->begin();  i_pv != endpv;  ++i_pv) {
	
    double z = i_pv->z();
    FillHist1D("pv_all_z", z);

    if ( !i_pv->isFake() && z <= 24 && i_pv->position().rho() <= 2){
      nPV0++; //no cut on ndof

      if (i_pv->ndof() > 4){
        FillHist1D("pv_z", z);
        nPV++;
      }
    }
  }
  FillHist1D("nPV", nPV);
  FillHist1D("nPV0", nPV0);
/*
//------------ nPV by LS ------------//

  int run = int(iEvent.id().run());
  int ls = int(iEvent.getLuminosityBlock().luminosityBlock());

  if ( m_avgNPV.find(run) == m_avgNPV.end() || m_avgNPV.find(run)->second.find(ls) == m_avgNPV.find(run)->second.end() )
    m_avgNPV[run][ls] = make_pair(1, nPV);
  else{
    m_avgNPV[run][ls].first ++;
    m_avgNPV[run][ls].second += nPV;
  }
*/
//------------ Rho ------------//

  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(rhoTag_, rhoHandle);
  double rho = *rhoHandle;

  FillHist1D("rho", rho);

  edm::Handle<double> rhoCentral0Handle;
  iEvent.getByLabel("fixedGridRhoFastjetCentralNeutral", rhoCentral0Handle);
  double rhoCentral0 = *rhoCentral0Handle;

  FillHist1D("rhoCentral0", rhoCentral0);

  edm::Handle<double> rhoCentralChargedHandle;
  iEvent.getByLabel("fixedGridRhoFastjetCentralChargedPileUp", rhoCentralChargedHandle);
  double rhoCentralCharged = *rhoCentralChargedHandle;

  FillHist1D("rhoCentralCharged", rhoCentralCharged);

//------------ Histos ------------//

  int intPU = true_pileup + 0.5;
  FillHist1D("nPU", true_pileup);

  FillHist2D("2h_nPV_nPU", true_pileup, nPV);
  FillHist2D("2h_nPVall_nPU", true_pileup, nPV_all);
  FillHist2D("2h_nPV0_nPU", true_pileup, nPV0);
  FillHist2D("2h_rho_nPU", true_pileup, rho);
  FillHist2D("2h_rho_nPV", nPV, rho);

  FillProfile("p_nPV_nPU", true_pileup, nPV);
  FillProfile("p_nPVall_nPU", true_pileup, nPV_all);
  FillProfile("p_nPV0_nPU", true_pileup, nPV0);
  FillProfile("p_rho_nPU", true_pileup, rho);
  FillProfile("p_rho_nPV", nPV, rho);

  FillHist2D("2h_rhoCentral0_nPU", true_pileup, rhoCentral0);
  FillHist2D("2h_rhoCentral0_nPV", nPV, rhoCentral0);
  FillProfile("p_rhoCentral0_nPU", true_pileup, rhoCentral0);
  FillProfile("p_rhoCentral0_nPV", nPV, rhoCentral0);

  FillHist2D("2h_rhoCentralCharged_nPU", true_pileup, rhoCentralCharged);
  FillHist2D("2h_rhoCentralCharged_nPV", nPV, rhoCentralCharged);
  FillProfile("p_rhoCentralCharged_nPU", true_pileup, rhoCentralCharged);
  FillProfile("p_rhoCentralCharged_nPV", nPV, rhoCentralCharged);

//------------ PF Candidates ------------//

  edm::Handle< edm::View<reco::PFCandidate> > pfCandidates;
  iEvent.getByLabel(pfTag_, pfCandidates);

  double cone_phi = (rand->Uniform(-1.,1.))*M_PI;
  TString hname;

  for (int i=0; i<3; i++){

    cone_phi += 2.*M_PI*double(i)/3.;

    map< TString, vector< pair<double, double> > > m_pxpy_cone_PF;

    m_pxpy_cone_PF["all"].assign(100, make_pair(0., 0.) );
    m_pxpy_cone_PF["h"].assign(100, make_pair(0., 0.) );
    m_pxpy_cone_PF["e"].assign(100, make_pair(0., 0.) );
    m_pxpy_cone_PF["mu"].assign(100, make_pair(0., 0.) );
    m_pxpy_cone_PF["gamma"].assign(100, make_pair(0., 0.) );
    m_pxpy_cone_PF["h0"].assign(100, make_pair(0., 0.) );
    m_pxpy_cone_PF["h_HF"].assign(100, make_pair(0., 0.) );
    m_pxpy_cone_PF["egamma_HF"].assign(100, make_pair(0., 0.) );
    m_pxpy_cone_PF["chs"].assign(100, make_pair(0., 0.) );

    edm::View<reco::PFCandidate>::const_iterator i_pf, endpf = pfCandidates->end();
    for (i_pf = pfCandidates->begin();  i_pf != endpf;  ++i_pf) {

      TString pf_id = getPFCandidateId( i_pf->particleId() );    
      //double pf_pt  = i_pf->pt();
      double pf_px  = i_pf->px();
      double pf_py  = i_pf->py();
      double pf_eta = i_pf->eta();
      double pf_phi = i_pf->phi();

      bool attached = false;  //if false, include in charged hadron subtraction
      reco::TrackRef pftrack(i_pf->trackRef());

      if (!pftrack.isNull() ) {
      
        edm::View<reco::Vertex>::const_iterator i_pv, endpv = primaryVertices->end();
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
      for (int i=0; i<100; i++) {
        double cone_eta = -4.95 + double(i)*0.1;
        double dPhi = deltaPhi(pf_phi, cone_phi);
        double dEta = pf_eta - cone_eta ;
        double dR   = sqrt(dPhi*dPhi + dEta*dEta);

        if (dR < coneDR_) { 

          m_pxpy_cone_PF["all"][i].first += pf_px;
          m_pxpy_cone_PF["all"][i].second += pf_py;

          m_pxpy_cone_PF[pf_id][i].first += pf_px;
          m_pxpy_cone_PF[pf_id][i].second += pf_py;

          if (!attached){
            m_pxpy_cone_PF["chs"][i].first += pf_px;
            m_pxpy_cone_PF["chs"][i].second += pf_py;
          }
/*
          hname = Form("2h_cone_eta_pT_nPU%i_PFall", intPU);
          FillHist2D(hname, cone_eta, pf_pt);
          hname = Form("p_cone_eta_pT_nPU%i_PFall", intPU);
          FillProfile(hname, cone_eta, pf_pt);

          hname = Form("2h_cone_eta_pT_nPV%i_PFall", nPV);
          FillHist2D(hname, cone_eta, pf_pt);
          hname = Form("p_cone_eta_pT_nPV%i_PFall", nPV);
          FillProfile(hname, cone_eta, pf_pt);

          hname = Form("2h_cone_eta_pT_nPU%i_PF", intPU) + pf_id;
          FillHist2D(hname, cone_eta, pf_pt);
          hname = Form("p_cone_eta_pT_nPU%i_PF", intPU) + pf_id;
          FillProfile(hname, cone_eta, pf_pt);

          hname = Form("2h_cone_eta_pT_nPV%i_PF", nPV) + pf_id;
          FillHist2D(hname, cone_eta, pf_pt);
          hname = Form("p_cone_eta_pT_nPV%i_PF", nPV) + pf_id;
          FillProfile(hname, cone_eta, pf_pt);

          if (!attached){
            hname = Form("2h_cone_eta_pT_nPU%i_PFchs", intPU);
            FillHist2D(hname, cone_eta, pf_pt);
            hname = Form("p_cone_eta_pT_nPU%i_PFchs", intPU);
            FillProfile(hname, cone_eta, pf_pt);

            hname = Form("2h_cone_eta_pT_nPV%i_PFchs", nPV);
            FillHist2D(hname, cone_eta, pf_pt);
            hname = Form("p_cone_eta_pT_nPV%i_PFchs", nPV);
            FillProfile(hname, cone_eta, pf_pt);
          }
*/
        } 
      } //end loop over cone        
    } //end loop over pf candidates

    // fill in histograms for total pT
    for (map< TString, vector< pair<double, double> > >::iterator i_pf = m_pxpy_cone_PF.begin(); i_pf != m_pxpy_cone_PF.end(); i_pf++){

      TString pf_id = i_pf->first;
      for (int i=0; i<100; i++) {
        double cone_eta = -4.95 + double(i)*0.1;

        double total_px = i_pf->second[i].first;
        double total_py = i_pf->second[i].second;
        double total_pt = sqrt( total_px*total_px + total_py*total_py );

        hname = Form("2h_cone_eta_totalpT_nPU%i_PF", intPU) + pf_id;
        FillHist2D( hname, cone_eta, total_pt );
        hname = Form("p_cone_eta_totalpT_nPU%i_PF", intPU) + pf_id;
        FillProfile( hname, cone_eta, total_pt );

        hname = Form("2h_cone_eta_totalpT_nPV%i_PF", nPV) + pf_id;
        FillHist2D( hname, cone_eta, total_pt );
        hname = Form("p_cone_eta_totalpT_nPV%i_PF", nPV) + pf_id;
        FillProfile( hname, cone_eta, total_pt );

        hname = Form("2p_cone_eta_nPU_totalpT_nPV%i_PF", nPV) + pf_id;
        FillProfile2D( hname, cone_eta, intPU, total_pt );
      }
    }

  } //end loop over three cone_phi's

}

// ------------ method called once each job just after ending the event loop  ------------
void OffsetAnalysis::endJob() 
{
/*
  ofstream writeFile( "avgNPV.txt" );
  
  writeFile << "Run" << "\t\t" << "LS" << "\t" << "avg_nPV" << endl;

  for (map<int, map<int, pair<int, int> > >::iterator i_run = m_avgNPV.begin(); i_run != m_avgNPV.end(); i_run++){

    map<int, pair<int, int> > m_ls = i_run->second;

    for (map<int, pair<int, int> >::iterator i_ls = m_ls.begin(); i_ls != m_ls.end(); i_ls++){

      writeFile << i_run->first << "\t" << i_ls->first << "\t" << double(i_ls->second.second) / i_ls->second.first << endl;
    }
  }

  writeFile.close();
*/
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

    TString ids[] = {"all", "h", "e", "mu", "gamma", "h0", "h_HF", "egamma_HF", "chs"};

    for(map<TString, TDirectory*>::iterator i_dir = directories.begin(); i_dir != directories.end(); ++i_dir){
      for (int i_id=0; i_id < 9; i_id++)
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

void OffsetAnalysis::FillHist1D(const TString& histName, const Double_t& value) 
{
  map<TString, TH1*>::iterator hid=m_Histos1D.find(histName);
  if (hid==m_Histos1D.end())
    cout << "%FillHist1D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value);
}

void OffsetAnalysis::FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2) 
{
  map<TString, TH2*>::iterator hid=m_Histos2D.find(histName);
  if (hid==m_Histos2D.end())
    cout << "%FillHist2D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2);
}

void OffsetAnalysis::FillProfile(const TString& histName, const Double_t& value1, const Double_t& value2) 
{
  map<TString, TProfile*>::iterator hid=m_Profiles.find(histName);
  if (hid==m_Profiles.end())
    cout << "%FillProfile -- Could not find profile with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2);
}

void OffsetAnalysis::FillProfile2D(const TString& histName, const Double_t& value1, const Double_t& value2, const Double_t& value3) 
{
  map<TString, TProfile2D*>::iterator hid=m_Profiles2D.find(histName);
  if (hid==m_Profiles2D.end())
    cout << "%FillProfile2D -- Could not find profile with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, value3);
}

double OffsetAnalysis::deltaPhi(double phi1, double phi2)
{
  double dPhi = fabs(phi1 - phi2);
  if (dPhi > M_PI) dPhi = 2.0*M_PI - dPhi;
  return dPhi;
}

//converts enum to string
TString OffsetAnalysis::getPFCandidateId(const int& id)
{
	if (id == 0)
		return "X";         //undefined
	else if (id == 1)
		return "h";         //charged hadron
	else if (id == 2)
		return "e";         //electron
	else if (id == 3)
		return "mu";        //muon
	else if (id == 4)
		return "gamma";     //photon
	else if (id == 5)
		return "h0";        //neutral hadron
	else if (id == 6)
		return "h_HF";      //HF hadron
	else if (id == 7)
		return "egamma_HF"; //HF EM particle
	else
		return "no Id";
}

//define this as a plug-in
DEFINE_FWK_MODULE(OffsetAnalysis);
