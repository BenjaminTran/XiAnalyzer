#ifndef XIANALYZER__V0CAS_CORRELATION_H
#define XIANALYZER__V0CAS_CORRELATION_H

#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


#define PI 3.1416

using namespace std;

class V0CasCorrelation : public edm::EDAnalyzer {
public:
    explicit V0CasCorrelation(const edm::ParameterSet&);
    ~V0CasCorrelation();


private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------

    TH1D* hMult;
    TH1D* hMult_accept;

    TH2F* effhisto;
    TH2D* effhisto_ks;
    TH2D* effhisto_la;

    TH1D* hMass_ks[18];
    TH1D* hMass_la[18];

    TH1D* hKET_ks[18];
    TH1D* hKET_la[18];

    TH1D* hPt_ks[18];
    TH1D* hPt_la[18];

    TH1D* hEta_ks[18];
    TH1D* hEta_la[18];
    TH1D* hRap_ks[18];
    TH1D* hRap_la[18];
    TH1D* hRap_ks_Lorentz[18];

    TH1D* hMult_ks[18];
    TH2D* hSignal_ks[18];
    TH2D* hBackground_ks[18];

    TH1D* hKET_ks_bkg[18];
    TH1D* hKET_la_bkg[18];

    TH1D* hPt_ks_bkg[18];
    TH1D* hPt_la_bkg[18];

    TH1D* hEta_ks_bkg[18];
    TH1D* hEta_la_bkg[18];
    TH1D* hRap_ks_bkg[18];
    TH1D* hRap_la_bkg[18];

    TH1D* hMult_ks_bkg[18];
    TH2D* hSignal_ks_bkg[18];
    TH2D* hBackground_ks_bkg[18];

    TH1D* hMult_ass;

    TH1D* hMult_la[18];
    TH2D* hSignal_la[18];
    TH2D* hBackground_la[18];

    vector<TLorentzVector> *pVect_trg_ks[18];
    vector< vector<TLorentzVector> > *pVectVect_trg_ks[18];
    vector<TLorentzVector> *pVect_trg_la[18];
    vector< vector<TLorentzVector> > *pVectVect_trg_la[18];

    vector<TVector3> *pVect_dau_ks[18];
    vector< vector<TVector3> > *pVectVect_dau_ks[18];
    vector<TVector3> *pVect_dau_la[18];
    vector< vector<TVector3> > *pVectVect_dau_la[18];

    TH1D* hMult_la_bkg[18];
    TH2D* hSignal_la_bkg[18];
    TH2D* hBackground_la_bkg[18];

    vector<TLorentzVector> *pVect_trg_ks_bkg[18];
    vector< vector<TLorentzVector> > *pVectVect_trg_ks_bkg[18];
    vector<TLorentzVector> *pVect_trg_la_bkg[18];
    vector< vector<TLorentzVector> > *pVectVect_trg_la_bkg[18];

    vector<TVector3> *pVect_dau_ks_bkg[18];
    vector< vector<TVector3> > *pVectVect_dau_ks_bkg[18];
    vector<TVector3> *pVect_dau_la_bkg[18];
    vector< vector<TVector3> > *pVectVect_dau_la_bkg[18];

    vector<TVector3> *pVect_ass;
    vector< vector<TVector3> > *pVectVect_ass;
    vector<double> *zvtxVect;

    vector<double> ptcut_ks_;
    vector<double> ptcut_la_;
    vector<double> sigma_ks_;
    vector<double> mean_ks_;
    vector<double> sigma_la_;
    vector<double> mean_la_;

    int ptbin_n_;
    int cent_bin_low_;
    int cent_bin_high_;
    double etaMin_trg_;
    double etaMax_trg_;
    double etaMin_ass_;
    double etaMax_ass_;
    double ptMin_ass_;
    double ptMax_ass_;
    double multMax_;
    double multMin_;
    int bkgnum_;
    double peakFactor_;
    double sideFactor_;
    double mis_ks_range_;
    double mis_la_range_;
    double mis_ph_range_;
    bool rejectDaughter_;
    bool useCentrality_;
    bool doRap_;
    bool doGenRef_;
    bool doThrowAway_;
    bool doKs_;
    bool doLa_;
    bool doXi_;
    bool doOm_;

    edm::EDGetTokenT<reco::VertexCollection> _vertexCollName;
    edm::EDGetTokenT<reco::TrackCollection> _trkSrc;
    edm::EDGetTokenT<CaloTowerCollection> _towerSrc;

    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> _ksCollection;
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> _laCollection;
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> _xiCollection;
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> _omCollection;
    edm::EDGetTokenT<reco::GenParticleCollection> _gnCollection;


    int ptbin_n_cas_;
    vector<TLorentzVector> *pepVect_xi_peak[18];
    vector<TLorentzVector> *pepVect_xi_side[18];
    vector<TVector3> *pepVect_dau_xi_peak[18];
    vector<TVector3> *pepVect_dau_xi_side[18];
    vector<TLorentzVector> *pepVect_om_peak[18];
    vector<TLorentzVector> *pepVect_om_side[18];
    vector<TVector3> *pepVect_dau_om_peak[18];
    vector<TVector3> *pepVect_dau_om_side[18];
    vector<double> ptcut_xi_;
    vector<double> ptcut_om_;
    vector<double> mean_xi_;
    vector<double> sigma_xi_;
    vector<double> mean_om_;
    vector<double> sigma_om_;
    vector< vector<TLorentzVector> > *PepVect2_xi_peak[18];
    vector< vector<TLorentzVector> > *PepVect2_xi_side[18];
    vector< vector<TLorentzVector> > *PepVect2_om_peak[18];
    vector< vector<TLorentzVector> > *PepVect2_om_side[18];
    vector<vector<TVector3> > *PepVect2_dau_xi_peak[18];
    vector<vector<TVector3> > *PepVect2_dau_xi_side[18];
    vector<vector<TVector3> > *PepVect2_dau_om_peak[18];
    vector<vector<TVector3> > *PepVect2_dau_om_side[18];

    TH2D* MassPtXi;
    TH2D* MassPtOm;
    TH2D* BackgroundXiPeak[18];
    TH2D* BackgroundXiSide[18];
    TH2D* SignalXiPeak[18];
    TH2D* SignalXiSide[18];

    TH1D* KET_xi[18];
    TH1D* KET_xi_bkg[18];
    TH1D* Mass_xi[18];
    TH1D* Pt_xi[18];
    TH1D* Pt_xi_bkg[18];
    TH1D* Eta_xi[18];
    TH1D* Eta_xi_bkg[18];
    TH1D* rap_xi[18];
    TH1D* rap_xi_bkg[18];
    TH1D* rap_xi_Lorentz[18];
    TH1D* mult_xi[18];
    TH1D* mult_xi_bkg[18];

    TH2D* BackgroundOmPeak[18];
    TH2D* BackgroundOmSide[18];
    TH2D* SignalOmPeak[18];
    TH2D* SignalOmSide[18];

    TH1D* KET_om[18];
    TH1D* KET_om_bkg[18];
    TH1D* Mass_om[18];
    TH1D* Pt_om[18];
    TH1D* Pt_om_bkg[18];
    TH1D* Eta_om[18];
    TH1D* Eta_om_bkg[18];
    TH1D* rap_om[18];
    TH1D* rap_om_bkg[18];
    TH1D* rap_om_Lorentz[18];
    TH1D* mult_om[18];
    TH1D* mult_om_bkg[18];
};

#endif
