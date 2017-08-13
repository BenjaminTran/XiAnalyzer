#ifndef XIANALYZER__V0_CORRELATION_H
#define XIANALYZER__V0_CORRELATION_H

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

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


#define PI 3.1416

using namespace std;

class V0Correlation : public edm::EDAnalyzer {
public:
    explicit V0Correlation(const edm::ParameterSet&);
    ~V0Correlation();
    
    
private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    
    // ----------member data ---------------------------
    
    TH1D* hMult;
    
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
    
    vector<TVector3> *pVect_trg_ks[18];
    vector< vector<TVector3> > *pVectVect_trg_ks[18];
    vector<TVector3> *pVect_trg_la[18];
    vector< vector<TVector3> > *pVectVect_trg_la[18];
    
    vector<TVector3> *pVect_dau_ks[18];
    vector< vector<TVector3> > *pVectVect_dau_ks[18];
    vector<TVector3> *pVect_dau_la[18];
    vector< vector<TVector3> > *pVectVect_dau_la[18];
    
    TH1D* hMult_la_bkg[18];
    TH2D* hSignal_la_bkg[18];
    TH2D* hBackground_la_bkg[18];
    
    vector<TVector3> *pVect_trg_ks_bkg[18];
    vector< vector<TVector3> > *pVectVect_trg_ks_bkg[18];
    vector<TVector3> *pVect_trg_la_bkg[18];
    vector< vector<TVector3> > *pVectVect_trg_la_bkg[18];
    
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
    
    edm::EDGetTokenT<reco::VertexCollection> _vertexCollName;
    edm::EDGetTokenT<reco::TrackCollection> _trkSrc;
    
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> _ksCollection;
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> _laCollection;
};

#endif
