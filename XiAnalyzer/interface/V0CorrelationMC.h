#ifndef XIANALYZER__V0_CORRELATION_MC_H
#define XIANALYZER__V0_CORRELATION_MC_H

// system include files
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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

#define PI 3.1416

#define effkss(x) ((1.09925e+11*TMath::Power(x,1.50511e+01)*TMath::Exp(6.58074*x*x-2.94487e+01*x+9.72194))/(6.71504e+11*TMath::Power(x,3.19081)*TMath::Exp(1.97717*x*x-9.90447*x-4.65781)))

#define effksb(x) ((6.47559e+02*TMath::Power(x,2.95369e-01)*TMath::Exp(9.18237e-02*x*x-1.89678*x+7.63891))/(1.20601e+01*TMath::Power(x,-1.86165)*TMath::Exp(4.17033e-02*x*x-9.39659e-01*x+1.30386e+01)))

#define efflas(x) ((4.76443e+06*TMath::Power(x,1.62889e+01)*TMath::Exp(2.74004*x*x-1.97581e+01*x+1.16309e+01))/(5.30771e+08*TMath::Power(x,3.59273)*TMath::Exp(1.50703*x*x-8.45701*x+9.43797e-01)))

#define efflab(x) ((3.86297e-01*TMath::Power(x,1.91207)*TMath::Exp(8.37588e-02*x*x-2.15583*x+1.33689e+01))/(7.01220*TMath::Power(x,-4.80662e-01)*TMath::Exp(7.33837e-02*x*x-1.53854*x+1.35978e+01)))

using namespace std;

class V0CorrelationMC : public edm::EDAnalyzer {
public:
    explicit V0CorrelationMC(const edm::ParameterSet&);
    ~V0CorrelationMC();


private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------

    edm::EDGetTokenT<reco::GenParticleCollection> _gnCollection;
    edm::EDGetTokenT<reco::TrackCollection> _trkSrc;
    edm::EDGetTokenT<reco::VertexCollection> _vertexCollName;

    TH1D* hMult;

    TH1D* hRap_ks[18];
    TH1D* hRap_la[18];
    TH1D* hRap_xi[18];

    TH1D* hKET_ks[18];
    TH1D* hKET_la[18];
    TH1D* hKET_xi[18];

    TH1D* hPt_ks[18];
    TH1D* hPt_la[18];
    TH1D* hPt_xi[18];

    TH1D* hMult_ks[18];
    TH2D* hSignal_ks[18];
    TH2D* hBackground_ks[18];

    TH1D* hMult_ass;

    TH1D* h2Daughter;

    TH1D* hMult_la[18];
    TH2D* hSignal_la[18];
    TH2D* hBackground_la[18];

    TH1D* hMult_xi[18];
    TH2D* hSignal_xi[18];
    TH2D* hBackground_xi[18];

    vector<TVector3> *pVect_trg_ks[18];
    vector< vector<TVector3> > *pVectVect_trg_ks[18];
    vector<TVector3> *pVect_trg_la[18];
    vector< vector<TVector3> > *pVectVect_trg_la[18];
    vector<TVector3> *pVect_trg_xi[18];
    vector< vector<TVector3> > *pVectVect_trg_xi[18];

    vector<TVector3> *pVect_dau_ks[18];
    vector<TVector3> *pVect_dau_la[18];

    vector<TVector3> *pVect_ass;
    vector< vector<TVector3> > *pVectVect_ass;
    vector<double> *zvtxVect;
    vector<double> ptcut_ks_;
    vector<double> ptcut_la_;
    vector<double> ptcut_xi_;

    double etaMin_trg_;
    double etaMax_trg_;
    double etaMin_ass_;
    double etaMax_ass_;
    double ptMin_ass_;
    double ptMax_ass_;
    double rapMin_;
    double rapMax_;
    double multMax_;
    double multMin_;
    int bkgFactor_;
    bool doRap_;
    bool doKs_;
    bool doLa_;
    bool doXi_;
};
#endif
