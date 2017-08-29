#ifndef XIANALYZER__HAD_CORR_H
#define XIANALYZER__HAD_CORR_H
// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

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
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

//
// class declaration
//

#define PI 3.1416
using namespace std;

class HadronCorrelationGen : public edm::EDAnalyzer {
public:
    explicit HadronCorrelationGen(const edm::ParameterSet&);
    ~HadronCorrelationGen();
    
    
private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    
    // ----------member data ---------------------------

    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> _vertexCollName;
    edm::EDGetTokenT<reco::TrackCollection> _trkSrc;
    edm::EDGetTokenT<reco::GenParticleCollection> _gnCollection;

    TH1D* hMult_selected;
    TH1D* HadPerEvt;
    
    TH1D* hPt[18];
    
    TH1D* hMult[18];
    TH2D* hSignal[18];
    TH2D* hBackground[18];
    TH1D* hMult_good[18];
    TH1D* hMult_assoc[18];
    
    TH2F* effhisto;
    TH2D* SignalHad;
    TH2D* BackgroundHad;
    
    vector<TVector3> *pVect_trg[18];
    vector< vector<TVector3> > *pVectVect_trg[18];
    
    vector<TVector3> *pVect_ass;
    vector< vector<TVector3> > *pVectVect_ass;
    vector< vector<TVector3> > *pVect2_ass;
    vector<double> *zvtxVect;
    
    double etaMin_trg_;
    double etaMax_trg_;
    double etaMin_ass_;
    double etaMax_ass_;
    double ptMin_ass_;
    double ptMax_ass_;
    int bkgFactor_;
    double multMax_;
    double multMin_;
    double rapMax_;
    double rapMin_;
    edm::Service<TFileService> fs;
};

#endif
