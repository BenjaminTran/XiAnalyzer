// -*- C++ -*-
//
// Package:    HighMultFilter/HighMultFilter
// Class:      HighMultFilter
//
/**\class HighMultFilter HighMultFilter.cc HighMultFilter/HighMultFilter/plugins/HighMultFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhoudunming Tu
//         Created:  Sun, 15 May 2016 13:11:33 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>
#include <sstream>


#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TGraph.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//////////////////////////////////////////////
// CMSSW user include files
#include "DataFormats/Common/interface/DetSetAlgorithm.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"


#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

// Vertex significance
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

// Root include files
#include "TTree.h"
//
// Track Matching and fake rate calculations
//#include "RiceHIG/V0Analysis/interface/V0Validator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class HighMultFilter : public edm::EDFilter {
public:
    explicit HighMultFilter(const edm::ParameterSet&);
    ~HighMultFilter();
    virtual void endJob() ;
    virtual bool filter(edm::Event&, const edm::EventSetup&);
private:

    edm::EDGetTokenT<reco::GenParticleCollection> genSrc_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
    edm::EDGetTokenT<edm::View<reco::Track> > trackSrc_;

    double multMax_;
    double multMin_;
    double etaMax_;
    double etaMin_;

    bool doGenParticle_;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HighMultFilter::HighMultFilter(const edm::ParameterSet& iConfig) :
genSrc_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genSrc"))),
vertexSrc_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
trackSrc_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("trackSrc"))),
multMax_(iConfig.getParameter<double>("multMax")),
multMin_(iConfig.getParameter<double>("multMin")),
etaMax_(iConfig.getParameter<double>("etaMax")),
etaMin_(iConfig.getParameter<double>("etaMin")),
doGenParticle_(iConfig.getParameter<bool>("doGenParticle"))
{

}


HighMultFilter::~HighMultFilter()
{

}

bool
HighMultFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    bool accepted = false;

    double nMult_ass_good = 0;

    if( doGenParticle_ ){

    edm::Handle<reco::GenParticleCollection> genpars;
    iEvent.getByToken(genSrc_,genpars);

    for(unsigned it=0; it<genpars->size(); ++it){

        const reco::GenParticle & trk = (*genpars)[it];

        double eta = trk.eta();
        if(eta>etaMax_ || eta<etaMin_) continue;
        double pt  = trk.pt();
        if(pt<=0.4) continue;
        if(trk.status()!=1) continue;
        if(fabs(trk.charge())!=1) continue;

        nMult_ass_good++;
    }

    }
    else{

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vertexSrc_,vertices);
    double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();

    Handle<edm::View<reco::Track>> tracks;
    iEvent.getByToken(trackSrc_, tracks);

      for(unsigned it = 0; it < tracks->size(); it++){

       const reco::Track & trk = (*tracks)[it];

       math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

          double dzvtx = trk.dz(bestvtx);
          double dxyvtx = trk.dxy(bestvtx);
          double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
          double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);

          if(!trk.quality(reco::TrackBase::highPurity)) continue;
          if(fabs(trk.ptError())/trk.pt() > 0.1 ) continue;
          if(fabs(dzvtx/dzerror) > 3.0) continue;
          if(fabs(dxyvtx/dxyerror) > 3.0) continue;
          if(fabs(trk.eta()) < 2.4 && trk.pt() > 0.4 ){nMult_ass_good++;}// NtrkOffline

      }
    }

    if(nMult_ass_good>=multMin_ && nMult_ass_good<multMax_) accepted = true;

    return accepted;
}
void
HighMultFilter::endJob()
{
}
DEFINE_FWK_MODULE(HighMultFilter);