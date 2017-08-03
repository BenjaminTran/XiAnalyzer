#ifndef XIANALYZER__XI_OMTTREE_H
#define XIANALYZER__XI_OMTTREE_H

// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/IPTools/interface/IPTools.h"


#include <TString.h>
#include <TVector3.h>
#include <TMatrixD.h>

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//ROOT includes
#include "TTree.h"

using namespace std;
using namespace reco;
using std::vector;

struct XiData{
                //Skim Cuts
                float trkchi2cut_;
                float dauTransImpactSig_; //La_pion
                float dauLongImpactSig_; //La_pion
                float batDauTransImpactSig_; //Xi_pion
                float batDauLongImpactSig_; //Xi_pion
                float xiVtxSignificance3D_; //Xi_dls
                float vtxSignificance3D_; //La_dls
                //Hong Cuts
                float xi3DIpSigValue_;
                float xiPi3DIpSigValue_;
                float VTrkPi3DIpSigValue_;
                float VTrkP3DIpSigValue_;
                float xiFlightSigValue_;
                float distanceSigValue_; //lambda
                //Candidate parameters
                float mass_;
                float pt_;
                int n_;
        };

class XiOmTTree : public edm::EDAnalyzer
{
public:
	explicit XiOmTTree(const edm::ParameterSet&);
	~XiOmTTree();

private:
	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();

  // ----------member data ---------------------------

    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> _xiCollection;
    edm::EDGetTokenT<reco::VertexCollection> _vertexCollName;
    edm::EDGetTokenT<reco::TrackCollection> _trkSrc;

	TTree* XiTree;
	TTree* OmTree;

	XiData xistruct_;
	//OmData omstruct_;

    std::string v0CollName_;
    std::string v0IDName_;
};

#endif
