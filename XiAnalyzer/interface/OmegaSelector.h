// -*- C++ -*-
//
// Package:    OmegaSelector
// Class:      OmegaSelector
//
/**\class OmegaSelector OmegaSelector.h RiceHIG/V0Analysis/interface/OmegaSelector.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author: Benjamin Tran
//
//

#ifndef XIANALYZER__OMEGA_SELECTOR_H
#define XIANALYZER__OMEGA_SELECTOR_H

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

class OmegaSelector : public edm::EDProducer {
public:
  explicit OmegaSelector(const edm::ParameterSet&);
  ~OmegaSelector();

private:
  //virtual void beginJob() ;
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::EDGetTokenT<reco::VertexCollection> _vertexCollName;
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> _OmegaCollection;
  std::string v0CollName_;
  std::string v0IDName_;
  double etaCutMin_,etaCutMax_;
  double zVertexLow_, zVertexHigh_;
  double xi3DIpSigValue_;
  double xiPi3DIpSigValue_;
  double VTrkPi3DIpSigValue_;
  double VTrkP3DIpSigValue_;
  double xiFlightSigValue_;
  double distanceSigValue_;
  double vtxChi2Cut_;
  double cosThetaCut_;
  double misIDMassCut_;
  double misIDMassCutEE_;
};

#endif
