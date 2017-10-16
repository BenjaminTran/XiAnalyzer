// -*- C++ -*-
//
// Package:    V0Selector
// Class:      V0Selector
// 
/**\class V0Selector V0Selector.h RiceHIG/V0Analysis/interface/V0Selector.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li
//
//

#ifndef RICEHIG__V0_SELECTOR_H
#define RICEHIG__V0_SELECTOR_H

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
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <TString.h>
#include <TVector3.h>
#include <TMatrixD.h>
  
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

class V0Selector : public edm::EDProducer {
public:
  explicit V0Selector(const edm::ParameterSet&);
  ~V0Selector();

private:
  //virtual void beginJob() ;
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::EDGetTokenT<reco::VertexCollection> _vertexCollName;
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> _V0Collection;

  edm::InputTag vertexCollName_;
  std::string v0CollName_;
  std::string v0IDName_;
  double ptCut1_,ptCut2_;
  int nHitCut1_,nHitCut2_;
  double etaCutMin_,etaCutMax_;
  double dxySigCut1_,dxySigCut2_,dzSigCut1_,dzSigCut2_;
  double vtxChi2Cut_;
  double cosThetaCut_;
  double cosThetaCut0_;
  double decayLSigCut_;
  double misIDMassCut_;
  double misIDMassCutEE_;
  bool doRap_;
  double rapMax_;
  double rapMin_;
};

#endif
