// -*- C++ -*-
//
// Package:    XiCorrelation.h
// Class:      XiCorrelation.h
//
/**\class XiCorrelation XiCorrelation.h
 * RiceHIG/V0Analysis/interface/XiCorrelation.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author: Benjamin Tran
//
//

//#ifndef RICEHIG__XI_CORRELATION_H
//#define RICEHIG__XI_CORRELATION_H
// For interactive
#ifndef XIANALYZER__XI_CORRELATION_H
#define XIANALYZER__XI_CORRELATION_H

// System include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>
#include <sstream>

// user include files
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
#include <TLorentzVector.h>

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

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
//#include "Geometry/TrackerGeometryBuilder/interface/TrackerLayerIdAccessor.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h" //TrackerLayerIDAccessor.h has migrated to TrackerTopolgy.h as of 6x

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"

//#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h" // Don't know why but this has been moved from interface to plugin
#include "L1Trigger/GlobalTrigger/plugins/L1GlobalTrigger.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
//#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h" //deprecated?
#include "SimTracker/TrackAssociatorProducers/plugins/TrackAssociatorByHitsImpl.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// Heavyion
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

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

using namespace std;

class XiCorrelation : public edm::EDAnalyzer {
    public :
        explicit XiCorrelation(const edm::ParameterSet&);
        ~XiCorrelation();

    private :
        virtual void beginJob();
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob();

        //virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        //virtual void endRun(edm::Run const&, edm::EventSetup const&);
        //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        //virtual void fillDescriptions(edm::ConfigurationDescriptions& descriptions const&);

        edm::Service<TFileService> fs;

        edm::EDGetTokenT<reco::TrackCollection> _trkSrc;
        edm::EDGetTokenT<reco::VertexCollection> _vertexCollName;
        edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> _xiCollection;

        double bkgnum_;
        double etaMax_ass_;
        double etaMin_ass_;
        double ptMax_ass_;
        double ptMin_ass_;
        double xiMassHigh_;
        double xiMassLow_;
        double zVtxHigh_;
        double zVtxLow_;
        int PtBinNum_;
        int multHigh_;
        int multLow_;
        int peakFactor_;
        int sideFactor_;
        bool doRap_;

        //vector<TVector3> *pepVect_trkhad;
        vector<TVector3> *pepVect_trkass;
        vector<TLorentzVector> *pepVect_Xi_peak[11];
        vector<TLorentzVector> *pepVect_Xi_side[11];
        vector<TVector3> *pepVect_dau_xi_peak[11];
        vector<TVector3> *pepVect_dau_xi_side[11];
        vector<double> ptBin_;
        vector<double> xiMassMean_;
        vector<double> xiMassSigma_;
        vector<double> *zvtxVect;
        vector< vector<TLorentzVector> > *PepVect2_Xi_peak[11];
        vector< vector<TLorentzVector> > *PepVect2_Xi_side[11];
        vector<vector<TVector3> > *PepVect2_dau_xi_peak[11];
        vector<vector<TVector3> > *PepVect2_dau_xi_side[11];
        vector< vector<TVector3> > *PepVect2_ass;
        //vector< vector<TVector3> > *PepVect2_had;


        TH1D* nTrk;
        TH1D* EtaPtCutnTrackHist;
        TH1D* nEvtCut;
        TH1D* HadPerEvt;
        TH1D* TrkassPerEvt;

        TH1D* nEvt;

        TH2D* MassPt;
        TH2D* BackgroundXiPeak[11];
        TH2D* BackgroundXiSide[11];
        TH2D* BackgroundHad;
        TH2D* BackgroundXiHad;
        TH2D* SignalXiPeak[11];
        TH2D* SignalXiSide[11];
        TH2D* SignalHad;
        TH2D* SignalXiHad;

        TH1D* KET_xi[11];
        TH1D* KET_xi_bkg[11];
        TH1D* Mass_xi[11];
        TH1D* Pt_xi[11];
        TH1D* Pt_xi_bkg[11];
        TH1D* Eta_xi[11];
        TH1D* Eta_xi_bkg[11];
        TH1D* rap_xi[11];
        TH1D* rap_xi_bkg[11];
        TH1D* rap_xi_Lorentz[11];
        TH1D* mult_xi[11];
        TH1D* mult_xi_bkg[11];

        TH2D* effhisto_xi;

};

#endif
