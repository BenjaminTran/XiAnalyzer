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
namespace Utils
{
    bool VtxCheckIsGood(edm::Handle<reco::VertexCollection> vertices, double zVtxHigh_, double zVtxLow_)
    {
        double bestvz = -999;
        const reco::Vertex & vtx = (*vertices)[0];
        bestvz = vtx.z();
        if(bestvz > zVtxHigh_ || bestvz < zVtxLow_){
            cout << "Bad Z vertex" << endl;
            return false;
        }
        else
            return true;
    }

    int TrackSelection(edm::Handle<reco::TrackCollection> tracks, edm::Handle<reco::VertexCollection> vertices)
    {
        double bestvz      = -999, bestvx        = -999, bestvy        = -999;
        double bestvzError = -999.9, bestvxError = -999.9, bestvyError = -999.9;
        const reco::Vertex & vtx = (*vertices)[0];
        bestvx      = vtx.x();
        bestvy      = vtx.y();
        bestvz      = vtx.z();
        bestvxError = vtx.xError();
        bestvyError = vtx.yError();
        bestvzError = vtx.zError();

        int EtaPtCutnTracks = 0;
        for(unsigned it = 0; it < tracks->size(); it++)
        {
            const reco::Track & trk = (*tracks)[it];
            math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

            double dzvtx    = trk.dz(bestvtx);
            double dxyvtx   = trk.dxy(bestvtx);
            double dzerror  = sqrt(trk.dzError()*trk.dzError() + bestvzError*bestvzError);
            double dxyerror = sqrt(trk.d0Error()*trk.d0Error() + bestvxError*bestvyError);

            if(!trk.quality(reco::TrackBase::highPurity)) continue;
            if(fabs(trk.ptError()/trk.pt() > 0.10))       continue;
            if(fabs(dzvtx/dzerror) > 3)                   continue;
            if(fabs(dxyvtx/dxyerror) > 3)                 continue;
            if(fabs(trk.eta()) > 2.4 || trk.pt() < 0.4)   continue;
            EtaPtCutnTracks++;
        }
        return EtaPtCutnTracks;
    }

}
