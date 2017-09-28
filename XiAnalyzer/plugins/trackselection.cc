#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
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

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

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
