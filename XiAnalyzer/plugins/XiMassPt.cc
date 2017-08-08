// -*- C++ -*-
//
// Package:    XiMassPt
// Class:      XiMassPt
//
/**\class XiMassPt XiMassPt.cc XiAnalyzer/XiMassPt/src/XiMassPt.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benjamin Lucky Tran,,,
//         Created:  Wed Apr 19 19:33:04 CEST 2017
// $Id$
//
//


#include "XiAnalyzer/XiAnalyzer/interface/XiMassPt.h"

//
// constructors and destructor
//
XiMassPt::XiMassPt(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

    using std::string;
    TH1::SetDefaultSumw2();

    _ksCollection   = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("ksCollection"));
    _laCollection   = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("laCollection"));
    _trkSrc         = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trkSrc"));
    _vertexCollName = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollName"));
    _xiCollection   = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("xiCollection"));

    multHigh_ = iConfig.getParameter<double>("multHigh");
    multLow_  = iConfig.getParameter<double>("multLow");
    zVtxHigh_ = iConfig.getParameter<double>("zVtxHigh");
    zVtxLow_  = iConfig.getParameter<double>("zVtxLow");

    ks_    = iConfig.getUntrackedParameter<bool>("ks",false);
    la_    = iConfig.getUntrackedParameter<bool>("la",false);
    xi_    = iConfig.getUntrackedParameter<bool>("xi",false);

}


XiMassPt::~XiMassPt()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
XiMassPt::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(_vertexCollName, vertices);

    double bestvz      = -999, bestvx        = -999, bestvy        = -999;
    double bestvzError = -999.9, bestvxError = -999.9, bestvyError = -999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvx = vtx.x();
    bestvy = vtx.y();
    bestvz = vtx.z();
    bestvxError = vtx.xError();
    bestvyError = vtx.yError();
    bestvzError = vtx.zError();

    if(bestvz > zVtxHigh_ || bestvz < zVtxLow_){
        cout << "Bad zvtx" << endl;
        return;
    }


    edm::Handle<reco::VertexCompositeCandidateCollection> xiCollection;
    iEvent.getByToken(_xiCollection, xiCollection);

    edm::Handle<reco::VertexCompositeCandidateCollection> ksCollection;
    iEvent.getByToken(_ksCollection, ksCollection);

    edm::Handle<reco::VertexCompositeCandidateCollection> laCollection;
    iEvent.getByToken(_laCollection, laCollection);

    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(_trkSrc, tracks);

    int nTracks         = 0;
    int EtaPtCutnTracks = 0;

    // Track selection
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

        nTracks++;
        if(fabs(trk.eta()) > 2.4 || trk.pt() < 0.4) continue;
        EtaPtCutnTracks++;
    }

    nTrk->Fill(nTracks);

    if(EtaPtCutnTracks >= multLow_ && EtaPtCutnTracks < multHigh_){
        nEvtCut->Fill(1);
        EtaPtCutnTrackHist->Fill(EtaPtCutnTracks);
        //XI
        if(!xiCollection.isValid()) cout << "Bad XiCollection" << endl;
        if(xi_ && xiCollection.isValid())
        {
            for(reco::VertexCompositeCandidateCollection::const_iterator xiCand =
                    xiCollection->begin(); xiCand != xiCollection->end(); xiCand++) {

                double rap_xi  = xiCand->rapidity();
                double mass_xi = xiCand->mass();
                double pT_xi   = xiCand->pt();
                double eta_xi  = xiCand->eta();

                XiMassPt         ->Fill(mass_xi,pT_xi);
                rapidity_xi      ->Fill(rap_xi);
                pseudorapidity_xi->Fill(eta_xi);

                cout<<"Fill Xi"<<endl;
            }
        }
        //KS
        if(ks_ && ksCollection.isValid())
        {
            for(reco::VertexCompositeCandidateCollection::const_iterator ksCand =
                    ksCollection->begin(); ksCand != ksCollection->end(); ksCand++) {

                double rap_ks  = ksCand->rapidity();
                double mass_ks = ksCand->mass();
                double pT_ks   = ksCand->pt();
                double eta_ks  = ksCand->eta();

                KsMassPt         -> Fill(mass_ks,pT_ks);
                rapidity_ks      -> Fill(rap_ks);
                pseudorapidity_ks-> Fill(eta_ks);

                cout<<"Fill Ks"<<endl;
            }
        }
        //LAMBDA
        if(la_ && laCollection.isValid())
        {
            for(reco::VertexCompositeCandidateCollection::const_iterator laCand =
                    laCollection->begin(); laCand != laCollection->end(); laCand++) {

                double rap_la  = laCand->rapidity();
                double mass_la = laCand->mass();
                double pT_la   = laCand->pt();
                double eta_la  = laCand->eta();

                LaMassPt         -> Fill(mass, pT_la);
                rapidity_la      -> Fill(rap_ks);
                pseudorapidity_la-> Fill(eta_la);

                cout<<"Fill La"<<endl;
            }
        }
    }
}


// ------------ method called once each job just before starting event loop  ------------
void
XiMassPt::beginJob()
{
    if(xi_) cout << "Will Access Xi" << endl;
    if(ks_) cout << "Will Access Ks" << endl;
    if(la_) cout << "Will Access La" << endl;

    XiMassPt           = fs->make<TH2D>("MassPt", "#Xi Mass and Pt", 150, 1.25, 1.40, 400, 0, 40);
    LaMassPt           = fs->make<TH2D>("LaMassPt", "#Lambda Mass and Pt", 160, 1.08, 1.160, 400, 0, 40);
    KsMassPt           = fs->make<TH2D>("KsMassPt", "Ks Mass and Pt", 270, 0.43, 0.565, 400, 0, 40);
    nTrk               = fs->make<TH1D>("nTrk", "nTrk", 400, 0, 400);
    nEvtCut            = fs->make<TH1D>("nEvtCut", "nEvtCut", 10,0,10);
    EtaPtCutnTrackHist = fs->make<TH1D>("EtaPtCutnTrackHist", "EtaPtCutnTrack",250,0,250);
    rapidity_xi        = fs->make<TH1D>("Rapidity","Rapidity",200,-10,10);
    rapidity_ks        = fs->make<TH1D>("Rapidity","Rapidity",200,-10,10);
    rapidity_la        = fs->make<TH1D>("Rapidity","Rapidity",200,-10,10);
    pseudorapidity_xi  = fs->make<TH1D>("eta","eta",200,-10,10);
    pseudorapidity_ks  = fs->make<TH1D>("eta","eta",200,-10,10);
    pseudorapidity_la  = fs->make<TH1D>("eta","eta",200,-10,10);
}

// ------------ method called once each job just after ending the event loop  ------------
void
XiMassPt::endJob()
{
}

/*
// ------------ method called when starting to processes a run  ------------
void
XiMassPt::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
XiMassPt::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
XiMassPt::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
XiMassPt::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
XiMassPt::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
*/

