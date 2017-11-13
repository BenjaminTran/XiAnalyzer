// -*- C++ -*-
// Package:    MassPtProducer
// Class:      MassPtProducer
//
/**\class MassPtProducer MassPtProducer.cc XiAnalyzer/MassPtProducer/src/MassPtProducer.cc

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


#include "XiAnalyzer/XiAnalyzer/interface/MassPtProducer.h"

//
// constructors and destructor
//
MassPtProducer::MassPtProducer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

    using std::string;
    TH1::SetDefaultSumw2();

    _ksCollection   = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("ksCollection"));
    _laCollection   = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("laCollection"));
    _trkSrc         = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trkSrc"));
    _vertexCollName = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollName"));
    _xiCollection   = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("xiCollection"));
    _omCollection = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("omCollection"));

    multHigh_ = iConfig.getParameter<double>("multHigh");
    multLow_  = iConfig.getParameter<double>("multLow");
    zVtxHigh_ = iConfig.getParameter<double>("zVtxHigh");
    zVtxLow_  = iConfig.getParameter<double>("zVtxLow");
    rapMin_ = iConfig.getParameter<double>("rapMin");
    rapMax_ = iConfig.getParameter<double>("rapMax");

    ks_ = iConfig.getUntrackedParameter<bool>("ks",false);
    la_ = iConfig.getUntrackedParameter<bool>("la",false);
    xi_ = iConfig.getUntrackedParameter<bool>("xi",false);
    om_ = iConfig.getUntrackedParameter<bool>("om",false);
    MC_ = iConfig.getUntrackedParameter<bool>("MC",false);

    if(MC_)
    {
        _gnCollection   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gnCollection"));
    }
}


MassPtProducer::~MassPtProducer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MassPtProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    nEvt->Fill(1);
    using namespace edm;
    int EtaPtCutnTracks = 0;

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


    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(_trkSrc, tracks);

    edm::Handle<reco::VertexCompositeCandidateCollection> xiCollection;
    if(xi_)
    {
        iEvent.getByToken(_xiCollection, xiCollection);
        if(!xiCollection.isValid())
        {
            cout << "Xi Collection invalid" << endl;
            return;
        }
    }

    edm::Handle<reco::VertexCompositeCandidateCollection> omCollection;
    if(om_)
    {
        iEvent.getByToken(_omCollection, omCollection);
        if(!omCollection.isValid()) 
        {
            cout << "Om Collection invalid" << endl;
            return;
        }
    }

    edm::Handle<reco::VertexCompositeCandidateCollection> ksCollection;
    if(ks_)
    {
        iEvent.getByToken(_ksCollection, ksCollection);
        if(!ksCollection.isValid()) 
        {
            cout << "Ks Collection invalid" << endl;
            return;
        }
    }

    edm::Handle<reco::VertexCompositeCandidateCollection> laCollection;
    if(la_)
    {
        iEvent.getByToken(_laCollection, laCollection);
        if(!laCollection.isValid()) 
        {
            cout << "La Collection invalid" << endl;
            return;
        }
    }


    // Track selection
    //if((ks_ && ksCollection->size() != 0) || (la_ && laCollection->size() != 0) || (xi_ && xiCollection->size() != 0) || om_ && omCollection->size() != 0) //If all collection sizes are zero then skip looping over tracks and exit function to save time
    //{
        int nTracks         = 0;
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
        nTrk->Fill(nTracks); //number of acceptable tracks
    //}


    if(EtaPtCutnTracks >= multLow_ && EtaPtCutnTracks < multHigh_){
        nEvtCut->Fill(1); //number of events that pass the multiplicity cut
        EtaPtCutnTrackHist->Fill(EtaPtCutnTracks); //number of tracks in the current event that passed multiplicity requirement. Should only be between multlow and multhigh
        //XI
        if(xi_ && xiCollection.isValid())
        {
            for(reco::VertexCompositeCandidateCollection::const_iterator xiCand =
                    xiCollection->begin(); xiCand != xiCollection->end(); xiCand++) {

                double rap_xi  = xiCand->rapidity();
                double mass_xi = xiCand->mass();
                double pT_xi   = xiCand->pt();
                double eta_xi  = xiCand->eta();

                XiMassPtRap       -> Fill(mass_xi,pT_xi,rap_xi);
                rapidity_xi       -> Fill(rap_xi);
                pseudorapidity_xi -> Fill(eta_xi);
            }
        }
        //OM
        if(om_ && omCollection.isValid())
        {
            for(reco::VertexCompositeCandidateCollection::const_iterator omCand =
                    omCollection->begin(); omCand != omCollection->end(); omCand++) {
                double rap_om  = omCand->rapidity();
                double mass_om = omCand->mass();
                double pT_om   = omCand->pt();
                double eta_om  = omCand->eta();
                OmMassPtRap       -> Fill(mass_om,pT_om,rap_om);
                rapidity_om       -> Fill(rap_om);
                pseudorapidity_om -> Fill(eta_om);
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

                KsMassPtRap       -> Fill(mass_ks,pT_ks,rap_ks);
                rapidity_ks       -> Fill(rap_ks);
                pseudorapidity_ks -> Fill(eta_ks);
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

                LaMassPtRap       -> Fill(mass_la,pT_la,rap_la);
                rapidity_la       -> Fill(rap_la);
                pseudorapidity_la -> Fill(eta_la);
            }
        }

        if(MC_)
        {
            edm::Handle<reco::GenParticleCollection> gnCollection;
            iEvent.getByToken(_gnCollection, gnCollection);

            for(reco::GenParticleCollection::const_iterator gnCand = gnCollection->begin(); gnCand != gnCollection->end(); gnCand++)
            {
                int id = gnCand->pdgId();
                int st = gnCand->status();
                double rapidity_gn = gnCand->rapidity();
                double eta_gn = gnCand->eta();
                double pt_gn = gnCand->pt();
                double mass_gn = gnCand->mass();

                if(st != 1) continue;
                if(rapidity_gn > rapMin_ && rapidity_gn < rapMax_)
                {
                    //Lambda and mother identification
                    int mid = 0;
                    if(TMath::Abs(id) == 3122)
                    {
                        if(gnCand->numberOfMothers() == 1)
                        {
                            const reco::Candidate* mom = gnCand->mother();
                            mid = mom->pdgId();
                            if(mom->numberOfMothers()==1)
                            {
                                const reco::Candidate* mom1 = mom->mother();
                                mid = mom1->pdgId();
                            }
                        }
                        //make sure lambda isnt from decay channel
                        if(TMath::Abs(mid) != 3322 && TMath::Abs(mid) != 3312 && TMath::Abs(mid) != 3324 && TMath::Abs(mid) != 3314 && TMath::Abs(mid) != 3334)
                            LaMassPtRap_Gen->Fill(mass_gn,pt_gn,rapidity_gn);
                    }

                    //KShort
                    if(TMath::Abs(id) == 310)
                    {
                        KsMassPtRap_Gen->Fill(mass_gn,pt_gn,rapidity_gn);
                    }

                    //Cascade
                    int midXi = 0;
                    if(TMath::Abs(id) == 3312)
                    {
                        if(gnCand->numberOfMothers() == 1)
                        {
                            const reco::Candidate* mom = gnCand->mother();
                            midXi = mom->pdgId();
                            if(mom->numberOfMothers() == 1)
                            {
                                const reco::Candidate* mom1 = mom->mother();
                                midXi = mom1->pdgId();
                            }
                        }
                        if(TMath::Abs(midXi) != 3334)
                            XiMassPtRap_Gen->Fill(mass_gn,pt_gn,rapidity_gn);
                    }

                    //Omega
                    int midOm = 0;
                    if(TMath::Abs(id) == 3334)
                    {
                        //if(gnCand->numberOfMothers() == 1)
                        //{
                            //const reco::Candidate* mom = gnCand->mother();
                            //midOm = mom->pdgId();
                            //if(mom->numberOfMothers() == 1)
                            //{
                                //const reco::Candidate* mom1 = mom->mother();
                                //midOm = mom1->pdgId();
                            //}
                        //}
                        OmMassPtRap_Gen->Fill(mass_gn,pt_gn,rapidity_gn);
                    }
                }
            }
        }
    }

}


// ------------ method called once each job just before starting event loop  ------------
void
MassPtProducer::beginJob()
{
    if(xi_) cout << "Will Access Xi" << endl;
    if(ks_) cout << "Will Access Ks" << endl;
    if(la_) cout << "Will Access La" << endl;
    if(om_) cout << "Will Access Om" << endl;
    if(MC_) cout << "Will Access MC" << endl;

    XiMassPtRap        = fs->make<TH3D>("XiMassPtRap", "#Xi Mass, Pt, y", 150, 1.25, 1.40, 400, 0, 40,22,-1.1,1.1);
    LaMassPtRap        = fs->make<TH3D>("LaMassPtRap", "#Lambda Mass, Pt, y", 160, 1.08, 1.160, 400, 0, 40,22,-1.1,1.1);
    KsMassPtRap        = fs->make<TH3D>("KsMassPtRap", "Ks Mass, Pt, y", 270, 0.43, 0.565, 400, 0, 40,22,-1.1,1.1);
    OmMassPtRap = fs->make<TH3D>("OmMassPt", "Om Mass, Pt, y",150,1.60,1.75,400,0,40,22,-1.1,1.1);
    nEvt               = fs->make<TH1D>("nEvt","nEvt",10,0,10);
    nTrk               = fs->make<TH1D>("nTrk", "nTrk", 400, 0, 400);
    nEvtCut            = fs->make<TH1D>("nEvtCut", "nEvtCut", 10,0,10);
    EtaPtCutnTrackHist = fs->make<TH1D>("EtaPtCutnTrackHist", "EtaPtCutnTrack",250,0,250);
    rapidity_xi        = fs->make<TH1D>("XiRapidity","XiRapidity",200,-10,10);
    rapidity_om        = fs->make<TH1D>("OmRapidity","OmRapidity",200,-10,10);
    rapidity_ks        = fs->make<TH1D>("KsRapidity","KsRapidity",200,-10,10);
    rapidity_la        = fs->make<TH1D>("LaRapidity","LaRapidity",200,-10,10);
    pseudorapidity_xi  = fs->make<TH1D>("XiEta","XiEta",200,-10,10);
    pseudorapidity_om  = fs->make<TH1D>("OmEta","OmEta",200,-10,10);
    pseudorapidity_ks  = fs->make<TH1D>("KsEta","KsEta",200,-10,10);
    pseudorapidity_la  = fs->make<TH1D>("LaEta","LaEta",200,-10,10);

    if(MC_)
    {
        XiMassPtRap_Gen = fs->make<TH3D>("XiMassPtRap_Gen", "#Xi Mass, Pt, y",150,1.25,1.40,400,0,40,22,-1.1,1.1);
        LaMassPtRap_Gen = fs->make<TH3D>("LaMassPtRap_Gen", "#Lambda Mass, Pt, y",160,1.08,1.160,400,0,40,22,-1.1,1.1);
        KsMassPtRap_Gen = fs->make<TH3D>("KsMassPtRap_Gen", "Ks Mass and Pt",270,0.43,0.565,400,0,40,22,-1.1,1.1);
        OmMassPtRap_Gen = fs->make<TH3D>("OmMassPtRap_Gen", "Omega Mass and Pt",150,1.60,1.75,400,0,40,22,-1.1,1.1);
    }
}

// ------------ method called once each job just after ending the event loop  ------------
void
MassPtProducer::endJob()
{
}

/*
// ------------ method called when starting to processes a run  ------------
void
MassPtProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
MassPtProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
MassPtProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
MassPtProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MassPtProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
*/

