// -*- C++ -*-
//
// Package:    XiCorrelation
// Class:      XiCorrelation
//
/**\class XiCorrelation XiCorrelation.cc XiAnalyzer/XiCorrelation/src/XiCorrelation.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benjamin Lucky Tran,,,
//         Created:  Sat Nov 12 08:57:22 CET 2016
// $Id$
//
//

//#include "RiceHIG/V0Analysis/interface/XiCorrelation.h"
#include "XiAnalyzer/XiAnalyzer/interface/XiCorrelation.h" // for interactive

#define PI 3.1416

XiCorrelation::XiCorrelation(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    using std::string;

    TH1::SetDefaultSumw2();
    _trkSrc         = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trkSrc"));
    _xiCollection   = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("xiCollection"));
    _vertexCollName = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollName"));
    PtBinNum_       = iConfig.getParameter<int>("PtBinNum");
    bkgnum_         = iConfig.getParameter<double>("bkgnum");
    dorap_          = iConfig.getParameter<bool>("dorap");
    etaMax_ass_     = iConfig.getParameter<double>("etaMax_ass");
    etaMin_ass_     = iConfig.getParameter<double>("etaMin_ass");
    multHigh_       = iConfig.getParameter<int>("multHigh");
    multLow_        = iConfig.getParameter<int>("multLow");
    peakFactor_     = iConfig.getParameter<int>("peakFactor");
    ptBin_          = iConfig.getParameter<std::vector<double> >("ptBin");
    ptMax_ass_      = iConfig.getParameter<double>("ptMax_ass");
    ptMin_ass_      = iConfig.getParameter<double>("ptMin_ass");
    rapMax_         = iConfig.getParameter<double>("rapMax");
    rapMin_         = iConfig.getParameter<double>("rapMin");
    sideFactor_     = iConfig.getParameter<int>("sideFactor");
    xiMassMean_     = iConfig.getParameter<std::vector<double> >("xiMassMean");
    xiMassSigma_    = iConfig.getParameter<std::vector<double> >("xiMassSigma");
    zVtxHigh_       = iConfig.getParameter<double>("zVtxHigh");
    zVtxLow_        = iConfig.getParameter<double>("zVtxLow");
}


XiCorrelation::~XiCorrelation()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
XiCorrelation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(_vertexCollName, vertices);

    double bestvz      = -999, bestvx        = -999, bestvy        = -999;
    double bestvzError = -999.9, bestvxError = -999.9, bestvyError = -999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvx      = vtx.x();
    bestvy      = vtx.y();
    bestvz      = vtx.z();
    bestvxError = vtx.xError();
    bestvyError = vtx.yError();
    bestvzError = vtx.zError();

    if(bestvz > zVtxHigh_ || bestvz < zVtxLow_){
        cout << "Bad Z vertex" << endl;
        return;
    }

    pepVect_trkass = new vector<TVector3>;
    pepVect_Xi     = new vector<TVector3>;
    xiMass         = new vector<double>;

    edm::Handle<reco::VertexCompositeCandidateCollection> xiCollection;
    iEvent.getByToken(_xiCollection, xiCollection);

    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(_trkSrc, tracks);

    if(!xiCollection.isValid()){
        cout << "Collection invalid" << endl;
        return;
    }

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
        if(fabs(trk.eta()) > 2.4 || trk.pt() < 0.4)   continue;
        EtaPtCutnTracks++;
    }


    if(EtaPtCutnTracks >= multLow_ && EtaPtCutnTracks < multHigh_){
        nEvtCut->Fill(1);
        nTrk->Fill(EtaPtCutnTracks);
        for(reco::VertexCompositeCandidateCollection::const_iterator xiCand =
                xiCollection->begin(); xiCand != xiCollection->end(); xiCand++)
        {
            if(dorap_)
            {
                double rap = xiCand->rapidity();
                if(rap < rapMin_ || rap > rapMax_) continue;
            }

            // Make 2D Mass v Pt
            double xi_eta = xiCand->eta();
            double xi_phi = xiCand->phi();
            double mass   = xiCand->mass();
            double xi_pT  = xiCand->pt();
            MassPt->Fill(mass,xi_pT);

            // Make vector of Xi Candidate parameters
            TVector3 xiPEPvector;
            xiPEPvector.SetPtEtaPhi(xi_pT,xi_eta,xi_phi);
            pepVect_Xi->push_back(xiPEPvector);
            xiMass->push_back(mass);
        }

        // Make vectors for pairing of primary charged tracks
        //
        for(unsigned it = 0; it < tracks->size(); it++)
        {
            const reco::Track & trk = (*tracks)[it];
            math::XYZPoint bestvtx(bestvx, bestvy, bestvz);

            double dzvtx    = trk.dz(bestvtx);
            double dxyvtx   = trk.dxy(bestvtx);
            double dzerror  = sqrt(trk.dzError()*trk.dzError() + bestvzError*bestvzError);
            double dxyerror = sqrt(trk.d0Error()*trk.d0Error() + bestvxError*bestvyError);

            if(!trk.quality(reco::TrackBase::highPurity)) continue;
            if(fabs(trk.ptError()/trk.pt() > 0.10))   continue;
            if(fabs(dzvtx/dzerror) > 3)                   continue;
            if(fabs(dxyvtx/dxyerror) > 3)                 continue;

            double eta = trk.eta();
            double phi = trk.phi();
            double pt  = trk.pt();

            // This vector is for the signal histogram pairing of charged
            // primary tracks and with Xi candidates
            TVector3 primPEPvector;
            primPEPvector.SetPtEtaPhi(pt,eta,phi);
            if(eta <= etaMax_ass_  &&
               eta >= etaMin_ass_ &&
               pt <= ptMax_ass_ && //pT associated window
               pt >= ptMin_ass_ //pT associated window
            )
            pepVect_trkass->push_back(primPEPvector);
        }

        // Make signal histogram between Xi candidates and charged primary
        // tracks
        for(int i=0; i<PtBinNum_; i++)
        {
            int pepVect_Xi_size     = (int)pepVect_Xi->size();
            int pepVect_trkass_size = (int)pepVect_trkass->size();
            TrkassPerEvt->Fill(pepVect_trkass_size);

            // PEAK REGION signal
            for(int xi_trg = 0; xi_trg < pepVect_Xi_size; xi_trg++)
            {
                if((*xiMass)[xi_trg] >= (xiMassMean_[i] - peakFactor_*xiMassSigma_[i]) && (*xiMass)[xi_trg] <= (xiMassMean_[i] + peakFactor_*xiMassSigma_[i]))
                {
                    TVector3 pepVect_trg = (*pepVect_Xi)[xi_trg];
                    double eta_trg       = pepVect_trg.Eta();
                    double phi_trg       = pepVect_trg.Phi();
                    double pT            = pepVect_trg.Pt();
                    if(pT >= ptBin_[i] && pT <= ptBin_[i+1]){

                        for(int assoc = 0; assoc < pepVect_trkass_size; assoc++)
                        {
                            TVector3 pepVect_ass = (*pepVect_trkass)[assoc];
                            double eta_ass       = pepVect_ass.Eta();
                            double phi_ass       = pepVect_ass.Phi();

                            double dEta = eta_ass - eta_trg;
                            double dPhi = phi_ass - phi_trg;

                            if(dPhi > PI)
                                dPhi=dPhi-2*PI;
                            if(dPhi < -PI)
                                dPhi=dPhi+2*PI;
                            if(dPhi > -PI && dPhi < -PI/2.0)
                                dPhi=dPhi+2*PI;

                            // To reduce jet fragmentation contributions
                            if(fabs(dEta) < 0.028 && fabs(dPhi) < 0.02) continue;
                            SignalXiPeak[i]->Fill(dEta, dPhi, 1.0/pepVect_Xi_size);
                        }
                    }
                }
            }



            // SIDEBAND REGION signal
            for(int xi_trg = 0; xi_trg < pepVect_Xi_size; xi_trg++)
            {
                if(((*xiMass)[xi_trg] <= (xiMassMean_[i] - sideFactor_*xiMassSigma_[i]) && (*xiMass)[xi_trg] >= 1.25) || ((*xiMass)[xi_trg] <= 1.40 && (*xiMass)[xi_trg] >= (xiMassMean_[i] + sideFactor_*xiMassSigma_[i])))
                {
                    TVector3 pepVect_trg = (*pepVect_Xi)[xi_trg];
                    double eta_trg       = pepVect_trg.Eta();
                    double phi_trg       = pepVect_trg.Phi();
                    double pT = pepVect_trg.Pt();
                    if(pT >= ptBin_[i] && pT <= ptBin_[i+1])
                    {
                        for(int assoc = 0; assoc < pepVect_trkass_size; assoc++)
                        {
                            TVector3 pepVect_ass = (*pepVect_trkass)[assoc];
                            double eta_ass       = pepVect_ass.Eta();
                            double phi_ass       = pepVect_ass.Phi();

                            double dEta = eta_ass - eta_trg;
                            double dPhi = phi_ass - phi_trg;

                            if(dPhi > PI)
                                dPhi=dPhi-2*PI;
                            if(dPhi < -PI)
                                dPhi=dPhi+2*PI;
                            if(dPhi > -PI && dPhi < -PI/2.0)
                                dPhi=dPhi+2*PI;

                            // To reduce jet fragmentation contributions
                            if(fabs(dEta) < 0.028 && fabs(dPhi) < 0.02) continue;
                            SignalXiSide[i]->Fill(dEta, dPhi, 1.0/pepVect_Xi_size);
                        }
                    }
                }
            }
        }

        PepVect2_Xi->push_back(*pepVect_Xi);
        xiMass2->push_back(*xiMass);
        PepVect2_ass->push_back(*pepVect_trkass);
        zvtxVect->push_back(bestvz);

        //
        // Make signal histogram of Xi-h but now with 1 < pT < 3 GeV and no mass
        // window
        //
        int pepVect_Xi_size     = (int)pepVect_Xi->size();
        int pepVect_trkass_size = (int)pepVect_trkass->size();

        for(int xi_trg = 0; xi_trg < pepVect_Xi_size; xi_trg++)
        {
            TVector3 pepVect_trg = (*pepVect_Xi)[xi_trg];
            double eta_trg       = pepVect_trg.Eta();
            double phi_trg       = pepVect_trg.Phi();
            double Pt = pepVect_trg.Pt();
            if(Pt > 1.0 && Pt < 3.0)
            {
                for(int assoc = 0; assoc < pepVect_trkass_size; assoc++)
                {
                    TVector3 pepVect_ass = (*pepVect_trkass)[assoc];
                    double eta_ass       = pepVect_ass.Eta();
                    double phi_ass       = pepVect_ass.Phi();

                    double dEta = eta_ass - eta_trg;
                    double dPhi = phi_ass - phi_trg;

                    if(dPhi > PI)
                        dPhi=dPhi-2*PI;
                    if(dPhi < -PI)
                        dPhi=dPhi+2*PI;
                    if(dPhi > -PI && dPhi < -PI/2.0)
                        dPhi=dPhi+2*PI;

                    // To reduce jet fragmentation contributions
                    if(fabs(dEta) < 0.028 && fabs(dPhi) < 0.02) continue;
                    SignalXiHad->Fill(dEta, dPhi, 1.0/pepVect_Xi_size);
                }
            }
        }

        //
        // Make signal histogram for pairing of two charged primary tracks
        //
        int pepVect_trkhad_size = (int)pepVect_trkass->size();
        HadPerEvt->Fill(pepVect_trkhad_size);

        for(int trktrg1 = 0; trktrg1 < pepVect_trkhad_size; trktrg1++)
        {
            TVector3 pepVect_had1 = (*pepVect_trkass)[trktrg1];
            double eta_trg = pepVect_had1.Eta();
            double phi_trg = pepVect_had1.Phi();

            for(int trktrg2 = 0; trktrg2 < pepVect_trkass_size; trktrg2++)
            {
                if(trktrg2 == trktrg1){
                    continue;
                }
                TVector3 pepVect_had2 = (*pepVect_trkass)[trktrg2];
                double eta_ass = pepVect_had2.Eta();
                double phi_ass = pepVect_had2.Phi();

                double dEta = eta_ass - eta_trg;
                double dPhi = phi_ass - phi_trg;

                if(dPhi > PI)
                    dPhi=dPhi-2*PI;
                if(dPhi < -PI)
                    dPhi=dPhi+2*PI;
                if(dPhi > -PI && dPhi < -PI/2.0)
                    dPhi=dPhi+2*PI;

                // To reduce jet fragmentation contributions
                if(fabs(dEta) < 0.028 && fabs(dPhi) < 0.02) continue;
                SignalHad->Fill(dEta, dPhi, 1.0/pepVect_trkhad_size);
            }
        }

        delete pepVect_Xi;
        delete pepVect_trkass;
        delete xiMass;

    }
}


// ------------ method called once each job just before starting event loop  ------------
void
XiCorrelation::beginJob()
{

    nTrk            = fs->make<TH1D>("nTrk", "nTrk", 250, 0, 250);
    nEvtCut         = fs->make<TH1D>("nEvtCut", "nEvtCut", 10, 0, 10);
    HadPerEvt       = fs->make<TH1D>("HadPerEvent", "Hadrons per Event", 1500, 0, 1500);
    TrkassPerEvt    = fs->make<TH1D>("TrkassPerEvent", "Associated trks per Event", 300,0, 300);

    MassPt          = fs->make<TH2D>("MassPt", "", 150, 1.25, 1.40, 100, 0 ,10);
    BackgroundHad   = fs->make<TH2D>("BackgroundHad", "BkgHad; #Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    SignalHad       = fs->make<TH2D>("SignalHad", "SigHad; #Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    SignalXiHad     = fs->make<TH2D>("SignalXiHad", ";#Delta#eta;#Delta#phi ", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    BackgroundXiHad = fs->make<TH2D>("BackgroundXiHad", ";#Delta#eta;#Delta#phi ", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    for(int i=0; i<PtBinNum_; i++)
    {
        BackgroundXiPeak[i]   = fs->make<TH2D>(Form("BackgroundPeak_pt%d",i), ";#Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        BackgroundXiSide[i]   = fs->make<TH2D>(Form("BackgroundSide_pt%d",i), ";#Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        SignalXiPeak[i]       = fs->make<TH2D>(Form("SignalPeak_pt%d",i) , ";#Delta#eta;#Delta#phi ", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        SignalXiSide[i]       = fs->make<TH2D>(Form("SignalSide_pt%d",i) , ";#Delta#eta;#Delta#phi ", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);


        /*
        Correlation        = fs->make<TH2D>("Correlation", "Correlation", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        CorrelationPeak        = fs->make<TH2D>("CorrelationPeak", "Correlation Peak", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        CorrelationSide        = fs->make<TH2D>("CorrelationSide", "Correlation Side", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        CorrelationHad     = fs->make<TH2D>("CorrelationHad", "CorrelationHad", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    */
    }

    // For Background calculations which must be done in the endJob function
    //
    PepVect2_Xi        = new vector< vector<TVector3> >;
    PepVect2_ass       = new vector< vector<TVector3> >;
    xiMass2            = new vector< vector<double> >;
    zvtxVect           = new vector<double>;



}

// ------------ method called once each job just after ending the event loop  ------------
void
XiCorrelation::endJob()
{
    for(int i=0; i<PtBinNum_; i++)
    {
        // Make background histograms
        // Xi paired with hadron
        int PepVect2_Xi_size = (int)PepVect2_Xi->size();
        int PepVect2_ass_size = (int)PepVect2_ass->size();

        // PEAK REGION Background
        for(int bkgnum = 0; bkgnum<bkgnum_; bkgnum++)
        {
            int ncount = 0;
            for(int nevt_ass=0; nevt_ass<PepVect2_ass_size; nevt_ass++)
            {
                int nevt_trg = gRandom->Integer(PepVect2_Xi_size);
                if(nevt_trg == nevt_ass)
                {
                    nevt_ass--;
                    continue;
                }
                if(fabs((*zvtxVect)[nevt_trg] - (*zvtxVect)[nevt_ass]) > 0.5)
                {
                    nevt_ass--;
                    ncount++;
                    if(ncount > 5000)
                    {
                        nevt_ass++;
                        ncount=0;
                    }
                    continue;
                }

                vector<double> xiMassTmp = (*xiMass2)[nevt_trg];
                vector<TVector3> pepVectTmp_trg = (*PepVect2_Xi)[nevt_trg];
                vector<TVector3> pepVectTmp_ass = (*PepVect2_ass)[nevt_ass];
                int nMult_trg = pepVectTmp_trg.size();
                int nMult_ass = pepVectTmp_ass.size();


                for(int ntrg=0; ntrg<nMult_trg; ntrg++)
                {
                    if(xiMassTmp[ntrg] >= (xiMassMean_[i] - peakFactor_*xiMassSigma_[i]) && xiMassTmp[ntrg] <= (xiMassMean_[i] + peakFactor_*xiMassSigma_[i]))
                    {
                        TVector3 pvectorTmp_trg = pepVectTmp_trg[ntrg];
                        double eta_trg = pvectorTmp_trg.Eta();
                        double phi_trg = pvectorTmp_trg.Phi();
                        double pT = pvectorTmp_trg.Pt();
                        if(pT >= ptBin_[i] && pT <= ptBin_[i+1])
                        {
                            for(int nass=0; nass<nMult_ass; nass++)
                            {
                                TVector3 pvectorTmp_ass = pepVectTmp_ass[nass];
                                double eta_ass = pvectorTmp_ass.Eta();
                                double phi_ass = pvectorTmp_ass.Phi();

                                double dEta = eta_ass - eta_trg;
                                double dPhi = phi_ass - phi_trg;

                                if(dPhi > PI)                    dPhi=dPhi-2*PI;
                                if(dPhi < -PI)                   dPhi=dPhi+2*PI;
                                if(dPhi > -PI && dPhi < -PI/2.0) dPhi=dPhi+2*PI;

                                if(fabs(dPhi) < 0.028 && fabs(dEta) < 0.02) continue;

                                BackgroundXiPeak[i]->Fill(dEta, dPhi, 1.0/nMult_trg);

                            }
                        }
                    }
                }
            }
        }


        // SIDEBAND REGION Background
        for(int bkgnum = 0; bkgnum<bkgnum_; bkgnum++)
        {
            int ncount = 0;
            for(int nevt_ass=0; nevt_ass<PepVect2_ass_size; nevt_ass++)
            {
                int nevt_trg = gRandom->Integer(PepVect2_Xi_size);
                if(nevt_trg == nevt_ass)
                {
                    nevt_ass--;
                    continue;
                }
                if(fabs((*zvtxVect)[nevt_trg] - (*zvtxVect)[nevt_ass]) > 0.5)
                {
                    nevt_ass--;
                    ncount++;
                    if(ncount > 5000)
                    {
                        nevt_ass++;
                        ncount=0;
                    }
                    continue;
                }

                vector<double> xiMassTmp = (*xiMass2)[nevt_trg];
                vector<TVector3> pepVectTmp_trg = (*PepVect2_Xi)[nevt_trg];
                vector<TVector3> pepVectTmp_ass = (*PepVect2_ass)[nevt_ass];
                int nMult_trg = pepVectTmp_trg.size();
                int nMult_ass = pepVectTmp_ass.size();


                for(int ntrg=0; ntrg<nMult_trg; ntrg++)
                {
                    if((xiMassTmp[ntrg] <= (xiMassMean_[i] - sideFactor_*xiMassSigma_[i]) && xiMassTmp[ntrg] >= 1.25) || (xiMassTmp[ntrg] <= 1.40 && xiMassTmp[ntrg] >= (xiMassMean_[i] + sideFactor_*xiMassSigma_[i])))
                    {
                        TVector3 pvectorTmp_trg = pepVectTmp_trg[ntrg];
                        double eta_trg = pvectorTmp_trg.Eta();
                        double phi_trg = pvectorTmp_trg.Phi();
                        double pT = pvectorTmp_trg.Pt();
                        if(pT >= ptBin_[i] && pT <= ptBin_[i+1])
                        {
                            for(int nass=0; nass<nMult_ass; nass++)
                            {
                                TVector3 pvectorTmp_ass = pepVectTmp_ass[nass];
                                double eta_ass = pvectorTmp_ass.Eta();
                                double phi_ass = pvectorTmp_ass.Phi();

                                double dEta = eta_ass - eta_trg;
                                double dPhi = phi_ass - phi_trg;

                                if(dPhi > PI)                    dPhi=dPhi-2*PI;
                                if(dPhi < -PI)                   dPhi=dPhi+2*PI;
                                if(dPhi > -PI && dPhi < -PI/2.0) dPhi=dPhi+2*PI;

                                if(fabs(dPhi) < 0.028 && fabs(dEta) < 0.02) continue;

                                BackgroundXiSide[i]->Fill(dEta, dPhi, 1.0/nMult_trg);
                            }
                        }

                    }

                }
            }
        }
    }

    //
    // Make for 1 < pT < 3 GeV and no mass window
    //

    int PepVect2_Xi_size = (int)PepVect2_Xi->size();
    int PepVect2_ass_size = (int)PepVect2_ass->size();

    for(int bkgnum = 0; bkgnum<bkgnum_; bkgnum++)
    {
        int ncount = 0;
        for(int nevt_ass=0; nevt_ass<PepVect2_ass_size; nevt_ass++)
        {
            int nevt_trg = gRandom->Integer(PepVect2_Xi_size);
            if(nevt_trg == nevt_ass)
            {
                nevt_ass--;
                continue;
            }
            if(fabs((*zvtxVect)[nevt_trg] - (*zvtxVect)[nevt_ass]) > 0.5)
            {
                nevt_ass--;
                ncount++;
                if(ncount > 5000)
                {
                    nevt_ass++;
                    ncount=0;
                }
                continue;
            }

            vector<double> xiMassTmp = (*xiMass2)[nevt_trg];
            vector<TVector3> pepVectTmp_trg = (*PepVect2_Xi)[nevt_trg];
            vector<TVector3> pepVectTmp_ass = (*PepVect2_ass)[nevt_ass];
            int nMult_trg = pepVectTmp_trg.size();
            int nMult_ass = pepVectTmp_ass.size();


            for(int ntrg=0; ntrg<nMult_trg; ntrg++)
            {
                TVector3 pvectorTmp_trg = pepVectTmp_trg[ntrg];
                double eta_trg = pvectorTmp_trg.Eta();
                double phi_trg = pvectorTmp_trg.Phi();
                double pT = pvectorTmp_trg.Pt();
                if(pT >= 1 && pT <= 3)
                {
                    for(int nass=0; nass<nMult_ass; nass++)
                    {
                        TVector3 pvectorTmp_ass = pepVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();

                        double dEta = eta_ass - eta_trg;
                        double dPhi = phi_ass - phi_trg;

                        if(dPhi > PI)                    dPhi=dPhi-2*PI;
                        if(dPhi < -PI)                   dPhi=dPhi+2*PI;
                        if(dPhi > -PI && dPhi < -PI/2.0) dPhi=dPhi+2*PI;

                        if(fabs(dPhi) < 0.028 && fabs(dEta) < 0.02) continue;

                        BackgroundXiHad->Fill(dEta, dPhi, 1.0/nMult_trg);

                    }
                }
            }
        }
    }

    //
    // hadron paired with hadron
    //

    for(int bkgnum = 0; bkgnum<bkgnum_; bkgnum++)
    {
        int ncount = 0;
        for(int nevt_ass=0; nevt_ass<PepVect2_ass_size; nevt_ass++)
        {
            int nevt_trg = gRandom->Integer(PepVect2_ass_size);
            if(nevt_trg == nevt_ass)
            {
                nevt_ass--;
                continue;
            }
            if(fabs((*zvtxVect)[nevt_trg] - (*zvtxVect)[nevt_ass]) > 0.5)
            {
                nevt_ass--;
                ncount++;
                if(ncount > 5000)
                {
                    nevt_ass++;
                    ncount=0;
                }
                continue;
            }

            vector<TVector3> pepVectTmp_trg = (*PepVect2_ass)[nevt_trg];
            vector<TVector3> pepVectTmp_ass = (*PepVect2_ass)[nevt_ass];
            int nMult_trg = pepVectTmp_trg.size();
            int nMult_ass = pepVectTmp_ass.size();

            for(int ntrg=0; ntrg<nMult_trg; ntrg++)
            {
                TVector3 pvectorTmp_trg = pepVectTmp_trg[ntrg];
                double eta_trg = pvectorTmp_trg.Eta();
                double phi_trg = pvectorTmp_trg.Phi();

                for(int nass=0; nass<nMult_ass; nass++)
                {
                    /*For avoiding correlation with the same particle although
                    this would already be accounted for in the continue
                    statement below*/
                    /*
                    if(ntrg == nass) continue;
                    */
                    TVector3 pvectorTmp_ass = pepVectTmp_ass[nass];
                    double eta_ass = pvectorTmp_ass.Eta();
                    double phi_ass = pvectorTmp_ass.Phi();

                    double dEta = eta_ass - eta_trg;
                    double dPhi = phi_ass - phi_trg;

                    if(dPhi > PI)                    dPhi=dPhi-2*PI;
                    if(dPhi < -PI)                   dPhi=dPhi+2*PI;
                    if(dPhi > -PI && dPhi < -PI/2.0) dPhi=dPhi+2*PI;

                    if(fabs(dPhi) < 0.028 && fabs(dEta) < 0.02) continue;

                    BackgroundHad->Fill(dEta, dPhi, 1.0/nMult_trg);

                }
            }
        }
    }

    /*
    //int nEvent = XiPerEvt->Integral(3, 10000);
    //double bz = BackgroundXi->GetBinContent(275);
    int nEventPeak = XiPerEvtPeak->Integral(3, 10000);
    if(nEventPeak == 0) nEventPeak = 1;
    CorrelationPeak->Add(SignalXiPeak);
    //Correlation->Scale(bz);
    CorrelationPeak->Divide(BackgroundXiPeak);
    CorrelationPeak->Scale(1.0/nEventPeak);

    int nEventSide = XiPerEvtSide->Integral(3, 10000);
    if(nEventSide == 0) nEventSide = 1;
    CorrelationSide->Add(SignalXiSide);
    //Correlation->Scale(bz);
    CorrelationSide->Divide(BackgroundXiSide);
    CorrelationSide->Scale(1.0/nEventSide);

    int nEventbkg = HadPerEvt->Integral(3, 10000);
    if(nEventbkg == 0) nEventbkg = 1;
    //double bz = BackgroundXi->GetBinContent(275);
    CorrelationHad->Add(SignalHad);
    //Correlation->Scale(bz);
    CorrelationHad->Divide(BackgroundHad);
    CorrelationHad->Scale(1.0/nEventbkg);
    */


}

