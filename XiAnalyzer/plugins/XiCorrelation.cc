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
    _casCollection   = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("casCollection"));
    _vertexCollName = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollName"));
    PtBinNum_       = iConfig.getParameter<int>("PtBinNum");
    bkgnum_         = iConfig.getParameter<double>("bkgnum");
    etaMax_ass_     = iConfig.getParameter<double>("etaMax_ass");
    etaMin_ass_     = iConfig.getParameter<double>("etaMin_ass");
    multHigh_       = iConfig.getParameter<double>("multHigh");
    multLow_        = iConfig.getParameter<double>("multLow");
    peakFactor_     = iConfig.getParameter<int>("peakFactor");
    ptBin_          = iConfig.getParameter<std::vector<double> >("ptBin");
    ptMax_ass_      = iConfig.getParameter<double>("ptMax_ass");
    ptMin_ass_      = iConfig.getParameter<double>("ptMin_ass");
    sideFactor_     = iConfig.getParameter<int>("sideFactor");
    xiMassMean_     = iConfig.getParameter<std::vector<double> >("xiMassMean");
    xiMassSigma_    = iConfig.getParameter<std::vector<double> >("xiMassSigma");
    zVtxHigh_       = iConfig.getParameter<double>("zVtxHigh");
    zVtxLow_        = iConfig.getParameter<double>("zVtxLow");
    doRap_          = iConfig.getParameter<bool>("doRap");
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

    nEvt->Fill(1);

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
    for(int i=0; i<PtBinNum_; i++)
    {
        pepVect_Xi_peak[i]    = new vector<TLorentzVector>;
        pepVect_Xi_side[i]    = new vector<TLorentzVector>;
        pepVect_dau_xi_peak[i]= new vector<TVector3>;
        pepVect_dau_xi_side[i]= new vector<TVector3>;
    }

    edm::Handle<reco::VertexCompositeCandidateCollection> casCollection;
    iEvent.getByToken(_casCollection, casCollection);

    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(_trkSrc, tracks);

    if(!casCollection.isValid()){
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
        //Make vector of cascade PEP
        for(reco::VertexCompositeCandidateCollection::const_iterator xiCand =
                casCollection->begin(); xiCand != casCollection->end(); xiCand++)
        {
            // Make 2D Mass v Pt
            double xi_eta = xiCand->eta();
            double xi_phi = xiCand->phi();
            double mass   = xiCand->mass();
            double xi_pT  = xiCand->pt();
            double xi_rap = xiCand->rapidity();
            MassPt->Fill(mass,xi_pT);
            double Ket = sqrt(mass*mass + xi_pT*xi_pT) - mass;
            double EffXchoice = 0;

            const reco::Candidate *dau1 = xiCand->daughter(0);
            const reco::Candidate *dau2 = xiCand->daughter(1);

            double pt_dau1 = dau1->pt();
            double eta_dau1 = dau1->eta();
            double phi_dau1 = dau1->phi();

            double pt_dau2 = dau2->pt();
            double eta_dau2 = dau2->eta();
            double phi_dau2 = dau2->phi();

            TVector3 dau1PEPvector;
            dau1PEPvector.SetPtEtaPhi(pt_dau1,eta_dau1,phi_dau1);

            TVector3 dau2PEPvector;
            dau2PEPvector.SetPtEtaPhi(pt_dau2,eta_dau2,phi_dau2);

            if(doRap_)
                EffXchoice = xi_rap;
            else
                EffXchoice = xi_eta;

            //efficiency
            //double effxi = effhisto_xi->GetBinContent(effhisto_xi->FindBin(EffXchoice,xi_pT));

            // Make vector of Xi Candidate parameters
            TLorentzVector xiPEPvector;
            xiPEPvector.SetPtEtaPhiE(xi_pT,xi_eta,xi_phi,xi_rap);
            for(int i=0; i<PtBinNum_;i++)
            {
                if(xi_pT <= ptBin_[i+1] && xi_pT >= ptBin_[i])
                {
                    Mass_xi[i]->Fill(mass);
                    //peak
                    if(mass >= (xiMassMean_[i] - peakFactor_*xiMassSigma_[i]) && mass <= (xiMassMean_[i] + peakFactor_*xiMassSigma_[i]))
                    {
                        pepVect_Xi_peak[i]->push_back(xiPEPvector);
                        pepVect_dau_xi_peak[i]->push_back(dau1PEPvector);
                        pepVect_dau_xi_peak[i]->push_back(dau2PEPvector);
                        KET_xi[i]->Fill(Ket);//,1.0/effxi);
                        Pt_xi[i]->Fill(xi_pT);//,1.0/effxi);
                        Eta_xi[i]->Fill(xi_eta);//,1.0/effxi);
                        rap_xi[i]->Fill(xi_rap);//,1.0/effxi);
                        rap_xi_Lorentz[i]->Fill(xiPEPvector.E());//,1.0/effxi);
                    }
                    //sideband
                    if((mass <= (xiMassMean_[i] - sideFactor_*xiMassSigma_[i]) && mass >= 1.25) || (mass <= 1.40 && mass >= (xiMassMean_[i] + sideFactor_*xiMassSigma_[i])))
                    {
                        pepVect_Xi_side[i]->push_back(xiPEPvector);
                        pepVect_dau_xi_side[i]->push_back(dau1PEPvector);
                        pepVect_dau_xi_side[i]->push_back(dau2PEPvector);
                        KET_xi_bkg[i]->Fill(Ket);//,1.0/effxi);
                        Pt_xi_bkg[i]->Fill(xi_pT);//,1.0/effxi);
                        Eta_xi_bkg[i]->Fill(xi_eta);//,1.0/effxi);
                        rap_xi_bkg[i]->Fill(xi_rap);//,1.0/effxi);
                    }
                }
            }
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
            int pepVect_Xi_peak_size     = (int)pepVect_Xi_peak[i]->size();
            int pepVect_Xi_side_size     = (int)pepVect_Xi_side[i]->size();
            int pepVect_trkass_size = (int)pepVect_trkass->size();
            TrkassPerEvt->Fill(pepVect_trkass_size);

            double nMult_trg_eff_xi = 0;

            for(int xi_ntrg = 0; xi_ntrg<pepVect_Xi_peak_size; xi_ntrg++)
            {
                TLorentzVector pepVect_trg = (*pepVect_Xi_peak[i])[xi_ntrg];
                double pt_trg     = pepVect_trg.Pt();
                double eta_trg    = pepVect_trg.Eta();
                double rap_trg    = pepVect_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                //double effxi = effhisto_xi->GetBinContent(effhisto_xi->FindBin(EffXchoice,pt_trg));

                //nMult_trg_eff_xi = nMult_trg_eff_xi + 1.0/effxi;
            }

            mult_xi[i]->Fill(nMult_trg_eff_xi);

            // PEAK REGION signal
            for(int xi_trg = 0; xi_trg < pepVect_Xi_peak_size; xi_trg++)
            {
                TLorentzVector pepVect_trg = (*pepVect_Xi_peak[i])[xi_trg];
                double pt_trg     = pepVect_trg.Pt();
                double eta_trg    = pepVect_trg.Eta();
                double phi_trg    = pepVect_trg.Phi();
                double rap_trg    = pepVect_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                //double effxi = effhisto_xi->GetBinContent(effhisto_xi->FindBin(EffXchoice,pt_trg));

                TVector3 pepvector_trg_dau1 = (*pepVect_dau_xi_peak[i])[2*xi_trg];

                double eta_trg_dau1 = pepvector_trg_dau1.Eta();
                double phi_trg_dau1 = pepvector_trg_dau1.Phi();

                TVector3 pepvector_trg_dau2 = (*pepVect_dau_xi_peak[i])[2*xi_trg+1];

                double eta_trg_dau2 = pepvector_trg_dau2.Eta();
                double phi_trg_dau2 = pepvector_trg_dau2.Phi();

                for(int assoc = 0; assoc < pepVect_trkass_size; assoc++)
                {
                    TVector3 pepVect_ass = (*pepVect_trkass)[assoc];
                    double eta_ass       = pepVect_ass.Eta();
                    double phi_ass       = pepVect_ass.Phi();

                    if(fabs(eta_ass-eta_trg_dau1) < 0.03 && fabs(phi_ass-phi_trg_dau1) < 0.03) continue;
                    if(fabs(eta_ass-eta_trg_dau2) < 0.03 && fabs(phi_ass-phi_trg_dau2) < 0.03) continue;

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
                    SignalXiPeak[i]->Fill(dEta, dPhi);//,1.0/nMult_trg_eff_xi/effxi);
                }
            }

            nMult_trg_eff_xi = 0;

            for(int xi_trg = 0; xi_trg<pepVect_Xi_side_size; xi_trg++)
            {
                TLorentzVector pepVect_trg = (*pepVect_Xi_side[i])[xi_trg];
                double pt_trg     = pepVect_trg.Pt();
                double eta_trg    = pepVect_trg.Eta();
                double rap_trg    = pepVect_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                //double effxi = effhisto_xi->GetBinContent(effhisto_xi->FindBin(EffXchoice,pt_trg));

                //nMult_trg_eff_xi = nMult_trg_eff_xi + 1.0/effxi;
            }

            mult_xi_bkg[i]->Fill(nMult_trg_eff_xi);

            // SIDEBAND REGION signal
            for(int xi_trg = 0; xi_trg < pepVect_Xi_side_size; xi_trg++)
            {
                TLorentzVector pepVect_trg = (*pepVect_Xi_side[i])[xi_trg];
                double pt_trg     = pepVect_trg.Pt();
                double eta_trg    = pepVect_trg.Eta();
                double phi_trg    = pepVect_trg.Phi();
                double rap_trg    = pepVect_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                //double effxi = effhisto_xi->GetBinContent(effhisto_xi->FindBin(EffXchoice,pt_trg));

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
                    SignalXiSide[i]->Fill(dEta, dPhi);//, 1.0/nMult_trg_eff_xi/effxi);
                }
            }
        }

        for(int i=0 ; i<PtBinNum_; i++)
        {
            PepVect2_Xi_peak[i]->push_back(*pepVect_Xi_peak[i]);
            PepVect2_Xi_side[i]->push_back(*pepVect_Xi_side[i]);
            PepVect2_dau_xi_peak[i]->push_back(*pepVect_dau_xi_peak[i]);
            PepVect2_dau_xi_side[i]->push_back(*pepVect_dau_xi_side[i]);
            delete pepVect_Xi_peak[i];
            delete pepVect_Xi_side[i];
            delete pepVect_dau_xi_peak[i];
            delete pepVect_dau_xi_side[i];
        }
        PepVect2_ass->push_back(*pepVect_trkass);
        zvtxVect->push_back(bestvz);

        /* SignalXiHad
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
        */

        /* Already have these from past jobs so can just use those...
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
        */

        delete pepVect_trkass;

    }
}


// ------------ method called once each job just before starting event loop  ------------
void
XiCorrelation::beginJob()
{
    TH1::SetDefaultSumw2();
    //edm::FileInPath fip("XiAnalyzer/XiAnalyzer/data/EffhistoXi.root");
    //TFile f(fip.fullPath().c_str(),"READ");
    //effhisto_xi = (TH2D*)f.Get("EffHistoXi");

    nTrk            = fs->make<TH1D>("nTrk", "nTrk", 300, 0, 300);
    nEvt            = fs->make<TH1D>("nEvt","nEvt",10,0,10);
    nEvtCut         = fs->make<TH1D>("nEvtCut", "nEvtCut", 10, 0, 10);
    HadPerEvt       = fs->make<TH1D>("HadPerEvent", "Hadrons per Event", 1500, 0, 1500);
    TrkassPerEvt    = fs->make<TH1D>("TrkassPerEvent", "Associated trks per Event", 300,0, 300);
    MassPt          = fs->make<TH2D>("MassPt", "", 150, 1.25, 1.40, 300, 0 ,30);
    //BackgroundHad   = fs->make<TH2D>("BackgroundHad", "BkgHad; #Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    //SignalHad       = fs->make<TH2D>("SignalHad", "SigHad; #Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    //SignalXiHad     = fs->make<TH2D>("SignalXiHad", ";#Delta#eta;#Delta#phi ", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    //BackgroundXiHad = fs->make<TH2D>("BackgroundXiHad", ";#Delta#eta;#Delta#phi ", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    for(int i=0; i<PtBinNum_; i++)
    {
        BackgroundXiPeak[i] = fs->make<TH2D>(Form("BackgroundPeak_pt%d",i), ";#Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        BackgroundXiSide[i] = fs->make<TH2D>(Form("BackgroundSide_pt%d",i), ";#Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        SignalXiPeak[i]     = fs->make<TH2D>(Form("SignalPeak_pt%d",i), ";#Delta#eta;#Delta#phi ", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        SignalXiSide[i]     = fs->make<TH2D>(Form("SignalSide_pt%d",i), ";#Delta#eta;#Delta#phi ", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        KET_xi[i]           = fs->make<TH1D>(Form("KET_xi_pt%d",i),";GeV",40000,0,20);
        KET_xi_bkg[i]       = fs->make<TH1D>(Form("KET_xi_bkg_pt%d",i),";GeV",40000,0,20);
        Mass_xi[i]          = fs->make<TH1D>(Form("Mass_xi_pt%d",i),";GeV",2000,0.8,1.8);
        Pt_xi[i]            = fs->make<TH1D>(Form("Pt_xi_pt%d",i),";GeV",40000,0,20);
        Pt_xi_bkg[i]        = fs->make<TH1D>(Form("Pt_xi_bkg_pt%d",i),";GeV",40000,0,20);
        Eta_xi[i]           = fs->make<TH1D>(Form("Eta_xi_pt%d",i),";#eta",30,-3.0,3.0);
        Eta_xi_bkg[i]       = fs->make<TH1D>(Form("Eta_xi_bkg_pt%d",i),";#eta",30,-3.0,3.0);
        rap_xi[i]           = fs->make<TH1D>(Form("rap_xi_pt%d",i),";y",100,-5,5);
        rap_xi_bkg[i]       = fs->make<TH1D>(Form("rap_xi_bkg_pt%d",i),";y",100,-5,5);
        rap_xi_Lorentz[i]   = fs->make<TH1D>(Form("rap_xi_Lorentz_pt%d",i),";y",100,-5,5);
        mult_xi[i] = fs->make<TH1D>(Form("mult_xi_pt%d",i),"mult_xi",250,0,250);
        mult_xi_bkg[i] = fs->make<TH1D>(Form("mult_xi_bkg_pt%d",i),"mult_xi_bkg",250,0,250);
        PepVect2_Xi_peak[i] = new vector< vector<TLorentzVector> >;
        PepVect2_Xi_side[i] = new vector< vector<TLorentzVector> >;
        PepVect2_dau_xi_peak[i] = new vector<vector<TVector3> >;
        PepVect2_dau_xi_side[i] = new vector<vector<TVector3> >;
    }

    // For Background calculations which must be done in the endJob function
    //
    PepVect2_ass       = new vector< vector<TVector3> >;
    zvtxVect           = new vector<double>;
}

// ------------ method called once each job just after ending the event loop  ------------
void
XiCorrelation::endJob()
{
    int PepVect2_ass_size = (int)PepVect2_ass->size();
    for(int i=0; i<PtBinNum_; i++)
    {
        // Make background histograms
        // Xi paired with hadron
        int PepVect2_Xi_peak_size = (int)PepVect2_Xi_peak[i]->size();
        int PepVect2_Xi_side_size = (int)PepVect2_Xi_side[i]->size();

        // PEAK REGION Background
        for(int bkgnum = 0; bkgnum<bkgnum_; bkgnum++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<PepVect2_Xi_peak_size; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(PepVect2_ass_size);
                if(nevt_trg == nevt_ass)
                {
                    nevt_trg--;
                    continue;
                }
                if(fabs((*zvtxVect)[nevt_trg] - (*zvtxVect)[nevt_ass]) > 0.5)
                {
                    nevt_trg--;
                    ncount++;
                    if(ncount > 5000)
                    {
                        nevt_trg++;
                        ncount=0;
                    }
                    continue;
                }

                vector<TLorentzVector> pepVectTmp_trg = (*PepVect2_Xi_peak[i])[nevt_trg];
                vector<TVector3> pepVectTmp_dau = (*PepVect2_dau_xi_peak[i])[nevt_trg];
                vector<TVector3> pepVectTmp_ass = (*PepVect2_ass)[nevt_ass];
                int nMult_trg = pepVectTmp_trg.size();
                int nMult_ass = pepVectTmp_ass.size();

                double nMult_trg_eff_xi = 0;

                for(int ntrg=0; ntrg<nMult_trg; ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pepVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    //double effxi = effhisto_xi->GetBinContent(effhisto_xi->FindBin(EffXchoice,pt_trg));

                    //nMult_trg_eff_xi = nMult_trg_eff_xi + 1.0/effxi;
                }

                for(int ntrg=0; ntrg<nMult_trg; ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pepVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    //double effxi = effhisto_xi->GetBinContent(effhisto_xi->FindBin(EffXchoice,pt_trg));

                    TVector3 pvectorTmp_dau1_trg = pepVectTmp_dau[2*ntrg];
                    double eta_trg_dau1 = pvectorTmp_dau1_trg.Eta();
                    double phi_trg_dau1 = pvectorTmp_dau1_trg.Phi();
                    TVector3 pvectorTmp_dau2_trg = pepVectTmp_dau[2*ntrg+1];
                    double eta_trg_dau2 = pvectorTmp_dau2_trg.Eta();
                    double phi_trg_dau2 = pvectorTmp_dau2_trg.Phi();


                    for(int nass=0; nass<nMult_ass; nass++)
                    {
                        TVector3 pvectorTmp_ass = pepVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                            if(fabs(eta_ass - eta_trg_dau1)<0.03 && fabs(phi_ass - phi_trg_dau1)<0.03) continue;
                            if(fabs(eta_ass - eta_trg_dau2)<0.03 && fabs(phi_ass - phi_trg_dau2)<0.03) continue;

                        double dEta = eta_ass - eta_trg;
                        double dPhi = phi_ass - phi_trg;

                        if(dPhi > PI)                    dPhi=dPhi-2*PI;
                        if(dPhi < -PI)                   dPhi=dPhi+2*PI;
                        if(dPhi > -PI && dPhi < -PI/2.0) dPhi=dPhi+2*PI;

                        if(fabs(dPhi) < 0.028 && fabs(dEta) < 0.02) continue;

                        BackgroundXiPeak[i]->Fill(dEta, dPhi);//, 1.0/nMult_trg_eff_xi/effxi);
                    }
                }
            }
        }

        // SIDEBAND REGION Background
        for(int bkgnum = 0; bkgnum<bkgnum_; bkgnum++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<PepVect2_Xi_side_size; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(PepVect2_ass_size);
                if(nevt_trg == nevt_ass)
                {
                    nevt_trg--;
                    continue;
                }
                if(fabs((*zvtxVect)[nevt_trg] - (*zvtxVect)[nevt_ass]) > 0.5)
                {
                    nevt_trg--;
                    ncount++;
                    if(ncount > 5000)
                    {
                        nevt_trg++;
                        ncount=0;
                    }
                    continue;
                }

                vector<TLorentzVector> pepVectTmp_trg = (*PepVect2_Xi_side[i])[nevt_trg];
                vector<TVector3> pepVectTmp_dau = (*PepVect2_dau_xi_side[i])[nevt_trg];
                vector<TVector3> pepVectTmp_ass = (*PepVect2_ass)[nevt_ass];
                int nMult_trg = pepVectTmp_trg.size();
                int nMult_ass = pepVectTmp_ass.size();

                double nMult_trg_eff_xi = 0;

                for(int ntrg=0; ntrg<nMult_trg; ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pepVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    //double effxi = effhisto_xi->GetBinContent(effhisto_xi->FindBin(EffXchoice,pt_trg));
                    //nMult_trg_eff_xi = nMult_trg_eff_xi + 1.0/effxi;
                }

                for(int ntrg=0; ntrg<nMult_trg; ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pepVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    //double effxi = effhisto_xi->GetBinContent(effhisto_xi->FindBin(EffXchoice,pt_trg));

                    TVector3 pvectorTmp_dau1_trg = pepVectTmp_dau[2*ntrg];
                    double eta_trg_dau1 = pvectorTmp_dau1_trg.Eta();
                    double phi_trg_dau1 = pvectorTmp_dau1_trg.Phi();
                    TVector3 pvectorTmp_dau2_trg = pepVectTmp_dau[2*ntrg+1];
                    double eta_trg_dau2 = pvectorTmp_dau2_trg.Eta();
                    double phi_trg_dau2 = pvectorTmp_dau2_trg.Phi();

                    for(int nass=0; nass<nMult_ass; nass++)
                    {
                        TVector3 pvectorTmp_ass = pepVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();

                            if(fabs(eta_ass - eta_trg_dau1)<0.03 && fabs(phi_ass - phi_trg_dau1)<0.03) continue;
                            if(fabs(eta_ass - eta_trg_dau2)<0.03 && fabs(phi_ass - phi_trg_dau2)<0.03) continue;

                        double dEta = eta_ass - eta_trg;
                        double dPhi = phi_ass - phi_trg;

                        if(dPhi > PI)                    dPhi=dPhi-2*PI;
                        if(dPhi < -PI)                   dPhi=dPhi+2*PI;
                        if(dPhi > -PI && dPhi < -PI/2.0) dPhi=dPhi+2*PI;

                        if(fabs(dPhi) < 0.028 && fabs(dEta) < 0.02) continue;

                        BackgroundXiSide[i]->Fill(dEta, dPhi);//, 1.0/nMult_trg_eff_xi/effxi);
                    }
                }
            }
        }
    }

    /* BackgroundXiHad
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
    */

    /* BackgroundHad
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
    */

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

