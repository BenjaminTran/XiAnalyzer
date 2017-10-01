//system include files
#include <memory>

#include "../interface/V0CasCorrelation.h"

using namespace std;

#define UNUSED(x) (void)(x)

V0CasCorrelation::V0CasCorrelation(const edm::ParameterSet& iConfig)
{

    //now do what ever initialization is needed
    etaMin_trg_   = iConfig.getUntrackedParameter<double>("etaMin_trg", -2.4);
    etaMax_trg_   = iConfig.getUntrackedParameter<double>("etaMax_trg", 2.4);
    etaMin_ass_   = iConfig.getUntrackedParameter<double>("etaMin_ass", -2.4);
    etaMax_ass_   = iConfig.getUntrackedParameter<double>("etaMax_ass", 2.4);
    ptMin_ass_    = iConfig.getUntrackedParameter<double>("ptMin_ass", 0.3);
    ptMax_ass_    = iConfig.getUntrackedParameter<double>("ptMax_ass", 3.0);
    multMax_      = iConfig.getUntrackedParameter<double>("multHigh", 220);
    multMin_      = iConfig.getUntrackedParameter<double>("multLow", 185);
    bkgnum_       = iConfig.getUntrackedParameter<int>("bkgnum", 20);
    ptbin_n_      = iConfig.getUntrackedParameter<int>("ptbin_n", 13);
    ptbin_n_cas_      = iConfig.getUntrackedParameter<int>("ptbin_n_cas", 13);
    peakFactor_   = iConfig.getUntrackedParameter<double>("peakFactor", 2.0);
    sideFactor_   = iConfig.getUntrackedParameter<double>("sideFactor", 3.0);
    mis_ks_range_ = iConfig.getUntrackedParameter<double>("mis_ks_range", 0.020);
    mis_la_range_ = iConfig.getUntrackedParameter<double>("mis_la_range", 0.010);
    mis_ph_range_ = iConfig.getUntrackedParameter<double>("mis_ph_range", 0.015);


    sigma_ks_ = iConfig.getUntrackedParameter<std::vector<double> >("sigma_ks");
    mean_ks_  = iConfig.getUntrackedParameter<std::vector<double> >("mean_ks");
    sigma_la_ = iConfig.getUntrackedParameter<std::vector<double> >("sigma_la");
    mean_la_  = iConfig.getUntrackedParameter<std::vector<double> >("mean_la");
    sigma_xi_ = iConfig.getUntrackedParameter<std::vector<double> >("sigma_xi");
    mean_xi_  = iConfig.getUntrackedParameter<std::vector<double> >("mean_xi");
    sigma_om_ = iConfig.getUntrackedParameter<std::vector<double> >("sigma_om");
    mean_om_  = iConfig.getUntrackedParameter<std::vector<double> >("mean_om");
    ptcut_ks_ = iConfig.getUntrackedParameter<std::vector<double> >("ptcut_ks");
    ptcut_la_ = iConfig.getUntrackedParameter<std::vector<double> >("ptcut_la");
    ptcut_xi_ = iConfig.getUntrackedParameter<std::vector<double> >("ptcut_xi");
    ptcut_om_ = iConfig.getUntrackedParameter<std::vector<double> >("ptcut_om");

    rejectDaughter_ = iConfig.getUntrackedParameter<bool>("rejectDaughter");
    doRap_          = iConfig.getUntrackedParameter<bool>("doRap");
    doGenRef_       = iConfig.getUntrackedParameter<bool>("doGenRef");

    _vertexCollName = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollName"));
    _trkSrc         = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trkSrc"));
    _ksCollection   = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("ksCollection"));
    _laCollection   = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("laCollection"));
    _xiCollection   = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("xiCollection"));
    _omCollection   = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("omCollection"));
    if(doGenRef_)
    {
        _gnCollection   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gnCollection"));
    }

}


V0CasCorrelation::~V0CasCorrelation()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
V0CasCorrelation::analyze(const edm::Event& iEvent, const edm::EventSetup&
                         iSetup)
{
    using namespace edm;

    // select on requirement of valid vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(_vertexCollName,vertices);
    double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();

    if(bestvz < -15.0 || bestvz>15.0) return;

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_ks;
    iEvent.getByToken(_ksCollection,v0candidates_ks);
    if(!v0candidates_ks.isValid()) return;

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_la;
    iEvent.getByToken(_laCollection,v0candidates_la);
    if(!v0candidates_la.isValid()) return;

    edm::Handle<reco::VertexCompositeCandidateCollection> xiCollection;
    iEvent.getByToken(_xiCollection, xiCollection);
    if(!xiCollection.isValid()) return;

    edm::Handle<reco::VertexCompositeCandidateCollection> omCollection;
    iEvent.getByToken(_omCollection, omCollection);
    if(!omCollection.isValid()) return;

    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(_trkSrc, tracks);

    edm::Handle<reco::GenParticleCollection> genpars;
    if(doGenRef_)
    {
        iEvent.getByToken(_gnCollection,genpars);
    }

    for(int i=0;i<ptbin_n_;i++)
    {
        pVect_trg_ks[i] = new vector<TLorentzVector>;
        pVect_trg_la[i] = new vector<TLorentzVector>;
        pVect_dau_ks[i] = new vector<TVector3>;
        pVect_dau_la[i] = new vector<TVector3>;
        pVect_trg_ks_bkg[i] = new vector<TLorentzVector>;
        pVect_trg_la_bkg[i] = new vector<TLorentzVector>;
        pVect_dau_ks_bkg[i] = new vector<TVector3>;
        pVect_dau_la_bkg[i] = new vector<TVector3>;
    }
    for(int i=0; i<ptbin_n_cas_; i++)
    {
        pepVect_xi_peak[i]    = new vector<TLorentzVector>;
        pepVect_xi_side[i]    = new vector<TLorentzVector>;
        pepVect_dau_xi_peak[i]= new vector<TVector3>;
        pepVect_dau_xi_side[i]= new vector<TVector3>;

        pepVect_om_peak[i]    = new vector<TLorentzVector>;
        pepVect_om_side[i]    = new vector<TLorentzVector>;
        pepVect_dau_om_peak[i]= new vector<TVector3>;
        pepVect_dau_om_side[i]= new vector<TVector3>;
    }

    pVect_ass = new vector<TVector3>;


    //track selection
    int nMult_ass_good = 0;
    for(unsigned it=0; it<tracks->size(); ++it){

        const reco::Track & trk = (*tracks)[it];

        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);

        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt()>0.10) continue;
        if(fabs(dzvtx/dzerror) > 3) continue;
        if(fabs(dxyvtx/dxyerror) > 3) continue;

        double eta = trk.eta();
        double pt  = trk.pt();

        if(fabs(eta)>2.4) continue;
        if(pt<=0.4) continue;
        nMult_ass_good++;
    }
    hMult->Fill(nMult_ass_good);

    if(nMult_ass_good<multMax_ && nMult_ass_good>=multMin_){
        hMult_accept->Fill(nMult_ass_good);
        for(unsigned it=0; it<v0candidates_ks->size(); ++it){

            const reco::VertexCompositeCandidate & trk = (*v0candidates_ks)[it];

            double secvz=-999.9, secvx=-999.9, secvy=-999.9;
            //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;

            secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
            //secvzError = trk.vertexCovariance(2,2); secvxError = trk.vertexCovariance(0,0); secvyError = trk.vertexCovariance(1,1);

            double eta = trk.eta();
            double phi = trk.phi();
            double pt  = trk.pt();
            double mass = trk.mass();
            double px = trk.px();
            double py = trk.py();
            double pz = trk.pz();
            double rapidity = trk.rapidity();
            double EffXchoice = 0;

            if(doRap_)
                EffXchoice = rapidity;
            else
                EffXchoice = eta;

            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(px,py,pz);

            double agl = cos(secvec.Angle(ptosvec));

            if(!doGenRef_)
            {
                if(agl<=0.999) continue;
            }

            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;

            SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);

            double dl = ROOT::Math::Mag(distanceVector);
            double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;

            double dlos = dl/dlerror;

            if(!doGenRef_)
            {
                if(dlos<=5) continue;
            }

            const reco::Candidate * dau1 = trk.daughter(0);
            const reco::Candidate * dau2 = trk.daughter(1);

            double pxd1 = dau1->px();
            double pyd1 = dau1->py();
            double pzd1 = dau1->pz();
            double pd1 = dau1->p();
            double pxd2 = dau2->px();
            double pyd2 = dau2->py();
            double pzd2 = dau2->pz();
            double pd2 = dau2->p();

            TVector3 dauvec1(pxd1,pyd1,pzd1);
            TVector3 dauvec2(pxd2,pyd2,pzd2);

            TVector3 dauvecsum(dauvec1+dauvec2);
            double v0masspiproton1 = sqrt((sqrt(0.93827*0.93827+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))*(sqrt(0.93827*0.93827+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))-dauvecsum.Mag2());

            double v0masspiproton2 = sqrt((sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.93827*0.93827+pd2*pd2))*(sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.93827*0.93827+pd2*pd2))-dauvecsum.Mag2());

            if(!doGenRef_)
            {
                if((v0masspiproton1>=(1.115683-mis_la_range_) && v0masspiproton1<=(1.115683+mis_la_range_)) || (v0masspiproton2>=(1.115683-mis_la_range_) && v0masspiproton2<=(1.115683+mis_la_range_)) ) continue;
            }

            //efficiency
            double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt));

            double eta_dau1 = dau1->eta();
            double phi_dau1 = dau1->phi();
            double pt_dau1 = dau1->pt();

            double eta_dau2 = dau2->eta();
            double phi_dau2 = dau2->phi();
            double pt_dau2 = dau2->pt();

            TLorentzVector pvector;
            pvector.SetPtEtaPhiE(pt,eta,phi,rapidity);

            TVector3 pvector_dau1;
            pvector_dau1.SetPtEtaPhi(pt_dau1,eta_dau1,phi_dau1);

            TVector3 pvector_dau2;
            pvector_dau2.SetPtEtaPhi(pt_dau2,eta_dau2,phi_dau2);

            for(int i=0;i<ptbin_n_;i++)
            {
                double rapOrEtaMaxCut = etaMax_trg_;
                double rapOrEtaMinCut = etaMin_trg_;
                double rapOrEta = eta;

                if(doRap_)
                {
                    rapOrEtaMaxCut = 1.0;
                    rapOrEtaMinCut = -1.0;
                    rapOrEta = rapidity;
                }
                if(rapOrEta<=rapOrEtaMaxCut && rapOrEta>=rapOrEtaMinCut && pt<=ptcut_ks_[i+1] && pt>=ptcut_ks_[i]){
                    hMass_ks[i]->Fill(mass);
                    if(mass<=(mean_ks_[i]+peakFactor_*sigma_ks_[i]) && mass>=(mean_ks_[i]-peakFactor_*sigma_ks_[i])){
                        pVect_trg_ks[i]->push_back(pvector);
                        pVect_dau_ks[i]->push_back(pvector_dau1);
                        pVect_dau_ks[i]->push_back(pvector_dau2);
                        hPt_ks[i]->Fill(pt,1.0/effks);
						hEta_ks[i]->Fill(eta,1.0/effks);
                        hRap_ks[i]->Fill(rapidity);
                        hRap_ks_Lorentz[i]->Fill(pvector.E());
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_ks[i]->Fill(KET,1.0/effks);
                    }
                    if((mass<=(mean_ks_[i]-sideFactor_*sigma_ks_[i]) && mass>=0.425) || (mass<=0.57 && mass>=(mean_ks_[i]+sideFactor_*sigma_ks_[i]))){
                        pVect_trg_ks_bkg[i]->push_back(pvector);
                        pVect_dau_ks_bkg[i]->push_back(pvector_dau1);
                        pVect_dau_ks_bkg[i]->push_back(pvector_dau2);
                        hPt_ks_bkg[i]->Fill(pt,1.0/effks);
						hEta_ks_bkg[i]->Fill(eta,1.0/effks);
                        hRap_ks_bkg[i]->Fill(rapidity);
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_ks_bkg[i]->Fill(KET,1.0/effks);
                    }
                }
            }
        }

        for(unsigned it=0; it<v0candidates_la->size(); ++it){

            const reco::VertexCompositeCandidate & trk = (*v0candidates_la)[it];

            double secvz=-999.9, secvx=-999.9, secvy=-999.9;
            //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;

            secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
            //secvzError = trk.vertexCovariance(2,2); secvxError = trk.vertexCovariance(0,0); secvyError = trk.vertexCovariance(1,1);

            double eta = trk.eta();
            double phi = trk.phi();
            double pt  = trk.pt();
            double mass = trk.mass();
            double px = trk.px();
            double py = trk.py();
            double pz = trk.pz();
            double rapidity = trk.rapidity();
            double EffXchoice = 0;

            if(doRap_)
                EffXchoice = rapidity;
            else
                EffXchoice = eta;

            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(px,py,pz);

            double agl = cos(secvec.Angle(ptosvec));

            if(!doGenRef_)
            {
                if(agl<=0.999) continue;
            }

            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;

            SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);

            double dl = ROOT::Math::Mag(distanceVector);
            double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;

            double dlos = dl/dlerror;

            if(!doGenRef_)
            {
                if(dlos<=5) continue;
            }

            const reco::Candidate * dau1 = trk.daughter(0);
            const reco::Candidate * dau2 = trk.daughter(1);

            double pxd1 = dau1->px();
            double pyd1 = dau1->py();
            double pzd1 = dau1->pz();
            double pd1 = dau1->p();
            double pxd2 = dau2->px();
            double pyd2 = dau2->py();
            double pzd2 = dau2->pz();
            double pd2 = dau2->p();

            TVector3 dauvec1(pxd1,pyd1,pzd1);
            TVector3 dauvec2(pxd2,pyd2,pzd2);

            TVector3 dauvecsum(dauvec1+dauvec2);
            double v0masspipi = sqrt((sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))*(sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))-dauvecsum.Mag2());
            double v0massee = sqrt((sqrt(0.000511*0.000511+pd1*pd1)+sqrt(0.000511*0.000511+pd2*pd2))*(sqrt(0.000511*0.000511+pd1*pd1)+sqrt(0.000511*0.000511+pd2*pd2))-dauvecsum.Mag2());

            if(!doGenRef_)
            {
                if( (v0masspipi>=(0.497614-mis_ks_range_) && v0masspipi<=(0.497614+mis_ks_range_)) || v0massee <= mis_ph_range_ ) continue;
            }

            //efficiency
            double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt));

            double eta_dau1 = dau1->eta();
            double phi_dau1 = dau1->phi();
            double pt_dau1 = dau1->pt();

            double eta_dau2 = dau2->eta();
            double phi_dau2 = dau2->phi();
            double pt_dau2 = dau2->pt();


            TLorentzVector pvector;
            pvector.SetPtEtaPhiE(pt,eta,phi,rapidity);

            TVector3 pvector_dau1;
            pvector_dau1.SetPtEtaPhi(pt_dau1,eta_dau1,phi_dau1);

            TVector3 pvector_dau2;
            pvector_dau2.SetPtEtaPhi(pt_dau2,eta_dau2,phi_dau2);

            for(int i=0;i<ptbin_n_;i++)
            {
                double rapOrEtaMaxCut = etaMax_trg_;
                double rapOrEtaMinCut = etaMin_trg_;
                double rapOrEta = eta;

                if(doRap_)
                {
                    rapOrEtaMaxCut = 1.0;
                    rapOrEtaMinCut = -1.0;
                    rapOrEta = rapidity;
                }
                if(rapOrEta<=rapOrEtaMaxCut && rapOrEta>=rapOrEtaMinCut && pt<=ptcut_la_[i+1] && pt>=ptcut_la_[i]){
                    hMass_la[i]->Fill(mass);
                    if(mass<=(mean_la_[i]+peakFactor_*sigma_la_[i]) && mass>=(mean_la_[i]-peakFactor_*sigma_la_[i])){
                        pVect_trg_la[i]->push_back(pvector);
                        pVect_dau_la[i]->push_back(pvector_dau1);
                        pVect_dau_la[i]->push_back(pvector_dau2);
                        hPt_la[i]->Fill(pt,1.0/effla);
                        hEta_la[i]->Fill(eta,1.0/effla);
                        hRap_la[i]->Fill(rapidity);
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_la[i]->Fill(KET,1.0/effla);
                    }
                    if((mass<=1.165 && mass>=(mean_la_[i]+sideFactor_*sigma_la_[i])) || (mass<=(mean_la_[i]-sideFactor_*sigma_la_[i]) && mass>=1.075)){
                        pVect_trg_la_bkg[i]->push_back(pvector);
                        pVect_dau_la_bkg[i]->push_back(pvector_dau1);
                        pVect_dau_la_bkg[i]->push_back(pvector_dau2);
                        hPt_la_bkg[i]->Fill(pt,1.0/effla);
                        hEta_la_bkg[i]->Fill(eta,1.0/effla);
                        hRap_la_bkg[i]->Fill(rapidity);
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_la_bkg[i]->Fill(KET,1.0/effla);
                    }
                }
            }
        }

        for(reco::VertexCompositeCandidateCollection::const_iterator xiCand =
                xiCollection->begin(); xiCand != xiCollection->end(); xiCand++)
        {
            // Make 2D Mass v Pt
            double xi_eta = xiCand->eta();
            double xi_phi = xiCand->phi();
            double mass   = xiCand->mass();
            double xi_pT  = xiCand->pt();
            double xi_rap = xiCand->rapidity();

            MassPtXi->Fill(mass,xi_pT);
            double Ket = sqrt(mass*mass + xi_pT*xi_pT) - mass;
            double EffXchoice = 0;
            UNUSED(EffXchoice);

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
            for(int i=0; i<ptbin_n_cas_;i++)
            {
                if(xi_pT <= ptcut_xi_[i+1] && xi_pT >= ptcut_xi_[i])
                {
                    Mass_xi[i]->Fill(mass);
                    //peak
                    if(mass >= (mean_xi_[i] - peakFactor_*sigma_xi_[i]) && mass <= (mean_xi_[i] + peakFactor_*sigma_xi_[i]))
                    {
                        pepVect_xi_peak[i]->push_back(xiPEPvector);
                        pepVect_dau_xi_peak[i]->push_back(dau1PEPvector);
                        pepVect_dau_xi_peak[i]->push_back(dau2PEPvector);
                        KET_xi[i]->Fill(Ket);//,1.0/effxi);
                        Pt_xi[i]->Fill(xi_pT);//,1.0/effxi);
                        Eta_xi[i]->Fill(xi_eta);//,1.0/effxi);
                        rap_xi[i]->Fill(xi_rap);//,1.0/effxi);
                    }
                    //sideband
                    if((mass <= (mean_xi_[i] - sideFactor_*sigma_xi_[i]) && mass >= 1.25) || (mass <= 1.40 && mass >= (mean_xi_[i] + sideFactor_*sigma_xi_[i])))
                    {
                        pepVect_xi_side[i]->push_back(xiPEPvector);
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

        for(reco::VertexCompositeCandidateCollection::const_iterator omCand =
                omCollection->begin(); omCand != omCollection->end(); omCand++)
        {
            // Make 2D Mass v Pt
            double om_eta = omCand->eta();
            double om_phi = omCand->phi();
            double mass   = omCand->mass();
            double om_pT  = omCand->pt();
            double om_rap = omCand->rapidity();

            MassPtOm->Fill(mass,om_pT);
            double Ket = sqrt(mass*mass + om_pT*om_pT) - mass;
            double EffXchoice = 0;
            UNUSED(EffXchoice);

            const reco::Candidate *dau1 = omCand->daughter(0);
            const reco::Candidate *dau2 = omCand->daughter(1);

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
                EffXchoice = om_rap;
            else
                EffXchoice = om_eta;

            //efficiency
            //double effom = effhisto_om->GetBinContent(effhisto_om->FindBin(EffXchoice,om_pT));

            // Make vector of Om Candidate parameters
            TLorentzVector omPEPvector;
            omPEPvector.SetPtEtaPhiE(om_pT,om_eta,om_phi,om_rap);
            for(int i=0; i<ptbin_n_cas_;i++)
            {
                cout << "beg for loop" << i << endl;
                if(om_pT <= ptcut_om_[i+1] && om_pT >= ptcut_om_[i])
                {
                    cout << "Inside pt bin" << i << endl;
                    Mass_om[i]->Fill(mass);
                    //peak
                    if(mass >= (mean_om_[i] - peakFactor_*sigma_om_[i]) && mass <= (mean_om_[i] + peakFactor_*sigma_om_[i]))
                    {
                        pepVect_om_peak[i]->push_back(omPEPvector);
                        pepVect_dau_om_peak[i]->push_back(dau1PEPvector);
                        pepVect_dau_om_peak[i]->push_back(dau2PEPvector);
                        KET_om[i]->Fill(Ket);//,1.0/effom);
                        Pt_om[i]->Fill(om_pT);//,1.0/effom);
                        Eta_om[i]->Fill(om_eta);//,1.0/effom);
                        rap_om[i]->Fill(om_rap);//,1.0/effom);
                    }
                    //sideband
                    if((mass <= (mean_om_[i] - sideFactor_*sigma_om_[i]) && mass >= 1.6) || (mass <= 1.75 && mass >= (mean_om_[i] + sideFactor_*sigma_om_[i])))
                    {
                        pepVect_om_side[i]->push_back(omPEPvector);
                        pepVect_dau_om_side[i]->push_back(dau1PEPvector);
                        pepVect_dau_om_side[i]->push_back(dau2PEPvector);
                        KET_om_bkg[i]->Fill(Ket);//,1.0/effom);
                        Pt_om_bkg[i]->Fill(om_pT);//,1.0/effom);
                        Eta_om_bkg[i]->Fill(om_eta);//,1.0/effom);
                        rap_om_bkg[i]->Fill(om_rap);//,1.0/effom);
                    }
                }
            }
        }


        if(!doGenRef_)
        {
            for(unsigned it=0; it<tracks->size(); ++it){

                const reco::Track & trk = (*tracks)[it];

                math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

                double dzvtx = trk.dz(bestvtx);
                double dxyvtx = trk.dxy(bestvtx);
                double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
                double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);

                if(!trk.quality(reco::TrackBase::highPurity)) continue;
                if(fabs(trk.ptError())/trk.pt()>0.10) continue;
                if(fabs(dzvtx/dzerror) > 3) continue;
                if(fabs(dxyvtx/dxyerror) > 3) continue;


                double eta = trk.eta();
                double phi = trk.phi();
                double pt  = trk.pt();

                TVector3 pvector;
                pvector.SetPtEtaPhi(pt,eta,phi);
                if(eta<=etaMax_ass_ && eta>=etaMin_ass_ && pt<=ptMax_ass_ && pt>=ptMin_ass_) pVect_ass->push_back(pvector);
            }
        }
        else
        {
            for(unsigned it=0; it<genpars->size(); ++it){

                const reco::GenParticle & trk = (*genpars)[it];

                double eta = trk.eta();
                double phi = trk.phi();
                double pt  = trk.pt();
                double mass = trk.mass();
                int id = trk.pdgId();
                int st = trk.status();

                TVector3 pvector;
                pvector.SetPtEtaPhi(pt,eta,phi);

                if(trk.eta()<=etaMax_ass_ && trk.eta()>=etaMin_ass_ && trk.pt()<=ptMax_ass_ && trk.pt()>=ptMin_ass_ && fabs(trk.charge())==1 && st==1) pVect_ass->push_back(pvector);
            }
        }

        //Calculating signal
        int nMult_ass = (int)pVect_ass->size();
        hMult_ass->Fill(nMult_ass);

        for(int i=0; i<ptbin_n_; i++)
        {
            int nMult_trg_ks = (int)pVect_trg_ks[i]->size();
            int nMult_trg_la = (int)pVect_trg_la[i]->size();

            double nMult_trg_eff_ks=0;
            double nMult_trg_eff_la=0;

            for(int ntrg=0;ntrg<nMult_trg_ks;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_ks[i])[ntrg];
                double pt_trg = pvector_trg.Pt();
                double eta_trg = pvector_trg.Eta();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

                nMult_trg_eff_ks = nMult_trg_eff_ks + 1.0/effks;
            }

			for(int ntrg=0;ntrg<nMult_trg_la;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_la[i])[ntrg];
                double pt_trg = pvector_trg.Pt();
                double eta_trg = pvector_trg.Eta();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                nMult_trg_eff_la = nMult_trg_eff_la + 1.0/effla;
            }

            hMult_ks[i]->Fill(nMult_trg_ks);
            hMult_la[i]->Fill(nMult_trg_la);

            for(int ntrg=0;ntrg<nMult_trg_ks;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_ks[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

                TVector3 pvector_trg_dau1 = (*pVect_dau_ks[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();

                TVector3 pvector_trg_dau2 = (*pVect_dau_ks[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();

                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    double pt_ass = pvector_ass.Pt();

                    double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

                    if(rejectDaughter_)
                    {
                        if(fabs(eta_ass-eta_trg_dau1)<0.03 && fabs(phi_ass-phi_trg_dau1)<0.03) continue;
                        if(fabs(eta_ass-eta_trg_dau2)<0.03 && fabs(phi_ass-phi_trg_dau2)<0.03) continue;
                    }

                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;

                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_ks[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff_ks/effks/effweight_ass);
                }
            }

            for(int ntrg=0;ntrg<nMult_trg_la;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_la[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                TVector3 pvector_trg_dau1 = (*pVect_dau_la[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();

                TVector3 pvector_trg_dau2 = (*pVect_dau_la[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();

                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    double pt_ass = pvector_ass.Pt();

                    double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

		if(rejectDaughter_)
                {
                    if(fabs(eta_ass-eta_trg_dau1)<0.03 && fabs(phi_ass-phi_trg_dau1)<0.03) continue;
                    if(fabs(eta_ass-eta_trg_dau2)<0.03 && fabs(phi_ass-phi_trg_dau2)<0.03) continue;
                }

                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;

                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_la[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff_la/effla/effweight_ass);
                }
            }

        }

        for(int i=0; i<ptbin_n_; i++)
        {
            int nMult_trg_ks = (int)pVect_trg_ks_bkg[i]->size();
            int nMult_trg_la = (int)pVect_trg_la_bkg[i]->size();

            double nMult_trg_eff_ks=0;
            double nMult_trg_eff_la=0;

            for(int ntrg=0;ntrg<nMult_trg_ks;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_ks_bkg[i])[ntrg];
                double pt_trg = pvector_trg.Pt();
                double eta_trg = pvector_trg.Eta();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

                nMult_trg_eff_ks = nMult_trg_eff_ks + 1.0/effks;
            }

            for(int ntrg=0;ntrg<nMult_trg_la;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_la_bkg[i])[ntrg];
                double pt_trg = pvector_trg.Pt();
                double eta_trg = pvector_trg.Eta();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                nMult_trg_eff_la = nMult_trg_eff_la + 1.0/effla;
            }

            hMult_ks_bkg[i]->Fill(nMult_trg_ks);
            hMult_la_bkg[i]->Fill(nMult_trg_la);

            for(int ntrg=0;ntrg<nMult_trg_ks;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_ks_bkg[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

                TVector3 pvector_trg_dau1 = (*pVect_dau_ks_bkg[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();

                TVector3 pvector_trg_dau2 = (*pVect_dau_ks_bkg[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();

                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    double pt_ass = pvector_ass.Pt();

                    double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

                    if(rejectDaughter_)
                    {
                        if(fabs(eta_ass-eta_trg_dau1)<0.03 && fabs(phi_ass-phi_trg_dau1)<0.03) continue;
                        if(fabs(eta_ass-eta_trg_dau2)<0.03 && fabs(phi_ass-phi_trg_dau2)<0.03) continue;
                    }

                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;

                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_ks_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff_ks/effks/effweight_ass);
                }
            }

            for(int ntrg=0;ntrg<nMult_trg_la;ntrg++)
            {
                TLorentzVector pvector_trg = (*pVect_trg_la_bkg[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                double rap_trg = pvector_trg.E();
                double EffXchoice = 0;

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                TVector3 pvector_trg_dau1 = (*pVect_dau_la_bkg[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();

                TVector3 pvector_trg_dau2 = (*pVect_dau_la_bkg[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();

                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    double pt_ass = pvector_ass.Pt();

                    double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

                    if(rejectDaughter_)
                    {
                        if(fabs(eta_ass-eta_trg_dau1)<0.03 && fabs(phi_ass-phi_trg_dau1)<0.03) continue;
                        if(fabs(eta_ass-eta_trg_dau2)<0.03 && fabs(phi_ass-phi_trg_dau2)<0.03) continue;
                    }

                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;

                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_la_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff_la/effla/effweight_ass);
                }
            }

        }

        for(int i=0; i<ptbin_n_cas_; i++)
        {
            int pepVect_xi_peak_size     = (int)pepVect_xi_peak[i]->size();
            int pepVect_xi_side_size     = (int)pepVect_xi_side[i]->size();

            double nMult_trg_eff_xi = 0;

            for(int xi_ntrg = 0; xi_ntrg<pepVect_xi_peak_size; xi_ntrg++)
            {
                TLorentzVector pepVect_trg = (*pepVect_xi_peak[i])[xi_ntrg];
                double pt_trg     = pepVect_trg.Pt();
                double eta_trg    = pepVect_trg.Eta();
                double rap_trg    = pepVect_trg.E();
                double EffXchoice = 0;
                UNUSED(EffXchoice);

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                //double effxi = effhisto_xi->GetBinContent(effhisto_xi->FindBin(EffXchoice,pt_trg));

                //nMult_trg_eff_xi = nMult_trg_eff_xi + 1.0/effxi;
            }

            mult_xi[i]->Fill(nMult_trg_eff_xi);

            // PEAK REGION signal
            for(int xi_trg = 0; xi_trg < pepVect_xi_peak_size; xi_trg++)
            {
                TLorentzVector pepVect_trg = (*pepVect_xi_peak[i])[xi_trg];
                double pt_trg     = pepVect_trg.Pt();
                double eta_trg    = pepVect_trg.Eta();
                double phi_trg    = pepVect_trg.Phi();
                double rap_trg    = pepVect_trg.E();
                double EffXchoice = 0;
                UNUSED(EffXchoice);

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

                for(int assoc = 0; assoc < nMult_ass; assoc++)
                {
                    TVector3 pepVect_ass = (*pVect_ass)[assoc];
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

            for(int xi_trg = 0; xi_trg<pepVect_xi_side_size; xi_trg++)
            {
                TLorentzVector pepVect_trg = (*pepVect_xi_side[i])[xi_trg];
                double pt_trg     = pepVect_trg.Pt();
                double eta_trg    = pepVect_trg.Eta();
                double rap_trg    = pepVect_trg.E();
                double EffXchoice = 0;
                UNUSED(EffXchoice);

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                //double effxi = effhisto_xi->GetBinContent(effhisto_xi->FindBin(EffXchoice,pt_trg));

                //nMult_trg_eff_xi = nMult_trg_eff_xi + 1.0/effxi;
            }

            mult_xi_bkg[i]->Fill(nMult_trg_eff_xi);

            // SIDEBAND REGION signal
            for(int xi_trg = 0; xi_trg < pepVect_xi_side_size; xi_trg++)
            {
                TLorentzVector pepVect_trg = (*pepVect_xi_side[i])[xi_trg];
                double pt_trg     = pepVect_trg.Pt();
                double eta_trg    = pepVect_trg.Eta();
                double phi_trg    = pepVect_trg.Phi();
                double rap_trg    = pepVect_trg.E();
                double EffXchoice = 0;
                UNUSED(EffXchoice);

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                //double effxi = effhisto_xi->GetBinContent(effhisto_xi->FindBin(EffXchoice,pt_trg));

                for(int assoc = 0; assoc < nMult_ass; assoc++)
                {
                    TVector3 pepVect_ass = (*pVect_ass)[assoc];
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

        for(int i=0; i<ptbin_n_cas_; i++)
        {
            int pepVect_om_peak_size     = (int)pepVect_om_peak[i]->size();
            int pepVect_om_side_size     = (int)pepVect_om_side[i]->size();

            double nMult_trg_eff_om = 0;

            for(int om_ntrg = 0; om_ntrg<pepVect_om_peak_size; om_ntrg++)
            {
                TLorentzVector pepVect_trg = (*pepVect_om_peak[i])[om_ntrg];
                double pt_trg     = pepVect_trg.Pt();
                double eta_trg    = pepVect_trg.Eta();
                double rap_trg    = pepVect_trg.E();
                double EffXchoice = 0;
                UNUSED(EffXchoice);

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                //double effom = effhisto_om->GetBinContent(effhisto_om->FindBin(EffXchoice,pt_trg));

                //nMult_trg_eff_om = nMult_trg_eff_om + 1.0/effom;
            }

            mult_om[i]->Fill(nMult_trg_eff_om);

            // PEAK REGION signal
            for(int om_trg = 0; om_trg < pepVect_om_peak_size; om_trg++)
            {
                TLorentzVector pepVect_trg = (*pepVect_om_peak[i])[om_trg];
                double pt_trg     = pepVect_trg.Pt();
                double eta_trg    = pepVect_trg.Eta();
                double phi_trg    = pepVect_trg.Phi();
                double rap_trg    = pepVect_trg.E();
                double EffXchoice = 0;
                UNUSED(EffXchoice);

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                //double effom = effhisto_om->GetBinContent(effhisto_om->FindBin(EffXchoice,pt_trg));

                TVector3 pepvector_trg_dau1 = (*pepVect_dau_om_peak[i])[2*om_trg];

                double eta_trg_dau1 = pepvector_trg_dau1.Eta();
                double phi_trg_dau1 = pepvector_trg_dau1.Phi();

                TVector3 pepvector_trg_dau2 = (*pepVect_dau_om_peak[i])[2*om_trg+1];

                double eta_trg_dau2 = pepvector_trg_dau2.Eta();
                double phi_trg_dau2 = pepvector_trg_dau2.Phi();

                for(int assoc = 0; assoc < nMult_ass; assoc++)
                {
                    TVector3 pepVect_ass = (*pVect_ass)[assoc];
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
                    SignalOmPeak[i]->Fill(dEta, dPhi);//,1.0/nMult_trg_eff_om/effom);
                }
            }

            nMult_trg_eff_om = 0;

            for(int om_trg = 0; om_trg<pepVect_om_side_size; om_trg++)
            {
                TLorentzVector pepVect_trg = (*pepVect_om_side[i])[om_trg];
                double pt_trg     = pepVect_trg.Pt();
                double eta_trg    = pepVect_trg.Eta();
                double rap_trg    = pepVect_trg.E();
                double EffXchoice = 0;
                UNUSED(EffXchoice);

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                //double effom = effhisto_om->GetBinContent(effhisto_om->FindBin(EffXchoice,pt_trg));

                //nMult_trg_eff_om = nMult_trg_eff_om + 1.0/effom;
            }

            mult_om_bkg[i]->Fill(nMult_trg_eff_om);

            // SIDEBAND REGION signal
            for(int om_trg = 0; om_trg < pepVect_om_side_size; om_trg++)
            {
                TLorentzVector pepVect_trg = (*pepVect_om_side[i])[om_trg];
                double pt_trg     = pepVect_trg.Pt();
                double eta_trg    = pepVect_trg.Eta();
                double phi_trg    = pepVect_trg.Phi();
                double rap_trg    = pepVect_trg.E();
                double EffXchoice = 0;
                UNUSED(EffXchoice);

                if(doRap_)
                    EffXchoice = rap_trg;
                else
                    EffXchoice = eta_trg;

                //double effom = effhisto_om->GetBinContent(effhisto_om->FindBin(EffXchoice,pt_trg));

                for(int assoc = 0; assoc < nMult_ass; assoc++)
                {
                    TVector3 pepVect_ass = (*pVect_ass)[assoc];
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
                    SignalOmSide[i]->Fill(dEta, dPhi);//, 1.0/nMult_trg_eff_om/effom);
                }
            }
        }


        for(int i=0; i<ptbin_n_; i++)
        {
            pVectVect_trg_ks[i]->push_back(*pVect_trg_ks[i]);
            pVectVect_trg_la[i]->push_back(*pVect_trg_la[i]);
            pVectVect_dau_ks[i]->push_back(*pVect_dau_ks[i]);
            pVectVect_dau_la[i]->push_back(*pVect_dau_la[i]);
            delete pVect_trg_ks[i];
            delete pVect_trg_la[i];
            delete pVect_dau_ks[i];
            delete pVect_dau_la[i];
            pVectVect_trg_ks_bkg[i]->push_back(*pVect_trg_ks_bkg[i]);
            pVectVect_trg_la_bkg[i]->push_back(*pVect_trg_la_bkg[i]);
            pVectVect_dau_ks_bkg[i]->push_back(*pVect_dau_ks_bkg[i]);
            pVectVect_dau_la_bkg[i]->push_back(*pVect_dau_la_bkg[i]);
            delete pVect_trg_ks_bkg[i];
            delete pVect_trg_la_bkg[i];
            delete pVect_dau_ks_bkg[i];
            delete pVect_dau_la_bkg[i];
        }
        for(int i=0 ; i<ptbin_n_cas_; i++)
        {
            PepVect2_xi_peak[i]->push_back(*pepVect_xi_peak[i]);
            PepVect2_xi_side[i]->push_back(*pepVect_xi_side[i]);
            PepVect2_dau_xi_peak[i]->push_back(*pepVect_dau_xi_peak[i]);
            PepVect2_dau_xi_side[i]->push_back(*pepVect_dau_xi_side[i]);
            delete pepVect_xi_peak[i];
            delete pepVect_xi_side[i];
            delete pepVect_dau_xi_peak[i];
            delete pepVect_dau_xi_side[i];

            PepVect2_om_peak[i]->push_back(*pepVect_om_peak[i]);
            PepVect2_om_side[i]->push_back(*pepVect_om_side[i]);
            PepVect2_dau_om_peak[i]->push_back(*pepVect_dau_om_peak[i]);
            PepVect2_dau_om_side[i]->push_back(*pepVect_dau_om_side[i]);
            delete pepVect_om_peak[i];
            delete pepVect_om_side[i];
            delete pepVect_dau_om_peak[i];
            delete pepVect_dau_om_side[i];
        }

        pVectVect_ass->push_back(*pVect_ass);
        zvtxVect->push_back(bestvz);
        delete pVect_ass;
    }
}


// ------------ method called once each job just before starting event
//loop  ------------
void
V0CasCorrelation::beginJob()
{
    edm::Service<TFileService> fs;

    TH1::SetDefaultSumw2();

    edm::FileInPath fip2("XiAnalyzer/XiAnalyzer/data/Hijing_8TeV_dataBS.root");
    TFile f2(fip2.fullPath().c_str(),"READ");
    effhisto = (TH2F*)f2.Get("rTotalEff3D_0");

	/*
    edm::FileInPath fip1("Demo/DemoAnalyzer/data/V0_pp13TeV_Efficiency.root");
    TFile f1(fip1.fullPath().c_str(),"READ");
    effhisto_ks = (TH2D*)f1.Get("ks_eff_1");
    effhisto_la = (TH2D*)f1.Get("la_eff_1");
    */

    //When go back to eta cut then need to add if statement on doRap_
    edm::FileInPath fip("XiAnalyzer/XiAnalyzer/data/Effhisto.root");
    TFile f(fip.fullPath().c_str(),"READ");
    effhisto_ks = (TH2D*)f.Get("EffHistoKs");
    effhisto_la = (TH2D*)f.Get("EffHistoLa");

    MassPtXi          = fs->make<TH2D>("MassPtXi", "", 750, 1.25, 2.00, 300, 0 ,30);
    MassPtOm          = fs->make<TH2D>("MassPtOm", "", 750, 1.25, 2.00, 300, 0 ,30);
    hMult = fs->make<TH1D>("mult",";N",600,0,600);
    hMult_ass = fs->make<TH1D>("mult_ass",";N",600,0,600);
    hMult_accept = fs->make<TH1D>("mult_acc",";N",600,0,600);

    for(int i=0; i<ptbin_n_; i++)
    {
        hKET_ks[i] = fs->make<TH1D>(Form("KETkshort_pt%d",i),";GeV",50000,0,25);
        hKET_la[i] = fs->make<TH1D>(Form("KETlambda_pt%d",i),";GeV",50000,0,25);
        hPt_ks[i] = fs->make<TH1D>(Form("Ptkshort_pt%d",i),";GeV",50000,0,25);
        hPt_la[i] = fs->make<TH1D>(Form("Ptlambda_pt%d",i),";GeV",50000,0,25);
        hEta_ks[i] = fs->make<TH1D>(Form("Etakshort_pt%d",i),";eta",24,-2.4,2.4);
        hRap_ks[i] = fs->make<TH1D>(Form("Rapkshort_pt%d",i),";y",30,-1.5,1.5);
        hRap_ks_Lorentz[i] = fs->make<TH1D>(Form("RapksLorentz_pt%d",i),";y",30,-1.5,1.5);
        hEta_la[i] = fs->make<TH1D>(Form("Etalambda_pt%d",i),";eta",24,-2.4,2.4);
        hRap_la[i] = fs->make<TH1D>(Form("Raplambda_pt%d",i),";y",30,-1.5,1.5);
        hMass_ks[i] = fs->make<TH1D>(Form("masskshort_pt%d",i),";GeV",2000,0,1.0);
        hMass_la[i] = fs->make<TH1D>(Form("masslambda_pt%d",i),";GeV",2000,0.5,1.5);
        hMult_ks[i] = fs->make<TH1D>(Form("mult_ks_pt%d",i),";N",250,0,250);
        hSignal_ks[i] = fs->make<TH2D>(Form("signalkshort_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_ks[i] = fs->make<TH2D>(Form("backgroundkshort_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hMult_la[i] = fs->make<TH1D>(Form("mult_la_pt%d",i),";N",250,0,250);
        hSignal_la[i] = fs->make<TH2D>(Form("signallambda_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_la[i] = fs->make<TH2D>(Form("backgroundlambda_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        pVectVect_trg_ks[i] = new vector< vector<TLorentzVector> >;
        pVectVect_trg_la[i] = new vector< vector<TLorentzVector> >;
        pVectVect_dau_ks[i] = new vector< vector<TVector3> >;
        pVectVect_dau_la[i] = new vector< vector<TVector3> >;
        hKET_ks_bkg[i] = fs->make<TH1D>(Form("KETkshort_bkg_pt%d",i),";GeV",50000,0,25);
        hKET_la_bkg[i] = fs->make<TH1D>(Form("KETlambda_bkg_pt%d",i),";GeV",50000,0,25);
        hPt_ks_bkg[i] = fs->make<TH1D>(Form("Ptkshort_bkg_pt%d",i),";GeV",50000,0,25);
        hPt_la_bkg[i] = fs->make<TH1D>(Form("Ptlambda_bkg_pt%d",i),";GeV",50000,0,25);
        hEta_ks_bkg[i] = fs->make<TH1D>(Form("Etakshort_bkg_pt%d",i),";GeV",24,-2.4,2.4);
        hRap_ks_bkg[i] = fs->make<TH1D>(Form("Rapkshort_bkg_pt%d",i),";y",30,-1.5,1.5);
        hEta_la_bkg[i] = fs->make<TH1D>(Form("Etalambda_bkg_pt%d",i),";GeV",24,-2.4,2.4);
        hRap_la_bkg[i] = fs->make<TH1D>(Form("Raplambda_bkg_pt%d",i),";y",30,-1.5,1.5);
        hMult_ks_bkg[i] = fs->make<TH1D>(Form("mult_ks_bkg_pt%d",i),";N",250,0,250);
        hSignal_ks_bkg[i] = fs->make<TH2D>(Form("signalkshort_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_ks_bkg[i] = fs->make<TH2D>(Form("backgroundkshort_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hMult_la_bkg[i] = fs->make<TH1D>(Form("mult_la_bkg_pt%d",i),";N",250,0,250);
        hSignal_la_bkg[i] = fs->make<TH2D>(Form("signallambda_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_la_bkg[i] = fs->make<TH2D>(Form("backgroundlambda_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        pVectVect_trg_ks_bkg[i] = new vector< vector<TLorentzVector> >;
        pVectVect_trg_la_bkg[i] = new vector< vector<TLorentzVector> >;
        pVectVect_dau_ks_bkg[i] = new vector< vector<TVector3> >;
        pVectVect_dau_la_bkg[i] = new vector< vector<TVector3> >;
    }

    for(int i=0; i<ptbin_n_cas_; i++)
    {
        BackgroundXiPeak[i] = fs->make<TH2D>(Form("BackgroundXiPeak_pt%d",i), ";#Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        BackgroundOmPeak[i] = fs->make<TH2D>(Form("BackgroundOmPeak_pt%d",i), ";#Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        BackgroundXiSide[i] = fs->make<TH2D>(Form("BackgroundXiSide_pt%d",i), ";#Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        BackgroundOmSide[i] = fs->make<TH2D>(Form("BackgroundOmSide_pt%d",i), ";#Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        SignalXiPeak[i]     = fs->make<TH2D>(Form("SignalXiPeak_pt%d",i), ";#Delta#eta;#Delta#phi ", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        SignalOmPeak[i]     = fs->make<TH2D>(Form("SignalOmPeak_pt%d",i), ";#Delta#eta;#Delta#phi ", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        SignalXiSide[i]     = fs->make<TH2D>(Form("SignalXiSide_pt%d",i), ";#Delta#eta;#Delta#phi ", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        SignalOmSide[i]     = fs->make<TH2D>(Form("SignalOmSide_pt%d",i), ";#Delta#eta;#Delta#phi ", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        KET_xi[i]           = fs->make<TH1D>(Form("KET_xi_pt%d",i),";GeV",40000,0,20);
        KET_om[i]           = fs->make<TH1D>(Form("KET_om_pt%d",i),";GeV",40000,0,20);
        KET_xi_bkg[i]       = fs->make<TH1D>(Form("KET_xi_bkg_pt%d",i),";GeV",40000,0,20);
        KET_om_bkg[i]       = fs->make<TH1D>(Form("KET_om_bkg_pt%d",i),";GeV",40000,0,20);
        Mass_xi[i]          = fs->make<TH1D>(Form("Mass_xi_pt%d",i),";GeV",2000,0.8,2.0);
        Mass_om[i]          = fs->make<TH1D>(Form("Mass_om_pt%d",i),";GeV",2000,0.8,2.0);
        Pt_xi[i]            = fs->make<TH1D>(Form("Pt_xi_pt%d",i),";GeV",40000,0,20);
        Pt_om[i]            = fs->make<TH1D>(Form("Pt_om_pt%d",i),";GeV",40000,0,20);
        Pt_xi_bkg[i]        = fs->make<TH1D>(Form("Pt_xi_bkg_pt%d",i),";GeV",40000,0,20);
        Pt_om_bkg[i]        = fs->make<TH1D>(Form("Pt_om_bkg_pt%d",i),";GeV",40000,0,20);
        Eta_xi[i]           = fs->make<TH1D>(Form("Eta_xi_pt%d",i),";#eta",30,-3.0,3.0);
        Eta_om[i]           = fs->make<TH1D>(Form("Eta_om_pt%d",i),";#eta",30,-3.0,3.0);
        Eta_xi_bkg[i]       = fs->make<TH1D>(Form("Eta_xi_bkg_pt%d",i),";#eta",30,-3.0,3.0);
        Eta_om_bkg[i]       = fs->make<TH1D>(Form("Eta_om_bkg_pt%d",i),";#eta",30,-3.0,3.0);
        rap_xi[i]           = fs->make<TH1D>(Form("rap_xi_pt%d",i),";y",100,-5,5);
        rap_om[i]           = fs->make<TH1D>(Form("rap_om_pt%d",i),";y",100,-5,5);
        rap_xi_bkg[i]       = fs->make<TH1D>(Form("rap_xi_bkg_pt%d",i),";y",100,-5,5);
        rap_om_bkg[i]       = fs->make<TH1D>(Form("rap_om_bkg_pt%d",i),";y",100,-5,5);
        mult_xi[i] = fs->make<TH1D>(Form("mult_xi_pt%d",i),"mult_xi",250,0,250);
        mult_om[i] = fs->make<TH1D>(Form("mult_om_pt%d",i),"mult_om",250,0,250);
        mult_xi_bkg[i] = fs->make<TH1D>(Form("mult_xi_bkg_pt%d",i),"mult_xi_bkg",250,0,250);
        mult_om_bkg[i] = fs->make<TH1D>(Form("mult_om_bkg_pt%d",i),"mult_om_bkg",250,0,250);
        PepVect2_xi_peak[i] = new vector< vector<TLorentzVector> >;
        PepVect2_om_peak[i] = new vector< vector<TLorentzVector> >;
        PepVect2_xi_side[i] = new vector< vector<TLorentzVector> >;
        PepVect2_om_side[i] = new vector< vector<TLorentzVector> >;
        PepVect2_dau_xi_peak[i] = new vector<vector<TVector3> >;
        PepVect2_dau_om_peak[i] = new vector<vector<TVector3> >;
        PepVect2_dau_xi_side[i] = new vector<vector<TVector3> >;
        PepVect2_dau_om_side[i] = new vector<vector<TVector3> >;
    }

    pVectVect_ass = new vector< vector<TVector3> >;
    zvtxVect = new vector<double>;

}

// ------------ method called once each job just after ending the event
//loop  ------------
void
V0CasCorrelation::endJob() {
    //Calculating background
    int nevttotal_ass = (int)pVectVect_ass->size();

    /*double etacut[7] = {-2.4,-1.6,-0.8,0,0.8,1.6,2.4};
    double kseff[6] = {1.2,1.2,1.2,1.2,1.2,1.2};
    double laeff[6] = {2,2,2,2,2,1.8};*/

    for(int i=0;i<ptbin_n_;i++)
    {
        int nevttotal_trg_ks = (int)pVectVect_trg_ks[i]->size();
        int nevttotal_trg_la = (int)pVectVect_trg_la[i]->size();

        for(int nround=0;nround<bkgnum_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_ks; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }

                vector<TLorentzVector> pVectTmp_trg = (*pVectVect_trg_ks[i])[nevt_trg];
                vector<TVector3> pVectTmp_dau = (*pVectVect_dau_ks[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();

                double nMult_trg_eff=0;

				for(int ntrg=0;ntrg<nMult_trg;ntrg++)
				{
						TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
						double pt_trg = pvectorTmp_trg.Pt();
						double eta_trg = pvectorTmp_trg.Eta();
                        double rap_trg = pvectorTmp_trg.E();
                        double EffXchoice = 0;

                        if(doRap_)
                            EffXchoice = rap_trg;
                        else
                            EffXchoice = eta_trg;

						double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

						nMult_trg_eff = nMult_trg_eff + 1.0/effks;
				}

                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

                    TVector3 pvector_trg_dau1 = pVectTmp_dau[2*ntrg];
                    double eta_trg_dau1 = pvector_trg_dau1.Eta();
                    double phi_trg_dau1 = pvector_trg_dau1.Phi();

                    TVector3 pvector_trg_dau2 = pVectTmp_dau[2*ntrg+1];
                    double eta_trg_dau2 = pvector_trg_dau2.Eta();
                    double phi_trg_dau2 = pvector_trg_dau2.Phi();

                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        double pt_ass = pvectorTmp_ass.Pt();

                        double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

                        if(rejectDaughter_)
                        {
                            if(fabs(eta_ass - eta_trg_dau1)<0.03 && fabs(phi_ass - phi_trg_dau1)<0.03) continue;
                            if(fabs(eta_ass - eta_trg_dau2)<0.03 && fabs(phi_ass - phi_trg_dau2)<0.03) continue;
                        }

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;

                        //if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_ks[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effks/effweight_ass);
                    }
                }
            }
        }

        for(int nround=0;nround<bkgnum_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_la; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }

                vector<TLorentzVector> pVectTmp_trg = (*pVectVect_trg_la[i])[nevt_trg];
                vector<TVector3> pVectTmp_dau = (*pVectVect_dau_la[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();

                double nMult_trg_eff=0;

                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                    nMult_trg_eff = nMult_trg_eff + 1.0/effla;
                }

                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                    TVector3 pvector_trg_dau1 = pVectTmp_dau[2*ntrg];
                    double eta_trg_dau1 = pvector_trg_dau1.Eta();
                    double phi_trg_dau1 = pvector_trg_dau1.Phi();

                    TVector3 pvector_trg_dau2 = pVectTmp_dau[2*ntrg+1];
                    double eta_trg_dau2 = pvector_trg_dau2.Eta();
                    double phi_trg_dau2 = pvector_trg_dau2.Phi();

                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        double pt_ass = pvectorTmp_ass.Pt();

                        double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

                        if(rejectDaughter_)
                        {
                            if(fabs(eta_ass - eta_trg_dau1)<0.03 && fabs(phi_ass - phi_trg_dau1)<0.03) continue;
                            if(fabs(eta_ass - eta_trg_dau2)<0.03 && fabs(phi_ass - phi_trg_dau2)<0.03) continue;
                        }

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;

                        //if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_la[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effla/effweight_ass);
                    }
                }
            }
        }

    }

    for(int i=0;i<ptbin_n_;i++)
    {
        int nevttotal_trg_ks = (int)pVectVect_trg_ks_bkg[i]->size();
        int nevttotal_trg_la = (int)pVectVect_trg_la_bkg[i]->size();

        for(int nround=0;nround<bkgnum_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_ks; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }

                vector<TLorentzVector> pVectTmp_trg = (*pVectVect_trg_ks_bkg[i])[nevt_trg];
                vector<TVector3> pVectTmp_dau = (*pVectVect_dau_ks_bkg[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();

                double nMult_trg_eff=0;

                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));
                    nMult_trg_eff = nMult_trg_eff + 1.0/effks;
                }

                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    double effks = effhisto_ks->GetBinContent(effhisto_ks->FindBin(EffXchoice,pt_trg));

                    TVector3 pvector_trg_dau1 = pVectTmp_dau[2*ntrg];
                    double eta_trg_dau1 = pvector_trg_dau1.Eta();
                    double phi_trg_dau1 = pvector_trg_dau1.Phi();

                    TVector3 pvector_trg_dau2 = pVectTmp_dau[2*ntrg+1];
                    double eta_trg_dau2 = pvector_trg_dau2.Eta();
                    double phi_trg_dau2 = pvector_trg_dau2.Phi();

                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        double pt_ass = pvectorTmp_ass.Pt();

                        double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

                        if(rejectDaughter_)
                        {
                            if(fabs(eta_ass - eta_trg_dau1)<0.03 && fabs(phi_ass - phi_trg_dau1)<0.03) continue;
                            if(fabs(eta_ass - eta_trg_dau2)<0.03 && fabs(phi_ass - phi_trg_dau2)<0.03) continue;
                        }

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;

                        //if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_ks_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effks/effweight_ass);
                    }
                }
            }
        }

        for(int nround=0;nround<bkgnum_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_la; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }

                vector<TLorentzVector> pVectTmp_trg = (*pVectVect_trg_la_bkg[i])[nevt_trg];
                vector<TVector3> pVectTmp_dau = (*pVectVect_dau_la_bkg[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();

                double nMult_trg_eff=0;

                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                    nMult_trg_eff = nMult_trg_eff + 1.0/effla;
                }

                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    double effla = effhisto_la->GetBinContent(effhisto_la->FindBin(EffXchoice,pt_trg));

                    TVector3 pvector_trg_dau1 = pVectTmp_dau[2*ntrg];
                    double eta_trg_dau1 = pvector_trg_dau1.Eta();
                    double phi_trg_dau1 = pvector_trg_dau1.Phi();

                    TVector3 pvector_trg_dau2 = pVectTmp_dau[2*ntrg+1];
                    double eta_trg_dau2 = pvector_trg_dau2.Eta();
                    double phi_trg_dau2 = pvector_trg_dau2.Phi();

                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        double pt_ass = pvectorTmp_ass.Pt();

                        double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

                        if(rejectDaughter_)
                        {
                            if(fabs(eta_ass - eta_trg_dau1)<0.03 && fabs(phi_ass - phi_trg_dau1)<0.03) continue;
                            if(fabs(eta_ass - eta_trg_dau2)<0.03 && fabs(phi_ass - phi_trg_dau2)<0.03) continue;
                        }

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;

                        //if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_la_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effla/effweight_ass);
                    }
                }
            }
        }

    }

    for(int i=0; i<ptbin_n_cas_; i++)
    {
        // Make background histograms
        // Xi paired with hadron
        int PepVect2_xi_peak_size = (int)PepVect2_xi_peak[i]->size();
        int PepVect2_xi_side_size = (int)PepVect2_xi_side[i]->size();

        // PEAK REGION Background
        for(int bkgnum = 0; bkgnum<bkgnum_; bkgnum++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<PepVect2_xi_peak_size; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
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

                vector<TLorentzVector> pepVectTmp_trg = (*PepVect2_xi_peak[i])[nevt_trg];
                vector<TVector3> pepVectTmp_dau = (*PepVect2_dau_xi_peak[i])[nevt_trg];
                vector<TVector3> pepVectTmp_ass = (*pVectVect_ass)[nevt_ass];
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
                    UNUSED(EffXchoice);

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
                    UNUSED(EffXchoice);

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
            for(int nevt_trg=0; nevt_trg<PepVect2_xi_side_size; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
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

                vector<TLorentzVector> pepVectTmp_trg = (*PepVect2_xi_side[i])[nevt_trg];
                vector<TVector3> pepVectTmp_dau = (*PepVect2_dau_xi_side[i])[nevt_trg];
                vector<TVector3> pepVectTmp_ass = (*pVectVect_ass)[nevt_ass];
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
                    UNUSED(EffXchoice);

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
                    UNUSED(EffXchoice);

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

    for(int i=0; i<ptbin_n_cas_; i++)
    {
        // Make background histograms
        // Xi paired with hadron
        int PepVect2_om_peak_size = (int)PepVect2_om_peak[i]->size();
        int PepVect2_om_side_size = (int)PepVect2_om_side[i]->size();

        // PEAK REGION Background
        for(int bkgnum = 0; bkgnum<bkgnum_; bkgnum++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<PepVect2_om_peak_size; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
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

                vector<TLorentzVector> pepVectTmp_trg = (*PepVect2_om_peak[i])[nevt_trg];
                vector<TVector3> pepVectTmp_dau = (*PepVect2_dau_om_peak[i])[nevt_trg];
                vector<TVector3> pepVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pepVectTmp_trg.size();
                int nMult_ass = pepVectTmp_ass.size();

                double nMult_trg_eff_om = 0;

                for(int ntrg=0; ntrg<nMult_trg; ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pepVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;
                    UNUSED(EffXchoice);

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    //double effom = effhisto_om->GetBinContent(effhisto_om->FindBin(EffXchoice,pt_trg));

                    //nMult_trg_eff_om = nMult_trg_eff_om + 1.0/effom;
                }

                for(int ntrg=0; ntrg<nMult_trg; ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pepVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;
                    UNUSED(EffXchoice);

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    //double effom = effhisto_om->GetBinContent(effhisto_om->FindBin(EffXchoice,pt_trg));

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

                        BackgroundOmPeak[i]->Fill(dEta, dPhi);//, 1.0/nMult_trg_eff_om/effom);
                    }
                }
            }
        }

        // SIDEBAND REGION Background
        for(int bkgnum = 0; bkgnum<bkgnum_; bkgnum++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<PepVect2_om_side_size; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
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

                vector<TLorentzVector> pepVectTmp_trg = (*PepVect2_om_side[i])[nevt_trg];
                vector<TVector3> pepVectTmp_dau = (*PepVect2_dau_om_side[i])[nevt_trg];
                vector<TVector3> pepVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pepVectTmp_trg.size();
                int nMult_ass = pepVectTmp_ass.size();

                double nMult_trg_eff_om = 0;

                for(int ntrg=0; ntrg<nMult_trg; ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pepVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;
                    UNUSED(EffXchoice);

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    //double effom = effhisto_om->GetBinContent(effhisto_om->FindBin(EffXchoice,pt_trg));
                    //nMult_trg_eff_om = nMult_trg_eff_om + 1.0/effom;
                }

                for(int ntrg=0; ntrg<nMult_trg; ntrg++)
                {
                    TLorentzVector pvectorTmp_trg = pepVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double rap_trg = pvectorTmp_trg.E();
                    double EffXchoice = 0;
                    UNUSED(EffXchoice);

                    if(doRap_)
                        EffXchoice = rap_trg;
                    else
                        EffXchoice = eta_trg;

                    //double effom = effhisto_om->GetBinContent(effhisto_om->FindBin(EffXchoice,pt_trg));

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

                        BackgroundOmSide[i]->Fill(dEta, dPhi);//, 1.0/nMult_trg_eff_om/effom);
                    }
                }
            }
        }
    }
}
