#include "XiAnalyzer/XiAnalyzer/interface/V0CorrelationMC.h"

using namespace std;

V0CorrelationMC::V0CorrelationMC(const edm::ParameterSet& iConfig)
{

    //now do what ever initialization is needed
    _gnCollection = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gnCollection"));
    _trkSrc         = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trkSrc"));
    _vertexCollName = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollName"));
    etaMin_trg_   = iConfig.getUntrackedParameter<double>("etaMin_trg", -2.4);
    etaMax_trg_   = iConfig.getUntrackedParameter<double>("etaMax_trg", 2.4);
    etaMin_ass_   = iConfig.getUntrackedParameter<double>("etaMin_ass", -2.4);
    etaMax_ass_   = iConfig.getUntrackedParameter<double>("etaMax_ass", 2.4);
    ptMin_ass_    = iConfig.getUntrackedParameter<double>("ptMin_ass", 0.3);
    ptMax_ass_    = iConfig.getUntrackedParameter<double>("ptMax_ass", 3.0);
    multMax_      = iConfig.getUntrackedParameter<double>("multHigh", 250);
    multMin_      = iConfig.getUntrackedParameter<double>("multLow", 185);
    bkgFactor_    = iConfig.getUntrackedParameter<int>("bkgnum", 10);
    rapMin_       = iConfig.getUntrackedParameter<double>("rapMin",0);
    rapMax_       = iConfig.getUntrackedParameter<double>("rapMax",0);
    doRap_        = iConfig.getUntrackedParameter<bool>("doRap",false);
    ptcut_ks_ = iConfig.getUntrackedParameter<std::vector<double> >("ptcut_ks");
    ptcut_la_ = iConfig.getUntrackedParameter<std::vector<double> >("ptcut_la");
    ptcut_xi_ = iConfig.getUntrackedParameter<std::vector<double> >("ptcut_xi");
}


V0CorrelationMC::~V0CorrelationMC()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
V0CorrelationMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    edm::Handle<reco::GenParticleCollection> genpars;
    iEvent.getByToken(_gnCollection,genpars);

    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(_trkSrc, tracks);

    for(int i=0;i<18;i++)
    {
        pVect_trg_ks[i] = new vector<TVector3>;
        pVect_trg_la[i] = new vector<TVector3>;
        pVect_trg_xi[i] = new vector<TVector3>;
        pVect_dau_ks[i] = new vector<TVector3>;
        pVect_dau_la[i] = new vector<TVector3>;
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
        //loop over tracks
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

            int nm;

            if(trk.numberOfMothers()==1){
                const reco::Candidate * mom = trk.mother();
                if(mom->numberOfMothers()!=0){
                    nm = 1;
                }
            }

            if(fabs(id)==310){
                double eta_dau1 = 0;
                double phi_dau1 = 0;
                double pt_dau1 = 999.999;

                double eta_dau2 = 0;
                double phi_dau2 = 0;
                double pt_dau2 = 999.999;

                if(trk.numberOfDaughters()==2){
                    const reco::Candidate * dau1 = trk.daughter(0);
                    const reco::Candidate * dau2 = trk.daughter(1);

                    eta_dau1 = dau1->eta();
                    phi_dau1 = dau1->phi();
                    pt_dau1 = dau1->pt();

                    eta_dau2 = dau2->eta();
                    phi_dau2 = dau2->phi();
                    pt_dau2 = dau2->pt();
                }

                TVector3 pvector;
                pvector.SetPtEtaPhi(pt,eta,phi);

                TVector3 pvector_dau1;
                pvector_dau1.SetPtEtaPhi(pt_dau1,eta_dau1,phi_dau1);

                TVector3 pvector_dau2;
                pvector_dau2.SetPtEtaPhi(pt_dau2,eta_dau2,phi_dau2);

                for(int i=0;i<18;i++)
                {
                    if(doRap_)
                    {
                        if(trk.rapidity()<rapMax_ && trk.rapidity()>rapMin_ && trk.pt()<=ptcut_ks_[i+1] && trk.pt()>=ptcut_ks_[i]){
                            pVect_trg_ks[i]->push_back(pvector);
                            pVect_dau_ks[i]->push_back(pvector_dau1);
                            pVect_dau_ks[i]->push_back(pvector_dau2);
                            hRap_ks[i]->Fill(trk.rapidity());
                            hPt_ks[i]->Fill(pt);
                            double KET = sqrt(mass*mass + pt*pt) - mass;
                            hKET_ks[i]->Fill(KET);
                        }
                    }
                    else
                    {
                        if(trk.eta()<=etaMax_trg_ && trk.eta()>=etaMin_trg_ && trk.pt()<=ptcut_ks_[i+1] && trk.pt()>=ptcut_ks_[i]){
                            pVect_trg_ks[i]->push_back(pvector);
                            pVect_dau_ks[i]->push_back(pvector_dau1);
                            pVect_dau_ks[i]->push_back(pvector_dau2);
                            hRap_ks[i]->Fill(trk.rapidity());
                            hPt_ks[i]->Fill(pt);
                            double KET = sqrt(mass*mass + pt*pt) - mass;
                            hKET_ks[i]->Fill(KET);
                        }
                    }
                }
            }

            if(fabs(id)==3122){

                int mid = 0;

                if(trk.numberOfMothers()==1){
                    const reco::Candidate * mom = trk.mother();
                    mid = mom->pdgId();
                    if(mom->numberOfMothers()==1){
                        const reco::Candidate * mom1 = mom->mother();
                        mid = mom1->pdgId();
                    }
                }

                if(fabs(mid)==3322 || fabs(mid)==3312 || fabs(mid)==3324 || fabs(mid)==3314 || fabs(mid)==3334) continue;

                double eta_dau1 = 0;
                double phi_dau1 = 0;
                double pt_dau1 = 999.999;

                double eta_dau2 = 0;
                double phi_dau2 = 0;
                double pt_dau2 = 999.999;

                if(trk.numberOfDaughters()==2){
                    const reco::Candidate * dau1 = trk.daughter(0);
                    const reco::Candidate * dau2 = trk.daughter(1);

                    eta_dau1 = dau1->eta();
                    phi_dau1 = dau1->phi();
                    pt_dau1 = dau1->pt();

                    eta_dau2 = dau2->eta();
                    phi_dau2 = dau2->phi();
                    pt_dau2 = dau2->pt();
                }


                TVector3 pvector;
                pvector.SetPtEtaPhi(pt,eta,phi);

                TVector3 pvector_dau1;
                pvector_dau1.SetPtEtaPhi(pt_dau1,eta_dau1,phi_dau1);

                TVector3 pvector_dau2;
                pvector_dau2.SetPtEtaPhi(pt_dau2,eta_dau2,phi_dau2);

                for(int i=0;i<18;i++)
                {
                    if(doRap_)
                    {
                        if(trk.rapidity()<rapMax_ && trk.rapidity()>rapMin_ && trk.pt()<=ptcut_la_[i+1] && trk.pt()>=ptcut_la_[i]){
                            pVect_trg_la[i]->push_back(pvector);
                            pVect_dau_la[i]->push_back(pvector_dau1);
                            pVect_dau_la[i]->push_back(pvector_dau2);
                            hRap_la[i]->Fill(trk.rapidity());
                            hPt_la[i]->Fill(pt);
                            double KET = sqrt(mass*mass + pt*pt) - mass;
                            hKET_la[i]->Fill(KET);
                        }
                    }
                    else
                    {
                        if(trk.eta()<=etaMax_trg_ && trk.eta()>=etaMin_trg_ && trk.pt()<=ptcut_la_[i+1] && trk.pt()>=ptcut_la_[i]){
                            pVect_trg_la[i]->push_back(pvector);
                            pVect_dau_la[i]->push_back(pvector_dau1);
                            pVect_dau_la[i]->push_back(pvector_dau2);
                            hRap_la[i]->Fill(trk.rapidity());
                            hPt_la[i]->Fill(pt);
                            double KET = sqrt(mass*mass + pt*pt) - mass;
                            hKET_la[i]->Fill(KET);
                        }
                    }
                }
            }

            int midXi = 0;
            if(fabs(id) == 3312)
            {
                if(trk.numberOfMothers() == 1)
                {
                    const reco::Candidate * mom = trk.mother();
                    midXi = mom->pdgId();
                    if(mom->numberOfMothers()==1){
                        const reco::Candidate * mom1 = mom->mother();
                        midXi = mom1->pdgId();
                    }
                }

                TVector3 pvector;
                pvector.SetPtEtaPhi(pt,eta,phi);

                if(fabs(midXi) != 3334){
                    for(int i=0; i<10; i++)
                    {
                        if(doRap_)
                        {
                            if(trk.rapidity()<rapMax_ && trk.rapidity()>rapMin_ && trk.pt()<=ptcut_xi_[i+1] && trk.pt()>=ptcut_xi_[i]){
                                pVect_trg_xi[i]->push_back(pvector);
                                hRap_xi[i]->Fill(trk.rapidity());
                                hPt_xi[i]->Fill(pt);
                                double KET = sqrt(mass*mass + pt*pt) - mass;
                                hKET_xi[i]->Fill(KET);
                            }
                        }
                        else
                        {
                            if(trk.eta()<=etaMax_trg_ && trk.eta()>=etaMin_trg_ && trk.pt()<=ptcut_xi_[i+1] && trk.pt()>=ptcut_xi_[i]){
                                pVect_trg_xi[i]->push_back(pvector);
                                hRap_xi[i]->Fill(trk.rapidity());
                                hPt_xi[i]->Fill(pt);
                                double KET = sqrt(mass*mass + pt*pt) - mass;
                                hKET_xi[i]->Fill(KET);
                            }
                        }
                    }
                }
            }
        }
        //Calculating signal
        int nMult_ass = (int)pVect_ass->size();
        hMult_ass->Fill(nMult_ass);

        for(int i=0; i<18; i++)
        {
            int nMult_trg_ks = (int)pVect_trg_ks[i]->size();
            int nMult_trg_la = (int)pVect_trg_la[i]->size();
            int nMult_trg_xi = (int)pVect_trg_xi[i]->size();
            hMult_ks[i]->Fill(nMult_trg_ks);
            hMult_la[i]->Fill(nMult_trg_la);
            hMult_xi[i]->Fill(nMult_trg_xi);
            for(int ntrg=0;ntrg<nMult_trg_ks;ntrg++)
            {
                TVector3 pvector_trg = (*pVect_trg_ks[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                //double pt_trg = pvector_trg.Pt();
                /*double effks = 1.0;

                  if(pt_trg<1.2)
                  {
                  effks = effkss(pt_trg);
                  }
                  else
                  {
                  effks = effksb(pt_trg);
                  }*/

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

                    if(eta_ass==eta_trg_dau1 && phi_ass==phi_trg_dau1) continue;
                    if(eta_ass==eta_trg_dau2 && phi_ass==phi_trg_dau2) continue;

                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;

                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_ks[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_ks);
                }
            }

            for(int ntrg=0;ntrg<nMult_trg_la;ntrg++)
            {
                TVector3 pvector_trg = (*pVect_trg_la[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                //double pt_trg = pvector_trg.Pt();
                /*double effla = 1.0;

                  if(pt_trg<1.6)
                  {
                  effla = efflas(pt_trg);
                  }
                  else
                  {
                  effla = efflab(pt_trg);
                  }*/


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

                    if(eta_ass==eta_trg_dau1 && phi_ass==phi_trg_dau1) continue;
                    if(eta_ass==eta_trg_dau2 && phi_ass==phi_trg_dau2) continue;

                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;

                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_la[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_la);
                }
            }

            for(int ntrg=0;ntrg<nMult_trg_xi;ntrg++)
            {
                TVector3 pvector_trg = (*pVect_trg_xi[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                //double pt_trg = pvector_trg.Pt();
                /*double effks = 1.0;

                  if(pt_trg<1.2)
                  {
                  effks = effkss(pt_trg);
                  }
                  else
                  {
                  effks = effksb(pt_trg);
                  }*/

                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();

                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;

                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_xi[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_xi);
                }
            }
        }


        for(int i=0; i<18; i++)
        {
            pVectVect_trg_ks[i]->push_back(*pVect_trg_ks[i]);
            pVectVect_trg_la[i]->push_back(*pVect_trg_la[i]);
            pVectVect_trg_xi[i]->push_back(*pVect_trg_xi[i]);
            delete pVect_trg_ks[i];
            delete pVect_trg_la[i];
            delete pVect_trg_xi[i];
            delete pVect_dau_ks[i];
            delete pVect_dau_la[i];
        }

        pVectVect_ass->push_back(*pVect_ass);
        zvtxVect->push_back(bestvz);
        delete pVect_ass;
    }
}


// ------------ method called once each job just before starting event
//loop  ------------
void
V0CorrelationMC::beginJob()
{
    edm::Service<TFileService> fs;

    TH1D::SetDefaultSumw2();


    hMult = fs->make<TH1D>("mult",";N",300,0,300);
    hMult_ass = fs->make<TH1D>("mult_ass",";N",600,0,600);

    for(int i=0; i<18; i++)
    {
        hRap_ks[i] = fs->make<TH1D>(Form("Rapkshort_pt%d",i),";y",30,-1.5,1.5);
        hRap_la[i] = fs->make<TH1D>(Form("Raplambda_pt%d",i),";y",30,-1.5,1.5);
        hRap_xi[i] = fs->make<TH1D>(Form("Rapxi_pt%d",i),";y",30,-1.5,1.5);
        hKET_ks[i] = fs->make<TH1D>(Form("KETkshort_pt%d",i),";GeV",25000,0,12.5);
        hKET_la[i] = fs->make<TH1D>(Form("KETlambda_pt%d",i),";GeV",25000,0,12.5);
        hKET_xi[i] = fs->make<TH1D>(Form("KETxi_pt%d",i),";GeV",25000,0,12.5);
        hPt_ks[i] = fs->make<TH1D>(Form("Ptkshort_pt%d",i),";GeV",25000,0,12.5);
        hPt_la[i] = fs->make<TH1D>(Form("Ptlambda_pt%d",i),";GeV",25000,0,12.5);
        hPt_xi[i] = fs->make<TH1D>(Form("Ptxi_pt%d",i),";GeV",25000,0,12.5);
        hMult_ks[i] = fs->make<TH1D>(Form("mult_ks_pt%d",i),";N",250,0,250);
        hSignal_ks[i] = fs->make<TH2D>(Form("signalkshort_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_ks[i] = fs->make<TH2D>(Form("backgroundkshort_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hMult_la[i] = fs->make<TH1D>(Form("mult_la_pt%d",i),";N",250,0,250);
        hSignal_la[i] = fs->make<TH2D>(Form("signallambda_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_la[i] = fs->make<TH2D>(Form("backgroundlambda_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hMult_xi[i] = fs->make<TH1D>(Form("mult_xi_pt%d",i),";N",250,0,250);
        hSignal_xi[i] = fs->make<TH2D>(Form("signalxi_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_xi[i] = fs->make<TH2D>(Form("backgroundxi_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        pVectVect_trg_ks[i] = new vector< vector<TVector3> >;
        pVectVect_trg_la[i] = new vector< vector<TVector3> >;
        pVectVect_trg_xi[i] = new vector< vector<TVector3> >;
    }

    pVectVect_ass = new vector< vector<TVector3> >;
    zvtxVect = new vector<double>;

}

// ------------ method called once each job just after ending the event
//loop  ------------
void
V0CorrelationMC::endJob() {
    //Calculating background
    int nevttotal_ass = (int)pVectVect_ass->size();

    for(int i=0;i<18;i++)
    {
        int nevttotal_trg_ks = (int)pVectVect_trg_ks[i]->size();
        int nevttotal_trg_la = (int)pVectVect_trg_la[i]->size();
        int nevttotal_trg_xi = (int)pVectVect_trg_xi[i]->size();

        for(int nround=0;nround<bkgFactor_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_ks; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>3000) {nevt_trg++; ncount = 0;}
                    continue; }

                vector<TVector3> pVectTmp_trg = (*pVectVect_trg_ks[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();

                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    //double pt_trg = pvectorTmp_trg.Pt();
                    /*double effks = 1.0;

                    if(pt_trg<1.2)
                    {
                        effks = effkss(pt_trg);
                    }
                    else
                    {
                        effks = effksb(pt_trg);
                    }*/

                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;

                        if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_ks[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg);
                    }
                }
            }
        }

        for(int nround=0;nround<bkgFactor_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_la; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>3000) {nevt_trg++; ncount = 0;}
                    continue; }

                vector<TVector3> pVectTmp_trg = (*pVectVect_trg_la[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();

                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    //double pt_trg = pvectorTmp_trg.Pt();
                    /*double effla = 1.0;

                    if(pt_trg<1.6)
                    {
                        effla = efflas(pt_trg);
                    }
                    else
                    {
                        effla = efflab(pt_trg);
                    }*/

                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;

                        if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_la[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg);
                    }
                }
            }
        }

        for(int nround=0;nround<bkgFactor_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_xi; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>3000) {nevt_trg++; ncount = 0;}
                    continue; }

                vector<TVector3> pVectTmp_trg = (*pVectVect_trg_xi[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();

                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    //double pt_trg = pvectorTmp_trg.Pt();
                    /*double effks = 1.0;

                    if(pt_trg<1.2)
                    {
                        effks = effkss(pt_trg);
                    }
                    else
                    {
                        effks = effksb(pt_trg);
                    }*/

                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;

                        if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_xi[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg);
                    }
                }
            }
        }
    }
}
