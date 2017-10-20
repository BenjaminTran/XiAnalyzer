#include "XiAnalyzer/XiAnalyzer/interface/HadronCorrelationGen.h"

HadronCorrelationGen::HadronCorrelationGen(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    etaMin_trg_ = iConfig.getParameter<double>("etaMin_trg");
    etaMax_trg_ = iConfig.getParameter<double>("etaMax_trg");
    etaMin_ass_ = iConfig.getParameter<double>("etaMin_ass");
    etaMax_ass_ = iConfig.getParameter<double>("etaMax_ass");
    ptMin_ass_ = iConfig.getParameter<double>("ptMin_ass");
    ptMax_ass_ = iConfig.getParameter<double>("ptMax_ass");
    bkgFactor_ = iConfig.getParameter<int>("bkgFactor");
    numPtBins_ = iConfig.getParameter<int>("numPtBins");
    multMax_ = iConfig.getParameter<double>("multMax");
    multMin_ = iConfig.getParameter<double>("multMin");
    rapMax_ = iConfig.getParameter<double>("rapMax");
    rapMin_ = iConfig.getParameter<double>("rapMin");
    ptcut_ = iConfig.getParameter<std::vector<double> >("ptcut");
    doGen_ = iConfig.getParameter<bool>("doGen");
    _trkSrc         = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trkSrc"));
    _vertexCollName = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollName"));
    _gnCollection   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gnCollection"));

}


HadronCorrelationGen::~HadronCorrelationGen()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HadronCorrelationGen::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    pepVect_trkass = new vector<TVector3>;

    // select on requirement of valid vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(_vertexCollName,vertices);
    double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    
    if(bestvz < -15 || bestvz>15) return;

    for(int i=0;i<numPtBins_;i++)
    {
        pVect_trg[i] = new vector<TVector3>;
    }
    
    pVect_ass = new vector<TVector3>;
    //----- loop over tracks -----
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(_trkSrc,tracks);
    
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
        double phi = trk.phi();

        nMult_ass_good++;

        TVector3 primPEPvector;
        primPEPvector.SetPtEtaPhi(pt,eta,phi);
        if(eta <= etaMax_ass_  &&
                eta >= etaMin_ass_ &&
                pt <= ptMax_ass_ && //pT associated window
                pt >= ptMin_ass_ //pT associated window
          )
            pepVect_trkass->push_back(primPEPvector);
    }
    hMult_selected->Fill(nMult_ass_good);
    
    edm::Handle<reco::GenParticleCollection> genpars;
    if(doGen_){
        iEvent.getByToken(_gnCollection,genpars);
    }
   
    if(nMult_ass_good<multMax_ && nMult_ass_good>=multMin_){
        if(doGen_){
            for(unsigned it=0; it<genpars->size(); ++it){

                const reco::GenParticle & trk = (*genpars)[it];

                double eta = trk.eta();
                double phi = trk.phi();
                double pt  = trk.pt();
                double rap = trk.rapidity();
                int st = trk.status();

                TVector3 pvector;
                pvector.SetPtEtaPhi(pt,eta,phi);

                double effweight = 1;

                for(int i=0;i<numPtBins_;i++)
                {
                    if(rap<=rapMax_ && rap>=rapMin_ && pt<=ptcut_[i+1] && pt>=ptcut_[i] && fabs(trk.charge())==1 && st==1){
                        hPt[i]->Fill(pt,1.0/effweight);
                        hEta[i]->Fill(eta);
                        hPhi[i]->Fill(phi);
                        hRap[i]->Fill(rap);
                        pVect_trg[i]->push_back(pvector);
                    }
                }

                if(eta<=etaMax_ass_ && eta>=etaMin_ass_ && pt<=ptMax_ass_ && pt>=ptMin_ass_ && fabs(trk.charge())==1 && st==1) pVect_ass->push_back(pvector);
            }
        }
        else{
            for(reco::TrackCollection::const_iterator it=tracks->begin(); it!=tracks->end(); ++it){
                double eta = it->eta();
                double phi = it->phi();
                double pt = it->pt();
                //double mass = it->mass();
                //double mom = ROOT::Math::Mag(it->momentum());
                //double enery = sqrt(mass*mass + mom*mom);
                //double rap = it->rapidity();

                TVector3 pvector;
                pvector.SetPtEtaPhi(pt,eta,phi);

                for(int i=0; i<numPtBins_; i++){
                    if(rap<=rapMax_ && rap>=rapMin_ && pt<=ptcut_[i+1] && pt>=ptcut_[i]){
                        hPt[i]->Fill(pt);
                        hEta[i]->Fill(eta);
                        hPhi[i]->Fill(phi);
                        hRap[i]->Fill(rap);
                        pVect_trg[i]->push_back(pvector);
                    }
                }

                if(eta<=etaMax_ass_ && eta>=etaMin_ass_ && pt<=ptMax_ass_ && pt>=ptMin_ass_) pVect_ass->push_back(pvector);
            }
        }

        int nMult_ass = (int)pVect_ass->size();
        // Calculating signal
        for(int i=0;i<numPtBins_;i++)
        {
            int nMult_trg = (int)pVect_trg[i]->size();
            hMult[i]->Fill(nMult_trg);

            int nMult_corr = 0;
            if(nMult_trg>0 && nMult_ass>0) nMult_corr = 1;

            double nMult_trg_eff=0;

            for(int ntrg=0;ntrg<nMult_trg;ntrg++)
            {
                TVector3 pvector_trg = (*pVect_trg[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                //double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();

                double effweight_trg = 1;

                nMult_trg_eff = nMult_trg_eff + 1.0/effweight_trg;
            }

            for(int ntrg=0;ntrg<nMult_trg;ntrg++)
            {
                int nMult_ass_pair = (int)pVect_ass->size();
                TVector3 pvector_trg = (*pVect_trg[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();

                double effweight_trg = 1;

                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    double pt_ass = pvector_ass.Pt();

                    double effweight_ass = 1;

                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;

                    if(fabs(deltaEta)==0 && fabs(deltaPhi)==0){
                        if(nMult_trg==1 && nMult_ass==1){nMult_corr = 0;}
                        nMult_ass_pair = nMult_ass_pair - 1;
                        continue;}
                    hSignal[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effweight_trg/effweight_ass);
                }
                if(nMult_corr > 0) hMult_assoc[i]->Fill(nMult_ass_pair);
            }
            hMult_good[i]->Fill(nMult_corr);
            pVectVect_trg[i]->push_back(*pVect_trg[i]);
        }
        /*
        //
        // Make signal histogram for pairing of two charged primary GEN tracks
        //
        int pepVect_trkhad_size = (int)pVect_ass->size();
        HadPerEvt->Fill(pepVect_trkhad_size);

        for(int trktrg1 = 0; trktrg1 < pepVect_trkhad_size; trktrg1++)
        {
            TVector3 pepVect_had1 = (*pVect_ass)[trktrg1];
            double eta_trg = pepVect_had1.Eta();
            double phi_trg = pepVect_had1.Phi();

            for(int trktrg2 = 0; trktrg2 < pepVect_trkhad_size; trktrg2++)
            {
                if(trktrg2 == trktrg1){
                    continue;
                }
                TVector3 pepVect_had2 = (*pVect_ass)[trktrg2];
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
                SignalHad->Fill(dEta, dPhi,1.0/pepVect_trkhad_size);
            }
        }
        */


        /*
        // Make signal histogram for pairing of two charged primary RECO tracks
        //
        int pepVect_trkass_size = (int)pepVect_trkass->size();

        for(int trktrg1 = 0; trktrg1 < pepVect_trkass_size; trktrg1++)
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
                SignalHadReco->Fill(dEta, dPhi);
            }
        }
        */


        pVect2_ass->push_back(*pepVect_trkass);
        pVectVect_ass->push_back(*pVect_ass);
        zvtxVect->push_back(bestvz);

        for(int i=0; i<numPtBins_; i++)
        {
            delete pVect_trg[i];
        }

        delete pepVect_trkass;
        delete pVect_ass;
    }
}

// ------------ method called once each job just before starting event loop  ------------
void
HadronCorrelationGen::beginJob(){
    
    TH1D::SetDefaultSumw2();
    
    TFileDirectory KineParam = fs->mkdir("KineParam");
    TFileDirectory Mult = fs->mkdir("Mult");
    hMult_selected = fs->make<TH1D>("mult_selected",";N",600,0,600);
    HadPerEvt       = fs->make<TH1D>("HadPerEvent", "Hadrons per Event", 1500, 0, 1500);
    //BackgroundHad   = fs->make<TH2D>("BackgroundHad", "BkgHad; #Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    //SignalHad       = fs->make<TH2D>("SignalHad", "SigHad; #Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    //BackgroundHadReco   = fs->make<TH2D>("BackgroundHadReco", "BkgHad; #Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    //SignalHadReco       = fs->make<TH2D>("SignalHadReco", "SigHad; #Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    
    for(int i=0;i<numPtBins_;i++)
    {
        hPt[i] = KineParam.make<TH1D>(Form("Pt%d",i),";GeV",25000,0,12.5);
        hRap[i] = KineParam.make<TH1D>(Form("rap_pt%d",i),";y",100,-5,5);
        hEta[i] = KineParam.make<TH1D>(Form("Eta_pt%d",i),";#eta",100,-5.0,5.0);
        hPhi[i] = KineParam.make<TH1D>(Form("Phi_pt%d",i),";#phi",132,-2.1*PI,2.1*PI);
        hMult[i] = Mult.make<TH1D>(Form("mult%d",i),";N",300,0,300);
        hMult_good[i] = Mult.make<TH1D>(Form("mult_good%d",i),";N",300,0,300);
        hMult_assoc[i] = Mult.make<TH1D>(Form("mult_assoc%d",i),";N",300,0,300);
        hSignal[i] = fs->make<TH2D>(Form("signal%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground[i] = fs->make<TH2D>(Form("background%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        pVectVect_trg[i] = new vector< vector<TVector3> >;
    }
    
    pVect2_ass = new vector< vector<TVector3> >;
    pVectVect_ass = new vector< vector<TVector3> >;
    zvtxVect = new vector<double>;
    cout << "End of beginJob" << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void
HadronCorrelationGen::endJob() {
    cout << "End event loop" << endl;
    int nevttotal_ass = (int)pVectVect_ass->size();
    int nevttotal_ass_Reco = (int)pVect2_ass->size();
    
    // Calculate background
    for(int i=0;i<numPtBins_;i++)
    {
        int nevttotal_trg = (int)pVectVect_trg[i]->size();
        
        for(int nround=0;nround<bkgFactor_;nround++)
        {
            int ncount = 0;
            for(int nevt_ass=0; nevt_ass<nevttotal_ass; nevt_ass++)
            {
                int nevt_trg = gRandom->Integer(nevttotal_trg);
                if(nevt_trg == nevt_ass) { nevt_ass--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_ass--;
                    ncount++;
                    if(ncount>5000) {nevt_ass++; ncount = 0;}
                    continue; }
                
                vector<TVector3> pVectTmp_trg = (*pVectVect_trg[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();
                
                double nMult_trg_eff=0;
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    //double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    
                    double effweight_trg = 1;
                    
                    nMult_trg_eff = nMult_trg_eff + 1.0/effweight_trg;
                }
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    
                    double effweight_trg = 1;
                    
                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        double pt_ass = pvectorTmp_ass.Pt();
                        
                        double effweight_ass = 1;

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;
                        
                        //if(fabs(deltaEta)<0.028 && fabs(deltaPhi)<0.02) continue;
                        
                        hBackground[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effweight_trg/effweight_ass);
                    }
                }
            }
        }
    }

    /*
    //
    // hadron paired with hadron GEN
    //

    for(int bkgnum = 0; bkgnum<bkgFactor_; bkgnum++)
    {
        int ncount = 0;
        for(int nevt_ass=0; nevt_ass<nevttotal_ass; nevt_ass++)
        {
            int nevt_trg = gRandom->Integer(nevttotal_ass);
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

            vector<TVector3> pepVectTmp_trg = (*pVectVect_ass)[nevt_trg];
            vector<TVector3> pepVectTmp_ass = (*pVectVect_ass)[nevt_ass];
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

    //
    // hadron paired with hadron RECO
    //

    for(int bkgnum = 0; bkgnum<bkgFactor_; bkgnum++)
    {
        int ncount = 0;
        for(int nevt_ass=0; nevt_ass<nevttotal_ass_Reco; nevt_ass++)
        {
            int nevt_trg = gRandom->Integer(nevttotal_ass_Reco);
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

            vector<TVector3> pepVectTmp_trg = (*pVect2_ass)[nevt_trg];
            vector<TVector3> pepVectTmp_ass = (*pVect2_ass)[nevt_ass];
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

                    BackgroundHadReco->Fill(dEta, dPhi);
                }
            }
        }
    }
*/
}
