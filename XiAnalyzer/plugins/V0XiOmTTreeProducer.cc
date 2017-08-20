#include "XiAnalyzer/XiAnalyzer/interface/V0XiOmTTreeProducer.h"



V0XiOmTTreeProducer::V0XiOmTTreeProducer(const edm::ParameterSet& iConfig)
{
    TH1::SetDefaultSumw2(  );
    using std::string;

    v0CollName_     = iConfig.getParameter<string>( "v0CollName" );
    v0IDName_       = iConfig.getParameter<string>( "v0IDName" );
    _vertexCollName = consumes<reco::VertexCollection>( iConfig.getParameter<edm::InputTag>( "vertexCollName"));
    _xiCollection   = consumes<reco::VertexCompositeCandidateCollection>(edm::InputTag(v0CollName_,v0IDName_,"ANASKIM"));
    _v0Collection   = consumes<reco::VertexCompositeCandidateCollection>(edm::InputTag(v0CollName_,v0IDName_,"ANASKIM"));
    _trkSrc         = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trkSrc"));
    multmin_        = iConfig.getParameter<double>("multmin");
    multmax_        = iConfig.getParameter<double>("multmax");
    doRap_          = iConfig.getParameter<bool>("doRap");
    doXi_           = iConfig.getParameter<bool>("doXi");
    doKs_           = iConfig.getParameter<bool>("doKs");
    doLa_           = iConfig.getParameter<bool>("doLa");
    doOm_           = iConfig.getParameter<bool>("doOm");
    rapMin_         = iConfig.getParameter<double>("rapMin");
    rapMax_         = iConfig.getParameter<double>("rapMax");
    //treeName_ = iConfig.getParameter<string>("treeName");
}

V0XiOmTTreeProducer::~V0XiOmTTreeProducer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
V0XiOmTTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    float piMass              = 0.13957018;
    float piMass_sigma        = piMass*1e-6;
    float piMassSquared       = piMass*piMass;
    float protonMass          = 0.938272013;
    float protonMass_sigma    = protonMass*1e-6;
    float xiMass              = 1.31486;
    float electronMass        = 0.000511;
    float electronMassSquared = electronMass*electronMass;
    float lambdaMass          = 1.115683;
    float lambdaMass_sigma    = 0.000006;
    float kshortMass          = 0.497614;

    using namespace edm;
    using namespace reco;
    using std::vector;

    // Get primary vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(_vertexCollName,vertices);
    const reco::Vertex & vtx = (*vertices)[0];
    double  bestvx=-999.9, bestvy=-999.9, bestvz=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    bestvx		 = vtx.x();
    bestvy		 = vtx.y();
    bestvz		 = vtx.z();
    bestvxError	 = vtx.xError();
    bestvyError	 = vtx.yError();
    bestvzError	 = vtx.zError();


    // Track selection

    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken( _trkSrc, tracks );
    int nTracks = 0;
    int EtaPtCutnTracks = 0;

    for( unsigned it = 0; it < tracks->size(  ); it++ )
    {
        const reco::Track & trk = ( *tracks )[it];
        math::XYZPoint bestvtx( bestvx, bestvy, bestvz );

        double dzvtx    = trk.dz( bestvtx );
        double dxyvtx   = trk.dxy( bestvtx );
        double dzerror  = sqrt( trk.dzError(  )*trk.dzError(  ) + bestvzError*bestvzError );
        double dxyerror = sqrt( trk.d0Error(  )*trk.d0Error(  ) + bestvxError*bestvyError );

        if( !trk.quality( reco::TrackBase::highPurity ) ) continue;
        if( fabs( trk.ptError(  )/trk.pt(  ) > 0.10 ) )   continue;
        if( fabs( dzvtx/dzerror ) > 3 )                   continue;
        if( fabs( dxyvtx/dxyerror ) > 3 )                 continue;

        nTracks++;
        if( fabs( trk.eta(  ) ) > 2.4 || trk.pt(  ) < 0.4 ) continue;
        EtaPtCutnTracks++;
    }

    if(doXi_)
    {
        edm::Handle< reco::VertexCompositeCandidateCollection > xicandidates;
        iEvent.getByToken(_xiCollection, xicandidates);
        if(!xicandidates.isValid()){
            return;
        }

        ESHandle<MagneticField> bFieldHandle;
        iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

        ESHandle<TransientTrackBuilder> theTTB;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTB);


        if(EtaPtCutnTracks >= multmin_ && EtaPtCutnTracks < multmax_)
        {
            //Determine number of Xi candidates with no cuts other than track selection

            int numCand_xi = 0;
            cout << "Good multiplicity" << endl;
            for( reco::VertexCompositeCandidateCollection::const_iterator xicand =
                    xicandidates->begin(); xicand != xicandidates->end();
                    xicand++)
            {
                cout << "Candidates exits" << endl;
                //double secvz = -999.9, secvx = -999.9, secvy = -999.9;

                // access daughters of Xi and Lambda

                const reco::Candidate * d1 = xicand->daughter(0); // Xi_Lambda
                const reco::Candidate * d2 = xicand->daughter(1); // Xi_Pion

                const reco::Candidate * lambda_d1 = d1->daughter(0); // Lambda_Proton
                const reco::Candidate * lambda_d2 = d1->daughter(1); // Lambda_Pion

                reco::TrackRef  dau1 = d1->get<reco::TrackRef>();
                reco::TrackRef  dau2 = d2->get<reco::TrackRef>();
                reco::TrackRef  lambda_dau1 = lambda_d1->get<reco::TrackRef>();
                reco::TrackRef  lambda_dau2 = lambda_d2->get<reco::TrackRef>();

                //pt,mass
                //double eta_xi      = xicand->eta();
                //double mass_xi     = xicand->mass();
                //double pt          = xicand->pt();
                //double rapidity_xi = xicand->rapidity();

                /*
                if(!doRap_)
                {
                    if(eta_xi > 2.4 || eta_xi < -2.4){
                        cout << "bad eta" << endl;
                        continue;
                    }
                }
                else
                {
                    if(rapidity_xi > rapMax_ || rapidity_xi < rapMin_)
                    {
                        cout << "bad rap" << endl;
                        continue;
                    }
                }
                */

                //secvz  = xicand->vz(); << endl;
                //secvx  = xicand->vx();
                //secvy  = xicand->vy();

                // Xi Impact Parameter Significance
                // Build Transient Tracks of Xi_pion, Lambda_Proton, Lambda_pion
                TransientTrack pion_lambdaTT(lambda_dau2, &(*bFieldHandle));
                TransientTrack pion_XiTT(dau2, &(*bFieldHandle));
                TransientTrack proton_lambdaTT(lambda_dau1, &(*bFieldHandle));

                // Create a kinematicParticleFactory
                KinematicParticleFactoryFromTransientTrack pFactory;

                // Initialize chi2 and ndf before kinematic fits
                float chi = 0.0;
                float ndf = 0.0;

                vector<RefCountedKinematicParticle> lamParticles;
                lamParticles.push_back(pFactory.particle(pion_lambdaTT,piMass,chi,ndf,piMass_sigma));
                lamParticles.push_back(pFactory.particle(proton_lambdaTT,protonMass,chi,ndf,protonMass_sigma));

                KinematicParticleVertexFitter fitter;
                RefCountedKinematicTree v0FitTree;
                v0FitTree = fitter.fit(lamParticles);
                if(!v0FitTree->isValid())
                {
                    cout<<"invalid v0 kinematic vertex fitting"<<endl;
                    continue;
                }
                v0FitTree->movePointerToTheTop();
                RefCountedKinematicParticle lamCand      = v0FitTree->currentParticle();
                RefCountedKinematicVertex lamDecayVertex = v0FitTree->currentDecayVertex();

                vector<RefCountedKinematicParticle> xiParticles;
                xiParticles.push_back(pFactory.particle(pion_XiTT,piMass,chi,ndf,piMass_sigma));
                xiParticles.push_back(lamCand);

                RefCountedKinematicTree xiFitTree = fitter.fit(xiParticles);
                if(!xiFitTree->isValid())
                {
                    cout<<"invalid xi kinematic vertex fitting"<<endl;
                    continue;
                }
                numCand_xi++;
            }

            //Loop over candidates ignoring track selection cuts
            for( reco::VertexCompositeCandidateCollection::const_iterator xicand =
                    xicandidates->begin(); xicand != xicandidates->end();
                    xicand++)
            {
                //double secvz = -999.9, secvx = -999.9, secvy = -999.9;

                // access daughters of Xi and Lambda

                const reco::Candidate * d1 = xicand->daughter(0); // Xi_Lambda
                const reco::Candidate * d2 = xicand->daughter(1); // Xi_Pion

                const reco::Candidate * lambda_d1 = d1->daughter(0); // Lambda_Proton
                const reco::Candidate * lambda_d2 = d1->daughter(1); // Lambda_Pion

                reco::TrackRef  dau1 = d1->get<reco::TrackRef>();
                reco::TrackRef  dau2 = d2->get<reco::TrackRef>();
                reco::TrackRef  lambda_dau1 = lambda_d1->get<reco::TrackRef>();
                reco::TrackRef  lambda_dau2 = lambda_d2->get<reco::TrackRef>();

                //pt,mass
                double eta_xi      = xicand->eta();
                double mass_xi     = xicand->mass();
                double pt_xi       = xicand->pt();
                double rapidity_xi = xicand->rapidity();

                //if(eta_xi > 2.4 || eta_xi < -2.4){
                    //cout << "bad eta" << endl;
                    //continue;
                //}

                //secvz  = xicand->vz(); << endl;
                //secvx  = xicand->vx();
                //secvy  = xicand->vy();

                // Xi Impact Parameter Significance
                // Build Transient Tracks of Xi_pion, Lambda_Proton, Lambda_pion
                TransientTrack pion_lambdaTT(lambda_dau2, &(*bFieldHandle));
                TransientTrack pion_XiTT(dau2, &(*bFieldHandle));
                TransientTrack proton_lambdaTT(lambda_dau1, &(*bFieldHandle));

                // Create a kinematicParticleFactory
                KinematicParticleFactoryFromTransientTrack pFactory;

                // Initialize chi2 and ndf before kinematic fits
                float chi = 0.0;
                float ndf = 0.0;

                vector<RefCountedKinematicParticle> lamParticles;
                lamParticles.push_back(pFactory.particle(pion_lambdaTT,piMass,chi,ndf,piMass_sigma));
                lamParticles.push_back(pFactory.particle(proton_lambdaTT,protonMass,chi,ndf,protonMass_sigma));

                KinematicParticleVertexFitter fitter;
                RefCountedKinematicTree v0FitTree;
                v0FitTree = fitter.fit(lamParticles);
                if(!v0FitTree->isValid())
                {
                    cout<<"invalid v0 kinematic vertex fitting"<<endl;
                    continue;
                }
                v0FitTree->movePointerToTheTop();
                RefCountedKinematicParticle lamCand      = v0FitTree->currentParticle();
                RefCountedKinematicVertex lamDecayVertex = v0FitTree->currentDecayVertex();

                vector<RefCountedKinematicParticle> xiParticles;
                xiParticles.push_back(pFactory.particle(pion_XiTT,piMass,chi,ndf,piMass_sigma));
                xiParticles.push_back(lamCand);

                RefCountedKinematicTree xiFitTree = fitter.fit(xiParticles);
                if(!xiFitTree->isValid())
                {
                    cout<<"invalid xi kinematic vertex fitting"<<endl;
                    continue;
                }
                xiFitTree->movePointerToTheTop();
                RefCountedKinematicParticle xiCand      = xiFitTree->currentParticle();
                RefCountedKinematicVertex xiDecayVertex = xiFitTree->currentDecayVertex();

                if(!xiCand->currentState().isValid())
                {
                    cout<<"invalid state from xi cand"<<endl;
                }

                KinematicState theCurrentXiCandKinematicState = xiCand->currentState();
                FreeTrajectoryState theXiFTS = theCurrentXiCandKinematicState.freeTrajectoryState();
                TransientTrack xiTT = (*theTTB).build(theXiFTS);


                //3D impact parameter of Xi wrt primary vertex
                float xi3DIpSigValue = -1000;
                if(xiTT.isValid())
                {
                    pair<bool,Measurement1D> xi3DIpPair = IPTools::absoluteImpactParameter3D(xiTT,vtx);
                    if(xi3DIpPair.first)
                    {
                        xi3DIpSigValue = xi3DIpPair.second.significance();
                    }
                }
                else cout << "bad xiTT" << endl;

                //3D impact parameter wrt to primary vertex for Lambda_pion,
                //Lambda_proton, Xi_Pion
                float VTrkP3DIpSigValue = -1000;
                pair<bool,Measurement1D> proton3DIpPair = IPTools::absoluteImpactParameter3D(proton_lambdaTT,vtx);
                if(proton3DIpPair.first)
                {
                    VTrkP3DIpSigValue = proton3DIpPair.second.significance();
                }
                else cout << "bad proton3dippair" << endl;

                float VTrkPi3DIpSigValue = -1000;
                pair<bool,Measurement1D> pion3DIpPair = IPTools::absoluteImpactParameter3D(pion_lambdaTT,vtx);
                if(pion3DIpPair.first)
                {
                    VTrkPi3DIpSigValue = pion3DIpPair.second.significance();
                }
                else cout << "bad pionlambda" << endl;

                float xiPi3DIpSigValue = -1000;
                pair<bool,Measurement1D> pionXi3DIpPair = IPTools::absoluteImpactParameter3D(pion_XiTT,vtx);
                if(pionXi3DIpPair.first)
                {
                    xiPi3DIpSigValue = pionXi3DIpPair.second.significance();
                }
                else cout << "bad pionxi" << endl;

                // Decay length

                typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
                typedef ROOT::Math::SVector<double, 3> SVector3;

                // Getting CovMatrix
                std::vector<double> vVtxEVec;
                vVtxEVec.push_back(lamDecayVertex->error().cxx());
                vVtxEVec.push_back(lamDecayVertex->error().cyx());
                vVtxEVec.push_back(lamDecayVertex->error().cyy());
                vVtxEVec.push_back(lamDecayVertex->error().czx());
                vVtxEVec.push_back(lamDecayVertex->error().czy());
                vVtxEVec.push_back(lamDecayVertex->error().czz());
                SMatrixSym3D vVtxCov(vVtxEVec.begin(),vVtxEVec.end() );

                std::vector<double> xiVtxEVec;
                xiVtxEVec.push_back(xiDecayVertex->error().cxx());
                xiVtxEVec.push_back(xiDecayVertex->error().cyx());
                xiVtxEVec.push_back(xiDecayVertex->error().cyy());
                xiVtxEVec.push_back(xiDecayVertex->error().czx());
                xiVtxEVec.push_back(xiDecayVertex->error().czy());
                xiVtxEVec.push_back(xiDecayVertex->error().czz());
                SMatrixSym3D xiVtxCov(xiVtxEVec.begin(),xiVtxEVec.end() );

                // Decay Lengths
                SMatrixSym3D totalCov = vVtxCov  + vtx.covariance();
                SMatrixSym3D xiCov    = xiVtxCov + vtx.covariance();

                //Xi dlos to primary vertex
                SVector3 xiFlightVector(
                        xiDecayVertex->position().x() - vtx.x(),
                        xiDecayVertex->position().y() - vtx.y(),
                        xiDecayVertex->position().z() - vtx.z()
                        );
                double xiFlightMag      = ROOT::Math::Mag(xiFlightVector);
                double xiFlightSigma    = sqrt(ROOT::Math::Similarity(xiCov, xiFlightVector))/xiFlightMag;
                double xiFlightSigValue = xiFlightMag/xiFlightSigma;

                //Lambda dlos between lam and primary vertex
                SVector3 distanceVector(
                        lamDecayVertex->position().x() - vtx.x(),
                        lamDecayVertex->position().y() - vtx.y(),
                        lamDecayVertex->position().z() - vtx.z()
                        );
                double distanceMag      = ROOT::Math::Mag(distanceVector);
                double distanceSigma    = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/distanceMag;
                double distanceSigValue = distanceMag/distanceSigma;

                xi.xi3DIpSigValue_     = xi3DIpSigValue;
                xi.xiPi3DIpSigValue_   = xiPi3DIpSigValue;
                xi.VTrkPi3DIpSigValue_ = VTrkPi3DIpSigValue;
                xi.VTrkP3DIpSigValue_  = VTrkP3DIpSigValue;
                xi.xiFlightSigValue_   = xiFlightSigValue;
                xi.distanceSigValue_   = distanceSigValue;
                xi.mass_               = mass_xi;
                xi.pt_                 = pt_xi;
                xi.n_                  = numCand_xi;
                xi.rapidity_           = rapidity_xi;
                xi.eta_                = eta_xi;

                /*
                //Skim cuts
                math::XYZPoint bestvtxpt(bestvx,bestvy,bestvz);
                //Impact Sig
                double lambda_dau2_dxy      = lambda_dau2->dxy(bestvtxpt);
                double lambda_dau2_dz       = lambda_dau2->dz(bestvtxpt);
                double lambda_dau2_dxyerror = sqrt(lambda_dau2->d0Error()*lambda_dau2->d0Error() + bestvxError*bestvyError);
                double lambda_dau2_dzerror  = sqrt(lambda_dau2->dzError()*lambda_dau2->dzError() + bestvzError*bestvzError);

                double xi_dau2_dxy      = dau2->dxy(bestvtxpt);
                double xi_dau2_dz       = dau2->dz(bestvtxpt);
                double xi_dau2_dxyerror = sqrt(dau2->d0Error()*dau2->d0Error() + bestvxError*bestvyError);
                double xi_dau2_dzerror  = sqrt(dau2->dzError()*dau2->dzError() + bestvzError*bestvzError);

                //dls

                double secvx_xi = -999, secvy_xi = -999, secvz_xi = -999;

                secvx_xi = xicand->vx();
                secvy_xi = xicand->vy();
                secvz_xi = xicand->vz();

                SVector3 distanceVector_xi(secvx_xi-bestvx,secvy_xi-bestvy,secvz_xi-bestvz);
                SMatrixSym3D totalCov_xi = bestvtx.covariance() + xicand->vertexCovariance();
                double dl_xi          = ROOT::Math::Mag(distanceVector_xi);
                double dlerror_xi     = sqrt(ROOT::Math::Similarity(totalCov_xi, distanceVector_xi))/dl_xi;

                double secvx_la = -999, secvy_la = -999, secvz_la = -999;

                secvx_la = dau1->vx();
                secvy_la = dau1->vy();
                secvz_la = dau1->vz();

                SVector3 distanceVector_la(secvx_la-bestvx,secvy_la-bestvy,secvz_la-bestvz);
                SMatrixSym3D totalCov_la = bestvtx.covariance() + d1->vertexCovariance();
                double dl_la          = ROOT::Math::Mag(distanceVector_la);
                double dlerror_la     = sqrt(ROOT::Math::Similarity(totalCov_la, distanceVector_la))/dl_la;

                xi.dauTransImpactSig_    = lambda_dau2_dxy/lambda_dau2_dxyerror;
                xi.dauLongImpactSig_     = lambda_dau2_dz/lambda_dau2_dzerror;
                xi.batDauTransImpactSig_ = xi_dau2_dxy/xi_dau2_dxyerror;
                xi.batDauLongImpactSig_  = xi_dau2_dz/xi_dau2_dzerror;
                xi.xiVtxSignificance3D_  = dl_xi/dlerror_xi;
                xi.vtxSignificance3D_    = dl_la/dlerror_la;
                */

                XiTree->Fill();
            }
        }
    }

    if(doKs_ || doLa_)
    {
        if(EtaPtCutnTracks >= multmin_ && EtaPtCutnTracks < multmax_)
        {
            edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
            iEvent.getByToken(_v0Collection, v0candidates);
            if(!v0candidates.isValid()) return;

            for( reco::VertexCompositeCandidateCollection::const_iterator v0cand = v0candidates->begin();
                    v0cand != v0candidates->end();
                    v0cand++)
            {
                double secvz=-999.9, secvx=-999.9, secvy=-999.9;

                const reco::Candidate * d1 = v0cand->daughter(0);
                const reco::Candidate * d2 = v0cand->daughter(1);

                reco::TrackRef dau1 = d1->get<reco::TrackRef>();
                reco::TrackRef dau2 = d2->get<reco::TrackRef>();

                //pt,mass
                double eta_v0 = v0cand->eta();
                double pt_v0 = v0cand->pt();
                double px_v0 = v0cand->px();
                double py_v0 = v0cand->py();
                double pz_v0 = v0cand->pz();
                double rapidity_v0 = v0cand->rapidity();
                double mass_v0 = v0cand->mass();

                /*
                if(dorap_ && (rapidity_v0 > rapMax_ || rapidity_v0 < rapMin_)) continue;
                else
                    if(eta_v0 > 2.4 || eta_v0 < -2.4 ) continue;
                    */

                secvz = v0cand->vz(); secvx = v0cand->vx(); secvy = v0cand->vy();

                //trkNHits
                int nhit1 = dau1->numberOfValidHits();
                int nhit2 = dau2->numberOfValidHits();

                //if(nhit1 <= nHitCut1_ || nhit2 <= nHitCut2_) continue;

                double pt1 = dau1->pt();
                double pt2 = dau2->pt();

                //if(pt1 <= ptCut1_ || pt2 <= ptCut2_) continue;

                //algo
                //       double algo1 = dau1->algo();
                //       double algo2 = dau2->algo();

                //dau eta
                //       double eta2 = dau2->eta();

                //DCA
                math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

                double dzbest1 = dau1->dz(bestvtx);
                double dxybest1 = dau1->dxy(bestvtx);
                double dzerror1 = sqrt(dau1->dzError()*dau1->dzError()+bestvzError*bestvzError);
                double dxyerror1 = sqrt(dau1->d0Error()*dau1->d0Error()+bestvxError*bestvyError);
                double dzos1 = dzbest1/dzerror1;
                double dxyos1 = dxybest1/dxyerror1;
                //if(fabs(dzos1) < dzSigCut1_ || fabs(dxyos1) < dxySigCut1_) continue;

                double dzbest2 = dau2->dz(bestvtx);
                double dxybest2 = dau2->dxy(bestvtx);
                double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
                double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
                double dzos2 = dzbest2/dzerror2;
                double dxyos2 = dxybest2/dxyerror2;
                //if(fabs(dzos2) < dzSigCut2_ || fabs(dxyos2) < dxySigCut2_) continue;

                //vtxChi2
                double vtxChi2 = v0cand->vertexChi2();
                //if(vtxChi2 > vtxChi2Cut_ ) continue;

                //PAngle
                TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
                TVector3 secvec(px_v0,py_v0,pz_v0);
                double agl = cos(secvec.Angle(ptosvec));
                //if(agl < cosThetaCut_) continue;

                //Decay length
                typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
                typedef ROOT::Math::SVector<double, 3> SVector3;
                SMatrixSym3D totalCov = vtx.covariance() + v0cand->vertexCovariance();
                SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
                double dl = ROOT::Math::Mag(distanceVector);
                double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
                double dlos = dl/dlerror;
                //if(dlos < decayLSigCut_) continue;

                double pd1 = d1->p();
                //       double charged1 = dau1->charge();
                double pd2 = d2->p();
                //       double charged2 = dau2->charge();

                TVector3 dauvec1(d1->px(),d1->py(),d1->pz());
                TVector3 dauvec2(d2->px(),d2->py(),d2->pz());
                TVector3 dauvecsum(dauvec1+dauvec2);

                double energyd1e = sqrt(electronMassSquared+pd1*pd1);
                double energyd2e = sqrt(electronMassSquared+pd2*pd2);
                double invmass_ee = sqrt((energyd1e+energyd2e)*(energyd1e+energyd2e)-dauvecsum.Mag2());
                //if(invmass_ee<misIDMassCutEE_) continue;

                double misIDMass_la = -999;
                double misIDMass_ks_pip = -999;
                double misIDMass_ks_ppi = -999;
                if(v0IDName_ == "Lambda")
                {
                    double massd1=piMass;
                    double massd2=piMass;
                    double energyd1 = sqrt(massd1*massd1+pd1*pd1);
                    double energyd2 = sqrt(massd2*massd2+pd2*pd2);
                    double invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
                    misIDMass_la = invmass - kshortMass;
                    //if(fabs(invmass-kshortMass)<misIDMassCut_) continue;
                }

                if(v0IDName_ == "Kshort")
                {
                    double massd1=piMass;
                    double massd2=protonMass;
                    double energyd1 = sqrt(massd1*massd1+pd1*pd1);
                    double energyd2 = sqrt(massd2*massd2+pd2*pd2);
                    double invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
                    misIDMass_ks_pip = invmass - lambdaMass;
                    //if(fabs(invmass-lambdaMass)<misIDMassCut_) continue;

                    massd2=piMass;
                    massd1=protonMass;
                    energyd1 = sqrt(massd1*massd1+pd1*pd1);
                    energyd2 = sqrt(massd2*massd2+pd2*pd2);
                    invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
                    misIDMass_ks_ppi = invmass - lambdaMass;
                    //if(fabs(invmass-lambdaMass)<misIDMassCut_) continue;
                }

                if(doKs_)
                {
                    ks.eta_ =  eta_v0;
                    ks.mass_ = mass_v0;
                    ks.pt_ = pt_v0;
                    ks.rapidity_ = rapidity_v0;
                    ks.nhit1_ = nhit1;
                    ks.nhit2_ = nhit2;
                    ks.dzSig1_ = dzos1;
                    ks.dzSig2_ = dzos2;
                    ks.dxySig1_ = dxyos1;
                    ks.dxySig2_ = dxyos2;
                    ks.vtxChi2_ = vtxChi2;
                    ks.cosTheta_ = agl;
                    ks.decayLSig_ = dlos;
                    ks.misIDMassEE_ = invmass_ee;
                    ks.misIDMassForward_ = misIDMass_ks_pip;
                    ks.misIDMassBackward_ = misIDMass_ks_ppi;

                    KsTree->Fill();
                }

                if(doLa_)
                {
                    la.eta_ =  eta_v0;
                    la.mass_ = mass_v0;
                    la.pt_ = pt_v0;
                    la.rapidity_ = rapidity_v0;
                    la.nhit1_ = nhit1;
                    la.nhit2_ = nhit2;
                    la.dzSig1_ = dzos1;
                    la.dzSig2_ = dzos2;
                    la.dxySig1_ = dxyos1;
                    la.dxySig2_ = dxyos2;
                    la.vtxChi2_ = vtxChi2;
                    la.cosTheta_ = agl;
                    la.decayLSig_ = dlos;
                    la.misIDMassEE_ = invmass_ee;
                    la.misIDMassForward_ = misIDMass_la;

                    LaTree->Fill();
                }
            }
        }
    }
}

// ------------ method called once each job just before starting event
//loop  ------------
void
V0XiOmTTreeProducer::beginJob()
{
    edm::Service<TFileService> fs;
    if(doXi_)
    {
        XiTree = fs->make<TTree>("XiTree","CutParameters");
        XiTree->Branch("xi3dipsig",&xi.xi3DIpSigValue_,"xi3dipsig/F");
        XiTree->Branch("xipi3dipsig",&xi.xiPi3DIpSigValue_,"xipi3dipsig/F");
        XiTree->Branch("vtrkpi3dipsig",&xi.VTrkPi3DIpSigValue_,"vtrkpi3dipsig/F");
        XiTree->Branch("vtrkp3dipsig",&xi.VTrkP3DIpSigValue_,"vtrkp3dipsigpt/F");
        XiTree->Branch("xiflightsig",&xi.xiFlightSigValue_,"xiflightsig/F");
        XiTree->Branch("distancesig",&xi.distanceSigValue_,"distancesig/F");
        XiTree->Branch("mass",&xi.mass_,"mass/F");
        XiTree->Branch("rapidity",&xi.rapidity_,"rapidity/F");
        XiTree->Branch("eta",&xi.eta_,"eta/F");
        XiTree->Branch("pt",&xi.pt_,"pt/F");
        XiTree->Branch("n",&xi.n_,"n/I");
        /*
           XiTree->Branch("dautransimpactsig",&xi.dauTransImpactSig_,"dautransimpactsig/F");
           XiTree->Branch("daulongimpactsig",&xi.dauLongImpactSig_,"daulongimpactsig/F");
           XiTree->Branch("batdautransimpactsig",&xi.batDauTransImpactSig_,"batdautransimpactsig/F");
           XiTree->Branch("batdaulongimpactsig",&xi.batDauTransImpactSig_,"batdaulongimpactsig/F");
           XiTree->Branch("xivtxsignificance3d",&xi.xiVtxSignificance3D_,"xivtxsignificance3d/F");
           XiTree->Branch("vtxsignificance3d",&xi.vtxSignificance3D_,"vtxsignificance3d/F");
           */
    }

    if(doKs_)
    {
        KsTree = fs->make<TTree>("KsTree","CutParameters");
        KsTree->Branch("eta",&ks.eta_,"eta/F");
        KsTree->Branch("mass",&ks.mass_,"mass/F");
        KsTree->Branch("pt",&ks.pt_,"pt/F");
        KsTree->Branch("rapidity",&ks.rapidity_,"rapidity/F");
        KsTree->Branch("dzSig1",&ks.dzSig1_,"dzSig1/F");
        KsTree->Branch("dzSig2",&ks.dzSig2_,"dzSig2/F");
        KsTree->Branch("dxySig1",&ks.dxySig1_,"dxySig1/F");
        KsTree->Branch("dxySig2",&ks.dxySig2_,"dxySig2/F");
        KsTree->Branch("vtxChi2",&ks.vtxChi2_,"vtxChi2/F");
        KsTree->Branch("cosTheta",&ks.cosTheta_,"cosTheta/F");
        KsTree->Branch("decayLSig",&ks.decayLSig_,"decayLSig/F");
        KsTree->Branch("misIDMassEE",&ks.misIDMassEE_,"misIDMassEE/F");
        KsTree->Branch("misIDMasspip",&ks.misIDMassForward_,"misIDMassForward/F");
        KsTree->Branch("misIDMassppi",&ks.misIDMassBackward_,"misIDMassBackward/F");
        KsTree->Branch("nhit1",&ks.nhit1_,"nhit1/I");
        KsTree->Branch("nhit2",&ks.nhit2_,"nhit2/I");
    }

    if(doLa_)
    {
        LaTree = fs->make<TTree>("LaTree","CutParameters");
        LaTree->Branch("eta",&la.eta_,"eta/F");
        LaTree->Branch("mass",&la.mass_,"mass/F");
        LaTree->Branch("pt",&la.pt_,"pt/F");
        LaTree->Branch("rapidity",&la.rapidity_,"rapidity/F");
        LaTree->Branch("dzSig1",&la.dzSig1_,"dzSig1/F");
        LaTree->Branch("dzSig2",&la.dzSig2_,"dzSig2/F");
        LaTree->Branch("dxySig1",&la.dxySig1_,"dxySig1/F");
        LaTree->Branch("dxySig2",&la.dxySig2_,"dxySig2/F");
        LaTree->Branch("vtxChi2",&la.vtxChi2_,"vtxChi2/F");
        LaTree->Branch("cosTheta",&la.cosTheta_,"cosTheta/F");
        LaTree->Branch("decayLSig",&la.decayLSig_,"decayLSig/F");
        LaTree->Branch("misIDMassEE",&la.misIDMassEE_,"misIDMassEE/F");
        LaTree->Branch("misIDMass",&la.misIDMassForward_,"misIDMassForward/F");
        LaTree->Branch("nhit1",&la.nhit1_,"nhit1/I");
        LaTree->Branch("nhit2",&la.nhit2_,"nhit2/I");
    }
}

// ------------ method called once each job just after ending the event
//loop  ------------
void
V0XiOmTTreeProducer::endJob(){}
