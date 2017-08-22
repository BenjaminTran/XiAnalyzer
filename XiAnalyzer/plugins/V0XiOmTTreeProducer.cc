#include "XiAnalyzer/XiAnalyzer/interface/V0XiOmTTreeProducer.h"



V0XiOmTTreeProducer::V0XiOmTTreeProducer(const edm::ParameterSet& iConfig)
{
    TH1::SetDefaultSumw2(  );
    using std::string;

    v0CollName_     = iConfig.getParameter<string>( "v0CollName" );
    v0IDName_       = iConfig.getParameter<string>( "v0IDName" );
    _vertexCollName = consumes<reco::VertexCollection>( iConfig.getParameter<edm::InputTag>( "vertexCollName"));
    _xiCollection   = consumes<reco::VertexCompositeCandidateCollection>(edm::InputTag(v0CollName_,v0IDName_,"ANASKIM"));
    _omCollection   = consumes<reco::VertexCompositeCandidateCollection>(edm::InputTag(v0CollName_,v0IDName_,"ANASKIM"));
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
    float kaonMass            = 0.493677;
    float kaonMass_sigma      = kaonMass*1.e-6;

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

                XI.xi3DIpSigValue_     = xi3DIpSigValue;
                XI.xiPi3DIpSigValue_   = xiPi3DIpSigValue;
                XI.VTrkPi3DIpSigValue_ = VTrkPi3DIpSigValue;
                XI.VTrkP3DIpSigValue_  = VTrkP3DIpSigValue;
                XI.xiFlightSigValue_   = xiFlightSigValue;
                XI.distanceSigValue_   = distanceSigValue;
                XI.mass_               = mass_xi;
                XI.pt_                 = pt_xi;
                XI.n_                  = numCand_xi;
                XI.rapidity_           = rapidity_xi;
                XI.eta_                = eta_xi;

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

                XI.dauTransImpactSig_    = lambda_dau2_dxy/lambda_dau2_dxyerror;
                XI.dauLongImpactSig_     = lambda_dau2_dz/lambda_dau2_dzerror;
                XI.batDauTransImpactSig_ = xi_dau2_dxy/xi_dau2_dxyerror;
                XI.batDauLongImpactSig_  = xi_dau2_dz/xi_dau2_dzerror;
                XI.xiVtxSignificance3D_  = dl_xi/dlerror_xi;
                XI.vtxSignificance3D_    = dl_la/dlerror_la;
                */

                XiTree->Fill();
            }
        }
    }

    if(doOm_)
    {
        if(EtaPtCutnTracks >= multmin_ && EtaPtCutnTracks < multmax_)
        {
            ESHandle<MagneticField> bFieldHandle;
            iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

            ESHandle<TransientTrackBuilder> theTTB;
            iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTB);

            // get vertex
            edm::Handle<reco::VertexCollection> vertices;
            iEvent.getByToken(_vertexCollName,vertices);
            //double bestvz = -999.9, bestvx = -999.9, bestvy = -999.9;
            //double bestvzError = -999.9, bestvxError = -999.9, bestvyError = -999.9;
            const reco::Vertex & bestvtx = (*vertices)[0];
            //bestvz = bestvtx.z(); bestvx = bestvtx.x(); bestvy = bestvtx.y();
            //bestvzError = bestvtx.zError(); bestvxError = bestvtx.xError(); bestvyError = bestvtx.yError();

            // Z Vertex cut
            //if(bestvz > zVertexHigh_ || bestvz < zVertexLow_ ) return;

            edm::Handle< reco::VertexCompositeCandidateCollection > v0candidates;
            iEvent.getByToken(_OmCollection, v0candidates);
            if(!v0candidates.isValid()) return;

            // Create auto_ptr for each collection to be stored in the Event
            std::auto_ptr< reco::VertexCompositeCandidateCollection >
                theNewOmCands( new reco::VertexCompositeCandidateCollection() );


            for( reco::VertexCompositeCandidateCollection::const_iterator v0cand =
                    v0candidates->begin(); v0cand != v0candidates->end();
                    v0cand++)
            {

                //double secvz = -999.9, secvx = -999.9, secvy = -999.9;

                // access daughters of Xi and Lambda

                const reco::Candidate * d1 = v0cand->daughter(0); // Om_Lambda
                const reco::Candidate * d2 = v0cand->daughter(1); // Om_Kaon

                const reco::Candidate * lambda_d1 = d1->daughter(0); // Lambda_Proton
                const reco::Candidate * lambda_d2 = d1->daughter(1); // Lambda_Pion

                reco::TrackRef  dau1 = d1->get<reco::TrackRef>();
                reco::TrackRef  dau2 = d2->get<reco::TrackRef>();
                reco::TrackRef  lambda_dau1 = lambda_d1->get<reco::TrackRef>();
                reco::TrackRef  lambda_dau2 = lambda_d2->get<reco::TrackRef>();



                //pt,mass
                double eta_xi = v0cand->eta();
                double rap_xi = v0cand->rapidity();
                double mass_xi = v0cand->mass();
                double pt_xi = v0cand->pt();
                //double px_xi   = v0cand->px();
                //double py_xi   = v0cand->py();
                //double pz_xi   = v0cand->pz();

                //if(dorap_ && (rap_xi > rapMax_ || rap_xi < rapMin_)) continue;
                //else
                //if(eta_xi > etaCutMax_ || eta_xi < etaCutMin_) continue;

                //secvz  = v0cand->vz();
                //secvx  = v0cand->vx();
                //secvy  = v0cand->vy();

                // Xi Impact Parameter Significance
                // Build Transient Tracks of Xi_pion, Lambda_Proton, Lambda_pion
                TransientTrack pion_lambdaTT(lambda_dau2, &(*bFieldHandle));
                TransientTrack kaon_OmTT(dau2, &(*bFieldHandle));
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
                xiParticles.push_back(pFactory.particle(kaon_OmTT,kaonMass,chi,ndf,kaonMass_sigma));
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
                    pair<bool,Measurement1D> xi3DIpPair = IPTools::absoluteImpactParameter3D(xiTT,bestvtx);
                    if(xi3DIpPair.first)
                    {
                        xi3DIpSigValue = xi3DIpPair.second.significance();
                    }
                }

                //3D impact parameter wrt to primary vertex for Lambda_pion,
                //Lambda_proton, Xi_Pion
                float VTrkP3DIpSigValue = -1000;
                pair<bool,Measurement1D> proton3DIpPair = IPTools::absoluteImpactParameter3D(proton_lambdaTT,bestvtx);
                if(proton3DIpPair.first)
                {
                    VTrkP3DIpSigValue = proton3DIpPair.second.significance();
                }

                float VTrkPi3DIpSigValue = -1000;
                pair<bool,Measurement1D> pion3DIpPair = IPTools::absoluteImpactParameter3D(pion_lambdaTT,bestvtx);
                if(pion3DIpPair.first)
                {
                    VTrkPi3DIpSigValue = pion3DIpPair.second.significance();
                }

                float xiPi3DIpSigValue = -1000;
                pair<bool,Measurement1D> kaonOm3DIpPair = IPTools::absoluteImpactParameter3D(kaon_OmTT,bestvtx);
                if(kaonOm3DIpPair.first)
                {
                    xiPi3DIpSigValue = kaonOm3DIpPair.second.significance();
                }

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
                SMatrixSym3D totalCov = vVtxCov  + bestvtx.covariance();
                SMatrixSym3D xiCov    = xiVtxCov + bestvtx.covariance();

                //Xi dlos to primary vertex
                SVector3 xiFlightVector(
                        xiDecayVertex->position().x() - bestvtx.x(),
                        xiDecayVertex->position().y() - bestvtx.y(),
                        xiDecayVertex->position().z() - bestvtx.z()
                        );
                double xiFlightMag      = ROOT::Math::Mag(xiFlightVector);
                double xiFlightSigma    = sqrt(ROOT::Math::Similarity(xiCov, xiFlightVector))/xiFlightMag;
                double xiFlightSigValue = xiFlightMag/xiFlightSigma;

                //Lambda dlos between lam and primary vertex
                SVector3 distanceVector(
                        lamDecayVertex->position().x() - bestvtx.x(),
                        lamDecayVertex->position().y() - bestvtx.y(),
                        lamDecayVertex->position().z() - bestvtx.z()
                        );
                double distanceMag      = ROOT::Math::Mag(distanceVector);
                double distanceSigma    = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/distanceMag;
                double distanceSigValue = distanceMag/distanceSigma;


                OM.xi3DIpSigValue_     = xi3DIpSigValue;
                OM.xiPi3DIpSigValue_   = xiPi3DIpSigValue;
                OM.VTrkPi3DIpSigValue_ = VTrkPi3DIpSigValue;
                OM.VTrkP3DIpSigValue_  = VTrkP3DIpSigValue;
                OM.xiFlightSigValue_   = xiFlightSigValue;
                OM.distanceSigValue_   = distanceSigValue;
                OM.mass_               = mass_xi;
                OM.pt_                 = pt_xi;
                OM.eta_                = eta_xi;
                OM.rapidity_           = rapidity_xi;
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
                        KS.eta_ =  eta_v0;
                        KS.mass_ = mass_v0;
                        KS.pt_ = pt_v0;
                        KS.rapidity_ = rapidity_v0;
                        KS.nhit1_ = nhit1;
                        KS.nhit2_ = nhit2;
                        KS.dzSig1_ = dzos1;
                        KS.dzSig2_ = dzos2;
                        KS.dxySig1_ = dxyos1;
                        KS.dxySig2_ = dxyos2;
                        KS.vtxChi2_ = vtxChi2;
                        KS.cosTheta_ = agl;
                        KS.decayLSig_ = dlos;
                        KS.misIDMassEE_ = invmass_ee;
                        KS.misIDMassForward_ = misIDMass_ks_pip;
                        KS.misIDMassBackward_ = misIDMass_ks_ppi;

                        KsTree->Fill();
                    }

                    if(doLa_)
                    {
                        LA.eta_ =  eta_v0;
                        LA.mass_ = mass_v0;
                        LA.pt_ = pt_v0;
                        LA.rapidity_ = rapidity_v0;
                        LA.nhit1_ = nhit1;
                        LA.nhit2_ = nhit2;
                        LA.dzSig1_ = dzos1;
                        LA.dzSig2_ = dzos2;
                        LA.dxySig1_ = dxyos1;
                        LA.dxySig2_ = dxyos2;
                        LA.vtxChi2_ = vtxChi2;
                        LA.cosTheta_ = agl;
                        LA.decayLSig_ = dlos;
                        LA.misIDMassEE_ = invmass_ee;
                        LA.misIDMassForward_ = misIDMass_la;

                        LaTree->Fill();
                    }
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
        XiTree = fs->make<TTree>("XiTree","XiCutParameters");
        XiTree->Branch("xi3dipsig",&XI.xi3DIpSigValue_,"xi3dipsig/F");
        XiTree->Branch("xipi3dipsig",&XI.xiPi3DIpSigValue_,"xipi3dipsig/F");
        XiTree->Branch("vtrkpi3dipsig",&XI.VTrkPi3DIpSigValue_,"vtrkpi3dipsig/F");
        XiTree->Branch("vtrkp3dipsig",&XI.VTrkP3DIpSigValue_,"vtrkp3dipsigpt/F");
        XiTree->Branch("xiflightsig",&XI.xiFlightSigValue_,"xiflightsig/F");
        XiTree->Branch("distancesig",&XI.distanceSigValue_,"distancesig/F");
        XiTree->Branch("mass",&XI.mass_,"mass/F");
        XiTree->Branch("rapidity",&XI.rapidity_,"rapidity/F");
        XiTree->Branch("eta",&XI.eta_,"eta/F");
        XiTree->Branch("pt",&XI.pt_,"pt/F");
        XiTree->Branch("n",&XI.n_,"n/I");
        /*
           XiTree->Branch("dautransimpactsig",&XI.dauTransImpactSig_,"dautransimpactsig/F");
           XiTree->Branch("daulongimpactsig",&XI.dauLongImpactSig_,"daulongimpactsig/F");
           XiTree->Branch("batdautransimpactsig",&XI.batDauTransImpactSig_,"batdautransimpactsig/F");
           XiTree->Branch("batdaulongimpactsig",&XI.batDauTransImpactSig_,"batdaulongimpactsig/F");
           XiTree->Branch("xivtxsignificance3d",&XI.xiVtxSignificance3D_,"xivtxsignificance3d/F");
           XiTree->Branch("vtxsignificance3d",&XI.vtxSignificance3D_,"vtxsignificance3d/F");
           */
    }

    if(doOm_)
    {
        OmTree = fs->make<TTree>("OmTree","OmCutParameters");
        OmTree->Branch("om3dipsig",&OM.xi3DIpSigValue_,"om3dipsig/F");
        OmTree->Branch("omKaon3dipsig",&OM.xiPi3DIpSigValue_,"ompi3dipsig/F");
        OmTree->Branch("vtrkpi3dipsig",&OM.VTrkPi3DIpSigValue_,"vtrkpi3dipsig/F");
        OmTree->Branch("vtrkp3dipsig",&OM.VTrkP3DIpSigValue_,"vtrkp3dipsigpt/F");
        OmTree->Branch("omflightsig",&OM.xiFlightSigValue_,"omflightsig/F");
        OmTree->Branch("distancesig",&OM.distanceSigValue_,"distancesig/F");
        OmTree->Branch("mass",&OM.mass_,"mass/F");
        OmTree->Branch("rapidity",&OM.rapidity_,"rapidity/F");
        OmTree->Branch("eta",&OM.eta_,"eta/F");
        OmTree->Branch("pt",&OM.pt_,"pt/F");
    }


    if(doKs_)
    {
        KsTree = fs->make<TTree>("KsTree","KsCutParameters");
        KsTree->Branch("eta",&KS.eta_,"eta/F");
        KsTree->Branch("mass",&KS.mass_,"mass/F");
        KsTree->Branch("pt",&KS.pt_,"pt/F");
        KsTree->Branch("rapidity",&KS.rapidity_,"rapidity/F");
        KsTree->Branch("dzSig1",&KS.dzSig1_,"dzSig1/F");
        KsTree->Branch("dzSig2",&KS.dzSig2_,"dzSig2/F");
        KsTree->Branch("dxySig1",&KS.dxySig1_,"dxySig1/F");
        KsTree->Branch("dxySig2",&KS.dxySig2_,"dxySig2/F");
        KsTree->Branch("vtxChi2",&KS.vtxChi2_,"vtxChi2/F");
        KsTree->Branch("cosTheta",&KS.cosTheta_,"cosTheta/F");
        KsTree->Branch("decayLSig",&KS.decayLSig_,"decayLSig/F");
        KsTree->Branch("misIDMassEE",&KS.misIDMassEE_,"misIDMassEE/F");
        KsTree->Branch("misIDMasspip",&KS.misIDMassForward_,"misIDMassForward/F");
        KsTree->Branch("misIDMassppi",&KS.misIDMassBackward_,"misIDMassBackward/F");
        KsTree->Branch("nhit1",&KS.nhit1_,"nhit1/I");
        KsTree->Branch("nhit2",&KS.nhit2_,"nhit2/I");
    }

    if(doLa_)
    {
        LaTree = fs->make<TTree>("LaTree","LaCutParameters");
        LaTree->Branch("eta",&LA.eta_,"eta/F");
        LaTree->Branch("mass",&LA.mass_,"mass/F");
        LaTree->Branch("pt",&LA.pt_,"pt/F");
        LaTree->Branch("rapidity",&LA.rapidity_,"rapidity/F");
        LaTree->Branch("dzSig1",&LA.dzSig1_,"dzSig1/F");
        LaTree->Branch("dzSig2",&LA.dzSig2_,"dzSig2/F");
        LaTree->Branch("dxySig1",&LA.dxySig1_,"dxySig1/F");
        LaTree->Branch("dxySig2",&LA.dxySig2_,"dxySig2/F");
        LaTree->Branch("vtxChi2",&LA.vtxChi2_,"vtxChi2/F");
        LaTree->Branch("cosTheta",&LA.cosTheta_,"cosTheta/F");
        LaTree->Branch("decayLSig",&LA.decayLSig_,"decayLSig/F");
        LaTree->Branch("misIDMassEE",&LA.misIDMassEE_,"misIDMassEE/F");
        LaTree->Branch("misIDMass",&LA.misIDMassForward_,"misIDMassForward/F");
        LaTree->Branch("nhit1",&LA.nhit1_,"nhit1/I");
        LaTree->Branch("nhit2",&LA.nhit2_,"nhit2/I");
    }
}

// ------------ method called once each job just after ending the event
//loop  ------------
void
V0XiOmTTreeProducer::endJob(){}
