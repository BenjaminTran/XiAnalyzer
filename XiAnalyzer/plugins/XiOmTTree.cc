#include "XiAnalyzer/XiAnalyzer/interface/XiOmTTree.h"



XiOmTTree::XiOmTTree(const edm::ParameterSet& iConfig)
{
    TH1::SetDefaultSumw2(  );
    using std::string;

    v0CollName_ = iConfig.getParameter<string>( "v0CollName" );
    v0IDName_ = iConfig.getParameter<string>( "v0IDName" );
    _vertexCollName = consumes<reco::VertexCollection>( iConfig.getParameter<edm::InputTag>( "vertexCollName"));
    _xiCollection   = consumes<reco::VertexCompositeCandidateCollection>(edm::InputTag(v0CollName_,v0IDName_,"ANASKIM"));
    _trkSrc = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trkSrc"));
}

XiOmTTree::~XiOmTTree()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
XiOmTTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    using namespace edm;
    using namespace reco;
    using std::vector;

    ESHandle<MagneticField> bFieldHandle;
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

    ESHandle<TransientTrackBuilder> theTTB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTB);

	// Get primary vertex
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(_vertexCollName,vertices);
	const reco::Vertex & bestvtx = (*vertices)[0];
	double  bestvx=-999.9, bestvy=-999.9, bestvz=-999.9;
	double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
	bestvx		 = bestvtx.x();
	bestvy		 = bestvtx.y();
	bestvz		 = bestvtx.z();
	bestvxError	 = bestvtx.xError();
	bestvyError	 = bestvtx.yError();
	bestvzError	 = bestvtx.zError();

    edm::Handle< reco::VertexCompositeCandidateCollection > v0candidates;
    iEvent.getByToken(_xiCollection, v0candidates);
    if(!v0candidates.isValid()){
        return;
    }

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

    if(EtaPtCutnTracks >= 185 && EtaPtCutnTracks < 220)
    {
        //Determine number of candidates

        int numv0cand = 0;
        cout << "Good multiplicity" << endl;
        for( reco::VertexCompositeCandidateCollection::const_iterator v0cand =
                v0candidates->begin(); v0cand != v0candidates->end();
                v0cand++)
        {
        cout << "Candidates exits" << endl;
            //double secvz = -999.9, secvx = -999.9, secvy = -999.9;

            // access daughters of Xi and Lambda

            const reco::Candidate * d1 = v0cand->daughter(0); // Xi_Lambda
            const reco::Candidate * d2 = v0cand->daughter(1); // Xi_Pion

            const reco::Candidate * lambda_d1 = d1->daughter(0); // Lambda_Proton
            const reco::Candidate * lambda_d2 = d1->daughter(1); // Lambda_Pion

            reco::TrackRef  dau1 = d1->get<reco::TrackRef>();
            reco::TrackRef  dau2 = d2->get<reco::TrackRef>();
            reco::TrackRef  lambda_dau1 = lambda_d1->get<reco::TrackRef>();
            reco::TrackRef  lambda_dau2 = lambda_d2->get<reco::TrackRef>();

            //pt,mass
            double eta_xi  = v0cand->eta();
            double mass = v0cand->mass();
            double pt = v0cand->pt();

            if(eta_xi > 2.4 || eta_xi < -2.4){
                cout << "bad eta" << endl;
                continue;
            }

            //secvz  = v0cand->vz(); << endl;
            //secvx  = v0cand->vx();
            //secvy  = v0cand->vy();

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
            numv0cand++;
        }

        //Loop over candidates
        for( reco::VertexCompositeCandidateCollection::const_iterator v0cand =
                v0candidates->begin(); v0cand != v0candidates->end();
                v0cand++)
        {
            //double secvz = -999.9, secvx = -999.9, secvy = -999.9;

            // access daughters of Xi and Lambda

            const reco::Candidate * d1 = v0cand->daughter(0); // Xi_Lambda
            const reco::Candidate * d2 = v0cand->daughter(1); // Xi_Pion

            const reco::Candidate * lambda_d1 = d1->daughter(0); // Lambda_Proton
            const reco::Candidate * lambda_d2 = d1->daughter(1); // Lambda_Pion

            reco::TrackRef  dau1 = d1->get<reco::TrackRef>();
            reco::TrackRef  dau2 = d2->get<reco::TrackRef>();
            reco::TrackRef  lambda_dau1 = lambda_d1->get<reco::TrackRef>();
            reco::TrackRef  lambda_dau2 = lambda_d2->get<reco::TrackRef>();

            //pt,mass
            double eta_xi  = v0cand->eta();
            double mass = v0cand->mass();
            double pt = v0cand->pt();

            if(eta_xi > 2.4 || eta_xi < -2.4){
                cout << "bad eta" << endl;
                continue;
            }

            //secvz  = v0cand->vz(); << endl;
            //secvx  = v0cand->vx();
            //secvy  = v0cand->vy();

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
                pair<bool,Measurement1D> xi3DIpPair = IPTools::absoluteImpactParameter3D(xiTT,bestvtx);
                if(xi3DIpPair.first)
                {
                    xi3DIpSigValue = xi3DIpPair.second.significance();
                }
            }
            else cout << "bad xiTT" << endl;

            //3D impact parameter wrt to primary vertex for Lambda_pion,
            //Lambda_proton, Xi_Pion
            float VTrkP3DIpSigValue = -1000;
            pair<bool,Measurement1D> proton3DIpPair = IPTools::absoluteImpactParameter3D(proton_lambdaTT,bestvtx);
            if(proton3DIpPair.first)
            {
                VTrkP3DIpSigValue = proton3DIpPair.second.significance();
            }
            else cout << "bad proton3dippair" << endl;

            float VTrkPi3DIpSigValue = -1000;
            pair<bool,Measurement1D> pion3DIpPair = IPTools::absoluteImpactParameter3D(pion_lambdaTT,bestvtx);
            if(pion3DIpPair.first)
            {
                VTrkPi3DIpSigValue = pion3DIpPair.second.significance();
            }
            else cout << "bad pionlambda" << endl;

            float xiPi3DIpSigValue = -1000;
            pair<bool,Measurement1D> pionXi3DIpPair = IPTools::absoluteImpactParameter3D(pion_XiTT,bestvtx);
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

            xistruct_.xi3DIpSigValue_     = xi3DIpSigValue;
            xistruct_.xiPi3DIpSigValue_   = xiPi3DIpSigValue;
            xistruct_.VTrkPi3DIpSigValue_ = VTrkPi3DIpSigValue;
            xistruct_.VTrkP3DIpSigValue_  = VTrkP3DIpSigValue;
            xistruct_.xiFlightSigValue_   = xiFlightSigValue;
            xistruct_.distanceSigValue_   = distanceSigValue;
            xistruct_.mass_               = mass;
            xistruct_.pt_                 = pt;
            xistruct_.n_                  = numv0cand;

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

            secvx_xi = v0cand->vx();
            secvy_xi = v0cand->vy();
            secvz_xi = v0cand->vz();

            SVector3 distanceVector_xi(secvx_xi-bestvx,secvy_xi-bestvy,secvz_xi-bestvz);
            SMatrixSym3D totalCov_xi = bestvtx.covariance() + v0cand->vertexCovariance();
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

            xistruct_.dauTransImpactSig_    = lambda_dau2_dxy/lambda_dau2_dxyerror;
            xistruct_.dauLongImpactSig_     = lambda_dau2_dz/lambda_dau2_dzerror;
            xistruct_.batDauTransImpactSig_ = xi_dau2_dxy/xi_dau2_dxyerror;
            xistruct_.batDauLongImpactSig_  = xi_dau2_dz/xi_dau2_dzerror;
            xistruct_.xiVtxSignificance3D_  = dl_xi/dlerror_xi;
            xistruct_.vtxSignificance3D_    = dl_la/dlerror_la;
             */

            XiTree->Fill();
        }
    }
}

// ------------ method called once each job just before starting event
//loop  ------------
void
XiOmTTree::beginJob()
{
    edm::Service<TFileService> fs;
    XiTree = fs->make<TTree>("XiTree","CutParameters");
    XiTree->Branch("xi3dipsig",&xistruct_.xi3DIpSigValue_,"xi3dipsig/F");
    XiTree->Branch("xipi3dipsig",&xistruct_.xiPi3DIpSigValue_,"xipi3dipsig/F");
    XiTree->Branch("vtrkpi3dipsig",&xistruct_.VTrkPi3DIpSigValue_,"vtrkpi3dipsig/F");
    XiTree->Branch("vtrkp3dipsig",&xistruct_.VTrkP3DIpSigValue_,"vtrkp3dipsigpt/F");
    XiTree->Branch("xiflightsig",&xistruct_.xiFlightSigValue_,"xiflightsig/F");
    XiTree->Branch("distancesig",&xistruct_.distanceSigValue_,"distancesig/F");
    XiTree->Branch( "mass",&xistruct_.mass_,"mass/F" );
    XiTree->Branch( "pt",&xistruct_.pt_,"pt/F" );
    XiTree->Branch( "n",&xistruct_.n_,"n/I" );
    /*
    XiTree->Branch("dautransimpactsig",&xistruct_.dauTransImpactSig_,"dautransimpactsig/F");
    XiTree->Branch("daulongimpactsig",&xistruct_.dauLongImpactSig_,"daulongimpactsig/F");
    XiTree->Branch("batdautransimpactsig",&xistruct_.batDauTransImpactSig_,"batdautransimpactsig/F");
    XiTree->Branch("batdaulongimpactsig",&xistruct_.batDauTransImpactSig_,"batdaulongimpactsig/F");
    XiTree->Branch("xivtxsignificance3d",&xistruct_.xiVtxSignificance3D_,"xivtxsignificance3d/F");
    XiTree->Branch("vtxsignificance3d",&xistruct_.vtxSignificance3D_,"vtxsignificance3d/F");
    */

	// XiOmTTree_dist = fs->make< TNtuple>("XiOmTTree","XiOmTTree","eta:phi:pt:chi2:nhits:3Ddca:pt_err_pt:qual");
}

// ------------ method called once each job just after ending the event
//loop  ------------
void
XiOmTTree::endJob(){}
