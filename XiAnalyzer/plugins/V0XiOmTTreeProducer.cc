#include "XiAnalyzer/XiAnalyzer/interface/V0XiOmTTreeProducer.h"



V0XiOmTTreeProducer::V0XiOmTTreeProducer(const edm::ParameterSet& iConfig)
{
    TH1::SetDefaultSumw2(  );
    using std::string;

    v0CollName_     = iConfig.getParameter<string>( "v0CollName" );
    v0IDName_       = iConfig.getParameter<string>( "v0IDName" );
    _vertexCollName = consumes<reco::VertexCollection>( iConfig.getParameter<edm::InputTag>( "vertexCollName"));
    _v0Collection   = consumes<reco::VertexCompositeCandidateCollection>(edm::InputTag(v0CollName_,v0IDName_,"ANASKIM"));
    _trkSrc         = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trkSrc"));
    _towerSrc       = consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("towerSrc"));
    multLow_        = iConfig.getParameter<double>("multLow");
    multHigh_        = iConfig.getParameter<double>("multHigh");
    rapMin_         = iConfig.getParameter<double>("rapMin");
    rapMax_         = iConfig.getParameter<double>("rapMax");
    ptCut1_         = iConfig.getParameter<double>("ptCut1");
    ptCut2_         = iConfig.getParameter<double>("ptCut2");
    misIDMassCut_   = iConfig.getParameter<double>("misIDMassCut");
    misIDMassCutEE_ = iConfig.getParameter<double>("misIDMassCutEE");
    cent_bin_low_   = iConfig.getParameter<int>("cent_bin_low");
    cent_bin_high_  = iConfig.getParameter<int>("cent_bin_high");
    nHitCut1_       = iConfig.getParameter<int>("nHitCut1");
    nHitCut2_       = iConfig.getParameter<int>("nHitCut2");
    doRap_          = iConfig.getParameter<bool>("doRap");
    useCentrality_  = iConfig.getParameter<bool>("useCentrality");

    if(useCentrality_){
        multLow_ = 0;
        multHigh_ = 999999;
    }
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
    nEv->Fill(1);
    float piMass              = 0.13957018;
    float piMass_sigma        = piMass*1e-6;
    float piMassSquared       = piMass*piMass;
    float protonMass          = 0.938272013;
    float protonMass_sigma    = protonMass*1e-6;
    float xiMass              = 1.32171;
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

    if(useCentrality_)
    {
        double etHFtowerSumPlus=0;
        double etHFtowerSumMinus=0;
        double etHFtowerSum=0;
        Handle<CaloTowerCollection> towers;
        iEvent.getByToken(_towerSrc,towers);
        for( size_t i = 0; i<towers->size(); ++ i){
            const CaloTower & tower = (*towers)[ i ];
            double eta = tower.eta();
            bool isHF = tower.ietaAbs() > 29;
            if(isHF && eta > 0){
                etHFtowerSumPlus += tower.pt();
            }
            if(isHF && eta < 0){
                etHFtowerSumMinus += tower.pt();
            }
        }
        etHFtowerSum=etHFtowerSumPlus + etHFtowerSumMinus;

        double binLowEdge[200]={4487.37, 4370.52, 4279.22, 4199.03, 4113.36, 4030.29, 3947.37, 3870.82, 3797.94, 3721.48, 3648.81, 3575.99, 3506.18, 3440.74, 3374.5, 3310.49, 3249.72, 3190.49, 3127.43, 3066.91, 3012.16, 2954.08, 2897.16, 2840.3, 2786.54, 2735.06, 2682.83, 2631.95, 2580.71, 2529.93, 2483.34, 2436.59, 2389.05, 2343.58, 2300.27, 2256.49, 2210.35, 2167.14, 2128.09, 2086.24, 2044.85, 2002.72, 1962.42, 1925.23, 1889.2, 1851.68, 1815.58, 1778.47, 1743.48, 1706.47, 1671.08, 1636.7, 1604.94, 1571.63, 1539.86, 1508.37, 1477.12, 1445.73, 1417.7, 1387.98, 1359.02, 1330.3, 1301.45, 1274.07, 1246.54, 1219.36, 1191.97, 1165.77, 1140.4, 1114.92, 1091.98, 1067.94, 1043.67, 1019.66, 995.39, 970.466, 947.786, 924.75, 902.723, 879.824, 859.262, 838.212, 817.18, 796.627, 776.494, 757.142, 737.504, 719.604, 701.142, 684.043, 665.89, 648.427, 630.224, 612.877, 596.435, 580.397, 565.396, 550.272, 535.204, 520.48, 505.854, 491.648, 477.531, 463.192, 449.773, 436.806, 423.944, 410.4, 397.962, 386.135, 374.47, 362.499, 351.17, 339.635, 328.402, 317.875, 307.348, 296.957, 287.002, 276.94, 267.822, 258.796, 249.366, 239.974, 231.563, 223.362, 214.902, 206.818, 199.417, 191.609, 184.184, 177.042, 169.839, 163.579, 157.186, 151.136, 145.165, 139.213, 133.218, 127.748, 122.445, 117.458, 112.715, 108.179, 103.713, 99.2518, 94.8864, 90.7892, 86.692, 82.819, 79.0331, 75.4791, 71.8774, 68.5738, 65.5363, 62.6369, 59.7441, 57.0627, 54.3838, 51.7242, 49.1577, 46.7914, 44.4615, 42.3374, 40.2863, 38.2674, 36.3979, 34.4769, 32.7274, 30.9911, 29.3998, 27.7739, 26.2442, 24.795, 23.3496, 21.8717, 20.5263, 19.2405, 18.08, 16.9542, 15.882, 14.8344, 13.8014, 12.7824, 11.8165, 10.8308, 9.94351, 9.08363, 8.20773, 7.40535, 6.57059, 5.81859, 5.0626, 4.32634, 3.57026, 2.83467, 2.09189, 1.36834, 0.673038, 0};

        int bin = -1;
        for(int i=0; i<200; i++){
            if(etHFtowerSum>=binLowEdge[i]){
                bin = i; break;
            }
        }

        if(bin<cent_bin_low_ || bin>cent_bin_high_) return;
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


    if(v0IDName_ == "Xi" || v0IDName_ == "Omega")
    {
        edm::Handle< reco::VertexCompositeCandidateCollection > CascadeCandidates;
        iEvent.getByToken(_v0Collection, CascadeCandidates);
        if(!CascadeCandidates.isValid()){
            return;
        }

        ESHandle<MagneticField> bFieldHandle;
        iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

        ESHandle<TransientTrackBuilder> theTTB;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTB);


        if(EtaPtCutnTracks >= multLow_ && EtaPtCutnTracks < multHigh_)
        {
            //Determine number of Xi candidates with no cuts other than track selection

            //Loop over candidates ignoring track selection cuts
            for( reco::VertexCompositeCandidateCollection::const_iterator CasCand =
                    CascadeCandidates->begin(); CasCand != CascadeCandidates->end();
                    CasCand++)
            {
                //double secvz = -999.9, secvx = -999.9, secvy = -999.9;

                // access daughters of Xi and Lambda

                const reco::Candidate * d1 = CasCand->daughter(0); // Xi_Lambda
                const reco::Candidate * d2 = CasCand->daughter(1); // Xi_Pion

                const reco::Candidate * lambda_d1 = d1->daughter(0); // Lambda_Proton
                const reco::Candidate * lambda_d2 = d1->daughter(1); // Lambda_Pion

                reco::TrackRef  dau1 = d1->get<reco::TrackRef>();
                reco::TrackRef  dau2 = d2->get<reco::TrackRef>();
                reco::TrackRef  lambda_dau1 = lambda_d1->get<reco::TrackRef>();
                reco::TrackRef  lambda_dau2 = lambda_d2->get<reco::TrackRef>();

                //pt,mass
                double eta  = CasCand->eta();
                double mass = CasCand->mass();
                double pt   = CasCand->pt();
                double rap  = CasCand->rapidity();

                //if(doRap_)
                //{
                    //if(rap < rapMin_ || rap > rapMax_) continue;
                //}
                //else
                //{
                    //if(eta > 2.4 || eta < -2.4) continue;
                //}

                //secvz  = CasCand->vz(); << endl;
                //secvx  = CasCand->vx();
                //secvy  = CasCand->vy();

                // Xi Impact Parameter Significance
                // Build Transient Tracks of Xi_pion, Lambda_Proton, Lambda_pion
                TransientTrack pion_lambdaTT(lambda_dau2, &(*bFieldHandle));
                TransientTrack bat_CasTT(dau2, &(*bFieldHandle));
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

                vector<RefCountedKinematicParticle> casParticles;
                if(v0IDName_ == "Xi") casParticles.push_back(pFactory.particle(bat_CasTT,piMass,chi,ndf,piMass_sigma));
                else if(v0IDName_ == "Omega") casParticles.push_back(pFactory.particle(bat_CasTT,kaonMass,chi,ndf,kaonMass_sigma));
                else cout << "unexpected particle option" << endl;
                casParticles.push_back(lamCand);

                RefCountedKinematicTree casFitTree = fitter.fit(casParticles);
                if(!casFitTree->isValid())
                {
                    cout<<"invalid cas kinematic vertex fitting"<<endl;
                    continue;
                }
                casFitTree->movePointerToTheTop();
                RefCountedKinematicParticle casCand      = casFitTree->currentParticle();
                RefCountedKinematicVertex casDecayVertex = casFitTree->currentDecayVertex();

                if(!casCand->currentState().isValid())
                {
                    cout<<"invalid state from cas cand"<<endl;
                }

                KinematicState theCurrentXiCandKinematicState = casCand->currentState();
                FreeTrajectoryState theCasFTS = theCurrentXiCandKinematicState.freeTrajectoryState();
                TransientTrack casTT = (*theTTB).build(theCasFTS);


                //3D impact parameter of Xi wrt primary vertex
                float cas3DIpSigValue = -1000;
                if(casTT.isValid())
                {
                    pair<bool,Measurement1D> cas3DIpPair = IPTools::absoluteImpactParameter3D(casTT,vtx);
                    if(cas3DIpPair.first)
                    {
                        cas3DIpSigValue = cas3DIpPair.second.significance();
                    }
                }
                else cout << "bad casTT" << endl;

                //3D impact parameter wrt to primary vertex for Lambda_pion,
                //Lambda_proton, Cas_batchelor
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

                float casBat3DIpSigValue = -1000;
                pair<bool,Measurement1D> casBat3DIpPair = IPTools::absoluteImpactParameter3D(bat_CasTT,vtx);
                if(casBat3DIpPair.first)
                {
                    casBat3DIpSigValue = casBat3DIpPair.second.significance();
                }
                else cout << "bad casBat" << endl;

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

                std::vector<double> casVtxEVec;
                casVtxEVec.push_back(casDecayVertex->error().cxx());
                casVtxEVec.push_back(casDecayVertex->error().cyx());
                casVtxEVec.push_back(casDecayVertex->error().cyy());
                casVtxEVec.push_back(casDecayVertex->error().czx());
                casVtxEVec.push_back(casDecayVertex->error().czy());
                casVtxEVec.push_back(casDecayVertex->error().czz());
                SMatrixSym3D casVtxCov(casVtxEVec.begin(),casVtxEVec.end() );

                // Decay Lengths
                SMatrixSym3D totalCov = vVtxCov  + vtx.covariance();
                SMatrixSym3D casCov    = casVtxCov + vtx.covariance();

                //Xi dlos to primary vertex
                SVector3 casFlightVector(
                        casDecayVertex->position().x() - vtx.x(),
                        casDecayVertex->position().y() - vtx.y(),
                        casDecayVertex->position().z() - vtx.z()
                        );
                double casFlightMag      = ROOT::Math::Mag(casFlightVector);
                double casFlightSigma    = sqrt(ROOT::Math::Similarity(casCov, casFlightVector))/casFlightMag;
                double casFlightSigValue = casFlightMag/casFlightSigma;

                //Lambda dlos between lam and primary vertex
                SVector3 distanceVector(
                        lamDecayVertex->position().x() - vtx.x(),
                        lamDecayVertex->position().y() - vtx.y(),
                        lamDecayVertex->position().z() - vtx.z()
                        );
                double distanceMag      = ROOT::Math::Mag(distanceVector);
                double distanceSigma    = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/distanceMag;
                double distanceSigValue = distanceMag/distanceSigma;

                double misIDMass_Om_pila = -999;
                double misIDMass_Om_lapi = -999;
                if(v0IDName_ == "Omega")
                {
                    double pd1 = d1->p();
                    double pd2 = d2->p();
                    TVector3 dauvec1(d1->px(),d1->py(),d1->pz());
                    TVector3 dauvec2(d2->px(),d2->py(),d2->pz());
                    TVector3 dauvecsum(dauvec1+dauvec2);
                    double massd1=piMass;
                    double massd2=lambdaMass;
                    double energyd1 = sqrt(massd1*massd1+pd1*pd1);
                    double energyd2 = sqrt(massd2*massd2+pd2*pd2);
                    double invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
                    misIDMass_Om_pila = invmass - xiMass;
                    if(fabs(misIDMass_Om_pila) < misIDMassCut_) continue;

                    massd1=lambdaMass;
                    massd2=piMass;
                    energyd1 = sqrt(massd1*massd1+pd1*pd1);
                    energyd2 = sqrt(massd2*massd2+pd2*pd2);
                    invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
                    misIDMass_Om_lapi = invmass - xiMass;
                    if(fabs(misIDMass_Om_lapi) < misIDMassCut_) continue;
                }

                if(v0IDName_ == "Xi")
                {
                    XI.cas3DIpSigValue_    = cas3DIpSigValue;
                    XI.casPi3DIpSigValue_  = casBat3DIpSigValue;
                    XI.VTrkPi3DIpSigValue_ = VTrkPi3DIpSigValue;
                    XI.VTrkP3DIpSigValue_  = VTrkP3DIpSigValue;
                    XI.casFlightSigValue_  = casFlightSigValue;
                    XI.distanceSigValue_   = distanceSigValue;
                    XI.mass_               = mass;
                    XI.pt_                 = pt;
                    XI.rapidity_           = rap;
                    XI.eta_                = eta;

                    XiTree->Fill();
                }
                else if(v0IDName_ == "Omega")
                {
                    OM.cas3DIpSigValue_    = cas3DIpSigValue;
                    OM.casPi3DIpSigValue_  = casBat3DIpSigValue;
                    OM.VTrkPi3DIpSigValue_ = VTrkPi3DIpSigValue;
                    OM.VTrkP3DIpSigValue_  = VTrkP3DIpSigValue;
                    OM.casFlightSigValue_  = casFlightSigValue;
                    OM.distanceSigValue_   = distanceSigValue;
                    OM.mass_               = mass;
                    OM.pt_                 = pt;
                    OM.rapidity_           = rap;
                    OM.eta_                = eta;

                    OmTree->Fill();
                }
                else cout << "unexpected particle option" << endl;
            }
        }
    }

    if(v0IDName_ == "Kshort" || v0IDName_ == "Lambda")
    {
        if(EtaPtCutnTracks >= multLow_ && EtaPtCutnTracks < multHigh_)
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

                //if(doRap_)
                //{
                    //if(rapidity_v0 > rapMax_ || rapidity_v0 < rapMin_) continue;
                //}

                secvz = v0cand->vz(); secvx = v0cand->vx(); secvy = v0cand->vy();

                //trkNHits
                int nhit1 = dau1->numberOfValidHits();
                int nhit2 = dau2->numberOfValidHits();

                if(nhit1 <= nHitCut1_ || nhit2 <= nHitCut2_) continue;

                double pt1 = dau1->pt();
                double pt2 = dau2->pt();

                if(pt1 <= ptCut1_ || pt2 <= ptCut2_) continue;

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
                       //double charged1 = dau1->charge();
                double pd2 = d2->p();
                       //double charged2 = dau2->charge();

                TVector3 dauvec1(d1->px(),d1->py(),d1->pz());
                TVector3 dauvec2(d2->px(),d2->py(),d2->pz());
                TVector3 dauvecsum(dauvec1+dauvec2);

                double energyd1e = sqrt(electronMassSquared+pd1*pd1);
                double energyd2e = sqrt(electronMassSquared+pd2*pd2);
                double invmass_ee = sqrt((energyd1e+energyd2e)*(energyd1e+energyd2e)-dauvecsum.Mag2());
                if(invmass_ee<misIDMassCutEE_) continue;

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
                    if(fabs(misIDMass_la)<misIDMassCut_) continue;
                }

                if(v0IDName_ == "Kshort")
                {
                    double massd1=piMass;
                    double massd2=protonMass;
                    double energyd1 = sqrt(massd1*massd1+pd1*pd1);
                    double energyd2 = sqrt(massd2*massd2+pd2*pd2);
                    double invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
                    misIDMass_ks_pip = invmass - lambdaMass;
                    if(fabs(misIDMass_ks_pip)<misIDMassCut_) continue;

                    massd2=piMass;
                    massd1=protonMass;
                    energyd1 = sqrt(massd1*massd1+pd1*pd1);
                    energyd2 = sqrt(massd2*massd2+pd2*pd2);
                    invmass = sqrt((energyd1+energyd2)*(energyd1+energyd2)-dauvecsum.Mag2());
                    misIDMass_ks_ppi = invmass - lambdaMass;
                    if(fabs(misIDMass_ks_ppi)<misIDMassCut_) continue;
                }

                if(v0IDName_ == "Kshort")
                {
                    KS.eta_               = eta_v0;
                    KS.mass_              = mass_v0;
                    KS.pt_                = pt_v0;
                    KS.rapidity_          = rapidity_v0;
                    KS.dzSig1_            = dzos1;
                    KS.dzSig2_            = dzos2;
                    KS.dxySig1_           = dxyos1;
                    KS.dxySig2_           = dxyos2;
                    KS.vtxChi2_           = vtxChi2;
                    KS.cosTheta_          = agl;
                    KS.decayLSig_         = dlos;

                    KsTree->Fill();
                }

                if(v0IDName_ == "Lambda")
                {
                    LA.eta_              = eta_v0;
                    LA.mass_             = mass_v0;
                    LA.pt_               = pt_v0;
                    LA.rapidity_         = rapidity_v0;
                    LA.dzSig1_           = dzos1;
                    LA.dzSig2_           = dzos2;
                    LA.dxySig1_          = dxyos1;
                    LA.dxySig2_          = dxyos2;
                    LA.vtxChi2_          = vtxChi2;
                    LA.cosTheta_         = agl;
                    LA.decayLSig_        = dlos;

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
    if(v0IDName_ == "Xi") cout << "Will Access Xi" << endl;
    if(v0IDName_ == "Omega") cout << "Will Access Om" << endl;
    if(v0IDName_ == "Kshort") cout << "Will Access Ks" << endl;
    if(v0IDName_ == "Lambda") cout << "Will Access La" << endl;

    edm::Service<TFileService> fs;
    if(v0IDName_ == "Xi")
    {
        XiTree = fs->make<TTree>("XiTree","XiCutParameters");
        XiTree->Branch("xi3dipsig",     &XI.cas3DIpSigValue_,"xi3dipsig/F");
        XiTree->Branch("xipi3dipsig",   &XI.casPi3DIpSigValue_,"xipi3dipsig/F");
        XiTree->Branch("vtrkpi3dipsig", &XI.VTrkPi3DIpSigValue_,"vtrkpi3dipsig/F");
        XiTree->Branch("vtrkp3dipsig",  &XI.VTrkP3DIpSigValue_,"vtrkp3dipsigpt/F");
        XiTree->Branch("xiflightsig",   &XI.casFlightSigValue_,"xiflightsig/F");
        XiTree->Branch("distancesig",   &XI.distanceSigValue_,"distancesig/F");
        XiTree->Branch("mass",          &XI.mass_,"mass/F");
        XiTree->Branch("rapidity",      &XI.rapidity_,"rapidity/F");
        XiTree->Branch("eta",           &XI.eta_,"eta/F");
        XiTree->Branch("pt",            &XI.pt_,"pt/F");
    }

    if(v0IDName_ == "Omega")
    {
        OmTree = fs->make<TTree>("OmTree","OmCutParameters");
        OmTree->Branch("om3dipsig",     &OM.cas3DIpSigValue_,"om3dipsig/F");
        OmTree->Branch("omKaon3dipsig", &OM.casPi3DIpSigValue_,"omKaon3dipsig/F");
        OmTree->Branch("vtrkpi3dipsig", &OM.VTrkPi3DIpSigValue_,"vtrkpi3dipsig/F");
        OmTree->Branch("vtrkp3dipsig",  &OM.VTrkP3DIpSigValue_,"vtrkp3dipsigpt/F");
        OmTree->Branch("omflightsig",   &OM.casFlightSigValue_,"omflightsig/F");
        OmTree->Branch("distancesig",   &OM.distanceSigValue_,"distancesig/F");
        OmTree->Branch("mass",          &OM.mass_,"mass/F");
        OmTree->Branch("rapidity",      &OM.rapidity_,"rapidity/F");
        OmTree->Branch("eta",           &OM.eta_,"eta/F");
        OmTree->Branch("pt",            &OM.pt_,"pt/F");
    }


    if(v0IDName_ == "Kshort")
    {
        KsTree = fs->make<TTree>("KsTree","KsCutParameters");
        KsTree->Branch("eta",          &KS.eta_,"eta/F");
        KsTree->Branch("mass",         &KS.mass_,"mass/F");
        KsTree->Branch("pt",           &KS.pt_,"pt/F");
        KsTree->Branch("rapidity",     &KS.rapidity_,"rapidity/F");
        KsTree->Branch("dzSig1",       &KS.dzSig1_,"dzSig1/F");
        KsTree->Branch("dzSig2",       &KS.dzSig2_,"dzSig2/F");
        KsTree->Branch("dxySig1",      &KS.dxySig1_,"dxySig1/F");
        KsTree->Branch("dxySig2",      &KS.dxySig2_,"dxySig2/F");
        KsTree->Branch("vtxChi2",      &KS.vtxChi2_,"vtxChi2/F");
        KsTree->Branch("cosTheta",     &KS.cosTheta_,"cosTheta/F");
        KsTree->Branch("decayLSig",    &KS.decayLSig_,"decayLSig/F");
    }

    if(v0IDName_ == "Lambda")
    {
        LaTree = fs->make<TTree>("LaTree","LaCutParameters");
        LaTree->Branch("eta",         &LA.eta_,"eta/F");
        LaTree->Branch("mass",        &LA.mass_,"mass/F");
        LaTree->Branch("pt",          &LA.pt_,"pt/F");
        LaTree->Branch("rapidity",    &LA.rapidity_,"rapidity/F");
        LaTree->Branch("dzSig1",      &LA.dzSig1_,"dzSig1/F");
        LaTree->Branch("dzSig2",      &LA.dzSig2_,"dzSig2/F");
        LaTree->Branch("dxySig1",     &LA.dxySig1_,"dxySig1/F");
        LaTree->Branch("dxySig2",     &LA.dxySig2_,"dxySig2/F");
        LaTree->Branch("vtxChi2",     &LA.vtxChi2_,"vtxChi2/F");
        LaTree->Branch("cosTheta",    &LA.cosTheta_,"cosTheta/F");
        LaTree->Branch("decayLSig",   &LA.decayLSig_,"decayLSig/F");
    }
    nEv = fs->make<TH1D>("nEv","nEv",10,0,10);
}

// ------------ method called once each job just after ending the event
//loop  ------------
void
V0XiOmTTreeProducer::endJob(){}

DEFINE_FWK_MODULE(V0XiOmTTreeProducer);
