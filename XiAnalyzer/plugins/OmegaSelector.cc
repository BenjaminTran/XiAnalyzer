// -*- C++ -*-
//
// Package:    OmegaSelector
// Class:      OmegaSelector
//
/**\class OmegaSelector OmegaSelector.cc RiceHIG/V0Analysis/src/OmegaSelector.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author: Benjamin Tran
//
//


// system include files
#include <memory>

//#include "RiceHIG/V0Analysis/interface/OmegaSelector.h"
#include "../interface/OmegaSelector.h" // For interactive run on lxplus

using namespace std;
using namespace edm;
using namespace reco;

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
float kaonMass            = 0.493677;
float kaonMass_sigma      = kaonMass*1.e-6;

// Constructor
OmegaSelector::OmegaSelector(const edm::ParameterSet& iConfig)
{
    using std::string;
/*  The following are the cut variables for Impact parameter and decay length
 *  used for Transient Track Method
 *  Impact Parameters
 *  Xi          : xi3DIpSigValue
 *  Lambda_pion : VTrkPi3DIpSigValue
 *  Lambda_prot : VTrkP3DIpSigValue
 *  Xi_pion     : xiPi3DIpSigValue
 *
 *  Decay lengths
 *  Xi     : xiFlightSigValue
 *  Lambda : distanceSigValue
*/

    VTrkP3DIpSigValue_  = iConfig.getParameter<double>("VTrkP3DIpSigValue");
    VTrkPi3DIpSigValue_ = iConfig.getParameter<double>("VTrkPi3DIpSigValue");
    distanceSigValue_   = iConfig.getParameter<double>("distanceSigValue");
    etaCutMax_          = iConfig.getParameter<double>("etaCutMax");
    etaCutMin_          = iConfig.getParameter<double>("etaCutMin");
    nHitCut1_           = iConfig.getParameter<int>("nHitCut1");
    v0CollName_         = iConfig.getParameter<string>("v0CollName");
    v0IDName_           = iConfig.getParameter<string>("v0IDName");
    xi3DIpSigValue_     = iConfig.getParameter<double>("xi3DIpSigValue");
    xiFlightSigValue_   = iConfig.getParameter<double>("xi3DIpSigValue");
    xiPi3DIpSigValue_   = iConfig.getParameter<double>("xiPi3DIpSigValue");
    zVertexHigh_        = iConfig.getParameter<double>("zVertexHigh");
    zVertexLow_         = iConfig.getParameter<double>("zVertexLow");
    dorap_              = iConfig.getParameter<bool>("dorap");
    rapMax_             = iConfig.getParameter<double>("rapMax");
    rapMin_             = iConfig.getParameter<double>("rapMin");
    _OmCollection       = consumes<reco::VertexCompositeCandidateCollection>( edm::InputTag(v0CollName_,v0IDName_,"ANASKIM"));
    _vertexCollName     = consumes<reco::VertexCollection>( iConfig.getParameter<edm::InputTag>( "vertexCollName"));

    // The argument gives the instance name of the collection
    produces< reco::VertexCompositeCandidateCollection >(v0IDName_);

}

// Destructor
OmegaSelector::~OmegaSelector() {

}

//
// Methods
//

// Producer Method
void OmegaSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;
    using namespace reco;
    using std::vector;

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
        //double px_xi   = v0cand->px();
        //double py_xi   = v0cand->py();
        //double pz_xi   = v0cand->pz();

        if(dorap_ && (rap_xi > rapMax_ || rap_xi < rapMin_)) continue;
        else
            if(eta_xi > etaCutMax_ || eta_xi < etaCutMin_) continue;

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


        // Apply Cuts
        if(xi3DIpSigValue     > xi3DIpSigValue_)     continue;
        if(xiPi3DIpSigValue   < xiPi3DIpSigValue_)   continue;
        if(VTrkPi3DIpSigValue < VTrkPi3DIpSigValue_) continue;
        if(VTrkP3DIpSigValue  < VTrkP3DIpSigValue_)  continue;
        if(xiFlightSigValue   < xiFlightSigValue_)   continue;
        if(distanceSigValue   < distanceSigValue_)   continue;

        theNewOmCands->push_back( *v0cand );

        cout<<"Successful"<<endl;
    }

    // Write the collections to the Event
    // Collection is stored as module label : instance
    // e.g. for this the InputTag should be selectV0CandidatesLowXi:Xi
    iEvent.put( theNewOmCands, std::string(v0IDName_) );
}

void OmegaSelector::beginJob() {

}


void OmegaSelector::endJob() {

}

