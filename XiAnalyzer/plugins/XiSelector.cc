// -*- C++ -*-
//
// Package:    XiSelector
// Class:      XiSelector
//
/**\class XiSelector XiSelector.cc RiceHIG/V0Analysis/src/XiSelector.cc

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

//#include "RiceHIG/V0Analysis/interface/XiSelector.h"
#include "../interface/XiSelector.h" // For interactive run on lxplus

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
XiSelector::XiSelector(const edm::ParameterSet& iConfig)
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
    cas3DIpSigValue_     = iConfig.getParameter<double>("cas3DIpSigValue");
    casFlightSigValue_   = iConfig.getParameter<double>("casFlightSigValue");
    casBat3DIpSigValue_   = iConfig.getParameter<double>("casBat3DIpSigValue");
    zVertexHigh_        = iConfig.getParameter<double>("zVertexHigh");
    zVertexLow_         = iConfig.getParameter<double>("zVertexLow");
    doRap_              = iConfig.getParameter<bool>("doRap");
    rapMax_             = iConfig.getParameter<double>("rapMax");
    rapMin_             = iConfig.getParameter<double>("rapMin");
    misIDMassCut_ = iConfig.getParameter<double>("misIDMassCut");
    _XiCollection       = consumes<reco::VertexCompositeCandidateCollection>( edm::InputTag(v0CollName_,v0IDName_,"ANASKIM"));
    _vertexCollName     = consumes<reco::VertexCollection>( iConfig.getParameter<edm::InputTag>( "vertexCollName"));

    // The argument gives the instance name of the collection
    produces< reco::VertexCompositeCandidateCollection >(v0IDName_);

}

// Destructor
XiSelector::~XiSelector() {

}

//
// Methods
//

// Producer Method
void XiSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

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
    const reco::Vertex & vtx = (*vertices)[0];
    //bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    //bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();

    // Z Vertex cut
    //if(bestvz > zVertexHigh_ || bestvz < zVertexLow_ ) return;

    edm::Handle< reco::VertexCompositeCandidateCollection > CascadeCandidates;
    iEvent.getByToken(_XiCollection, CascadeCandidates);
    if(!CascadeCandidates.isValid()) return;

    // Create auto_ptr for each collection to be stored in the Event
    std::auto_ptr< reco::VertexCompositeCandidateCollection >
        theNewCasCands( new reco::VertexCompositeCandidateCollection() );


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
        //double mass = CasCand->mass();
        //double pt   = CasCand->pt();
        double rap  = CasCand->rapidity();

        if(doRap_)
        {
        if(rap < rapMin_ || rap > rapMax_) continue;
        }
        else
        {
        if(eta > 2.4 || eta < -2.4) continue;
        }

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

        // Apply Cuts
        if(cas3DIpSigValue     > cas3DIpSigValue_)     continue;
        if(casBat3DIpSigValue   < casBat3DIpSigValue_)   continue;
        if(VTrkPi3DIpSigValue < VTrkPi3DIpSigValue_) continue;
        if(VTrkP3DIpSigValue  < VTrkP3DIpSigValue_)  continue;
        if(casFlightSigValue   < casFlightSigValue_)   continue;
        if(distanceSigValue   < distanceSigValue_)   continue;

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

        theNewCasCands->push_back( *CasCand );

    }

    // Write the collections to the Event
    // Collection is stored as module label : instance
    // e.g. for this the InputTag should be selectV0CandidatesLowXi:Xi
    iEvent.put( theNewCasCands, std::string(v0IDName_) );
}

void XiSelector::beginJob() {

}


void XiSelector::endJob() {

}

