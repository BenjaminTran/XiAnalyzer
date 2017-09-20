// -*- C++ -*-
//
// Package:    V0Selector
// Class:      V0Selector
// 
/**\class V0Selector V0Selector.cc RiceHIG/V0Analysis/src/V0Selector.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author: Wei Li 
//
//


// system include files
#include <memory>

#include "XiAnalyzer/XiAnalyzer/interface/MatchSelector.h"

const double piMass = 0.13957018;
const double piMassSquared = piMass*piMass;
const double protonMass = 0.93827203;
const double protonMassSquared = protonMass*protonMass;
const double electronMass = 0.000511;
const double electronMassSquared = electronMass*electronMass;
const double lambdaMass = 1.115683;
const double kshortMass = 0.497614;

// Constructor
V0Selector::V0Selector(const edm::ParameterSet& iConfig)
{
  using std::string;

  //vertexCollName_ = iConfig.getParameter<edm::InputTag>("vertexCollName");
  v0CollName_     = iConfig.getParameter<string>("v0CollName");
  v0IDName_       = iConfig.getParameter<string>("v0IDName");
  etaCutMin_      = iConfig.getParameter<double>("etaCutMin");
  etaCutMax_      = iConfig.getParameter<double>("etaCutMax");
  ptCut1_         = iConfig.getParameter<double>("ptCut1");
  ptCut2_         = iConfig.getParameter<double>("ptCut2");
  doRap_          = iConfig.getParameter<bool>("dorap");
  rapMax_         = iConfig.getParameter<double>("rapMax");
  rapMin_         = iConfig.getParameter<double>("rapMin");
  _vertexCollName = consumes<reco::VertexCollection>( iConfig.getParameter<edm::InputTag>( "vertexCollName" ) );
  _V0Collection = consumes<reco::VertexCompositeCandidateCollection>( edm::InputTag( v0CollName_,v0IDName_,"ANASKIM" ) );
  _gnV0Collection = consumes<reco::GenParticleCollection>(edm::InputTag("gnV0Collection"))
  //_gnCollection = consumes<reco::GenParticleCollection>(edm::InputTag("gnCollection"));
  // Trying this with Candidates instead of the simple reco::Vertex
  produces< reco::VertexCompositeCandidateCollection >(v0IDName_);

}

// (Empty) Destructor
V0Selector::~V0Selector() {

}

//
// Methods
//

// Producer Method
void V0Selector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

//   using std::vector;
   using namespace edm;
//   using namespace reco;
    
    // select on requirement of valid vertex
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(_vertexCollName,vertices);
   double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
   double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
   const reco::Vertex & vtx = (*vertices)[0];
   bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
   bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
//   if(bestvz < -15.0 || bestvz>15.0) return;

   edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
   iEvent.getByToken(_V0Collection, v0candidates);
   if(!v0candidates.isValid()) return;

   //edm::Handle<reco::GenParticleCollection> gencand;
   //iEvent.getByToken(_gnCollection, gencand);
   //if(!gencand.isValid()) return;
   edm::Handle<reco::GenParticleCollection> gencand;
   iEvent.getByToken(_gnV0Collection, gencand);
   if(!gencand.isValid()) return;

   // Create auto_ptr for each collection to be stored in the Event
   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     theNewV0Cands( new reco::VertexCompositeCandidateCollection() );

   for( reco::VertexCompositeCandidateCollection::const_iterator v0cand = v0candidates->begin();
         v0cand != v0candidates->end();
         v0cand++) {

       double secvz=-999.9, secvx=-999.9, secvy=-999.9;

       const reco::Candidate * d1 = v0cand->daughter(0);
       const reco::Candidate * d2 = v0cand->daughter(1);

       reco::TrackRef dau1 = d1->get<reco::TrackRef>();
       reco::TrackRef dau2 = d2->get<reco::TrackRef>();

       //pt,mass
       double eta = v0cand->eta();
//       double pt = v0cand->pt();
       double rap = v0cand->rapidity();
//       double mass = v0cand->mass();

        if(doRap_)
        {
            if(rap > rapMax_ || rap < rapMin_) continue;
        }
        else
        {
            if(eta > etaCutMax_ || eta < etaCutMin_) continue;
        }

       secvz = v0cand->vz(); secvx = v0cand->vx(); secvy = v0cand->vy();

       double pt1 = dau1->pt();
       double pt2 = dau2->pt();

       if(pt1 <= ptCut1_ || pt2 <= ptCut2_) continue;

       double phi1 = dau1->phi();
       double phi2 = dau2->phi();

       double eta1 = dau1->eta();
       double eta2 = dau2->eta();

       for(reco::GenParticleCollection::const_iterator gncand = gencand->begin(); gncand != gencand->end(); gncand++)
       {
           double eta = gncand->eta();
           double phi = gncand->phi();
           double pt  = gncand->pt();
           double mass = gncand->mass();
           int id = gncand->pdgId();
           int st = gncand->status();
           if(v0IDName_ == "Kshort")
           {
               if(fabs(id) == 310)
               {
                   if(gncand->numberOfDaughters()==2)
                   {
                       const reco::Candidate *gen_dau1 = gncand->daughter(0);
                       const reco::Candidate *gen_dau2 = gncand->daughter(1);

                       int id_dau1 = gen_dau1->pdgId();
                       int id_dau2 = gen_dau2->pdgId();

                       double eta_gen_dau1 = gen_dau1->eta();
                       double eta_gen_dau2 = gen_dau2->eta();

                       double phi_gen_dau1 = gen_dau1->phi();
                       double phi_gen_dau2 = gen_dau2->phi();

                       double pt_gen_dau1 = gen_dau1->pt();
                       double pt_gen_dau2 = gen_dau2->pt();
                       if(id_dau1 == 211)
                       {
                           double dphi1 = phi_gen_dau1 - phi1;
                           double deta1 = eta_gen_dau1 - eta1;
                           double dR1 = sqrt(dphi1*dphi1 + deta1*deta1);
                           double dpt1 = pt_gen_dau1 - pt1;
                           double sumpt1 = pt_gen_dau1 + pt1;
                           cout << "Mathc" << endl;
                           //dR Fill daughter1 hist
                           //sumpt Fill daughter1 hist
                           //
                       }
                   }
               }
           }
       }


       theNewV0Cands->push_back( *v0cand );
   }

   // Write the collections to the Event
   iEvent.put( theNewV0Cands, std::string(v0IDName_) );
//   iEvent.put( theNewV0Cands, std::string("") );
}


void V0Selector::beginJob() {
}


void V0Selector::endJob() {
}
