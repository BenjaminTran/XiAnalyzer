// -*- C++ -*-
//
// Package:    MatchSelector
// Class:      MatchSelector
//
/**\class MatchSelector MatchSelector.cc RiceHIG/V0Analysis/src/MatchSelector.cc

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
MatchSelector::MatchSelector(const edm::ParameterSet& iConfig)
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
  _gnV0Collection = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gnV0Collection"));
  //_gnCollection = consumes<reco::GenParticleCollection>(edm::InputTag("gnCollection"));
  // Trying this with Candidates instead of the simple reco::Vertex
  produces< reco::VertexCompositeCandidateCollection >(v0IDName_);

}

// (Empty) Destructor
MatchSelector::~MatchSelector() {

}

//
// Methods
//

// Producer Method
void MatchSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    //   using std::vector;
    using namespace edm;
    //   using namespace reco;

    // select on requirement of valid vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(_vertexCollName,vertices);
    //double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
    //double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    //const reco::Vertex & vtx = (*vertices)[0];
    //bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    //bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    //   if(bestvz < -15.0 || bestvz>15.0) return;

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
    iEvent.getByToken(_V0Collection, v0candidates);
    if(!v0candidates.isValid())
    {
        std::cout << "v0Collection invalid" << std::endl;
        return;
    }

    //edm::Handle<reco::GenParticleCollection> gencand;
    //iEvent.getByToken(_gnCollection, gencand);
    //if(!gencand.isValid()) return;
    edm::Handle<reco::GenParticleCollection> gencand;
    iEvent.getByToken(_gnV0Collection, gencand);
    if(!gencand.isValid())
    {
        std::cout << "genV0Collection invalid" << std::endl;
        return;
    }

    // Create auto_ptr for each collection to be stored in the Event
    std::auto_ptr< reco::VertexCompositeCandidateCollection >
        theNewV0Cands( new reco::VertexCompositeCandidateCollection() );

    for( reco::VertexCompositeCandidateCollection::const_iterator v0cand = v0candidates->begin();
            v0cand != v0candidates->end();
            v0cand++) {

        //double secvz=-999.9, secvx=-999.9, secvy=-999.9;

        double pt_V0 = v0cand->pt();
        double phi_V0 = v0cand->phi();
        double eta_V0 = v0cand->eta();
        double rap_V0 = v0cand->rapidity();

        if(doRap_)
        {
            if(rap_V0 > rapMax_ || rap_V0 < rapMin_) continue;
        }
        else
        {
            if(eta_V0 > etaCutMax_ || eta_V0 < etaCutMin_) continue;
        }

        //secvz = v0cand->vz(); secvx = v0cand->vx(); secvy = v0cand->vy();


        if(pt_V0 <= ptCut1_) continue;


        for(reco::GenParticleCollection::const_iterator gncand = gencand->begin(); gncand != gencand->end(); gncand++)
        {
            double eta_gen = gncand->eta();
            double phi_gen = gncand->phi();
            double pt_gen  = gncand->pt();
            int id = gncand->pdgId();
            if(v0IDName_ == "Kshort")
            {
                if(fabs(id) == 310)
                {
                    double dphi = phi_gen - phi_V0;
                    double deta = eta_gen - eta_V0;
                    double dR = sqrt(dphi*dphi + deta*deta);
                    double dpt = pt_gen - pt_V0;
                    dR_ks->Fill(dR);
                    dpt_ks->Fill(dpt/pt_gen);
                    if(dR < 0.5 && dpt/pt_gen < 0.5)
                    {
                        theNewV0Cands->push_back(*v0cand);
                        continue;
                    }
                }
            }

            if(v0IDName_ == "Lambda")
            {
                if(fabs(id) == 3122)
                {
                    double dphi = phi_gen - phi_V0;
                    double deta = eta_gen - eta_V0;
                    double dR = sqrt(dphi*dphi + deta*deta);
                    double dpt = pt_gen - pt_V0;
                    dR_la->Fill(dR);
                    dpt_la->Fill(dpt/pt_gen);
                    if(dR < 0.5 && dpt/pt_gen < 0.5)
                    {
                        theNewV0Cands->push_back(*v0cand);
                        continue;
                    }
                }
            }
        }
    }
    // Write the collections to the Event
    iEvent.put( theNewV0Cands, std::string(v0IDName_) );
    //   iEvent.put( theNewV0Cands, std::string("") );
}


void MatchSelector::beginJob() {

    edm::Service<TFileService> fs;

    TH1D::SetDefaultSumw2();

    dR_ks = fs->make<TH1D>("dR_ks","dR_ks",100,0,1);
    dR_la = fs->make<TH1D>("dR_la","dR_la",100,0,1);

    dpt_ks = fs->make<TH1D>("dpt_ks","dpt_ks",100,0,1);
    dpt_la = fs->make<TH1D>("dpt_la","dpt_la",100,0,1);
}


void MatchSelector::endJob() {
}
