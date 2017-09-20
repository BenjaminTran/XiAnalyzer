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

        //secvz = v0cand->vz(); secvx = v0cand->vx(); secvy = v0cand->vy();

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
            std::cout << "Status: " << st << std::endl;
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
                        std::cout << "id_dau1: " << id_dau1 << std::endl;
                        std::cout << "id_dau2: " << id_dau2 << std::endl;
                        if(id_dau1 == 211 && gen_dau1->charge() == 1)
                        {
                            double dphi1 = phi_gen_dau1 - phi1;
                            double deta1 = eta_gen_dau1 - eta1;
                            double dR1 = sqrt(dphi1*dphi1 + deta1*deta1);
                            double dpt1 = pt_gen_dau1 - pt1;
                            //dR Fill daughter1 hist
                            //dpt1 Fill daughter1 hist
                            if(dR1 < 0.5 && dpt1/pt_gen_dau1 < 0.5)
                            {
                                double dphi2 = phi_gen_dau2 - phi2;
                                double deta2 = eta_gen_dau2 - eta2;
                                double dR2 = sqrt(dphi2*dphi2 + deta2*deta2);
                                double dpt2 = pt_gen_dau2 - pt2;
                                //dR Fill daughter2 hist
                                //dpt2 Fill daughter2 hist
                                if(dR2 < 0.5 && dpt2/pt_gen_dau2 < 0.5)
                                {
                                    theNewV0Cands->push_back(*v0cand);
                                    continue;
                                }
                            }
                        }
                        else
                        {
                            double dphi1 = phi_gen_dau1 - phi2;
                            double deta1 = eta_gen_dau1 - eta2;
                            double dR1 = sqrt(dphi1*dphi1 + deta1*deta1);
                            double dpt1 = pt_gen_dau1 - pt2;
                            //dR Fill daughter1 hist
                            //dpt1 Fill daughter1 hist
                            if(dR1 < 0.5 && dpt1/pt_gen_dau1 < 0.5)
                            {
                                double dphi2 = phi_gen_dau2 - phi1;
                                double deta2 = eta_gen_dau2 - eta1;
                                double dR2 = sqrt(dphi2*dphi2 + deta2*deta2);
                                double dpt2 = pt_gen_dau2 - pt1;
                                //dR Fill daughter2 hist
                                //dpt2 Fill daughter2 hist
                                if(dR2 < 0.5 && dpt2/pt_gen_dau2 < 0.5)
                                {
                                    theNewV0Cands->push_back(*v0cand);
                                    continue;
                                }
                            }
                        }
                    }
                }
            }

            if(v0IDName_ == "Lambda")
            {
                if(fabs(id) == 3122)
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
                        std::cout << "id_dau1: " << id_dau1 << std::endl;
                        std::cout << "id_dau2: " << id_dau2 << std::endl;
                        if(id_dau1 == 2212 && gen_dau1->charge() == 1)
                        {
                            double dphi1 = phi_gen_dau1 - phi1;
                            double deta1 = eta_gen_dau1 - eta1;
                            double dR1 = sqrt(dphi1*dphi1 + deta1*deta1);
                            double dpt1 = pt_gen_dau1 - pt1;
                            //dR Fill daughter1 hist
                            //dpt1 Fill daughter1 hist
                            if(dR1 < 0.5 && dpt1/pt_gen_dau1 < 0.5)
                            {
                                double dphi2 = phi_gen_dau2 - phi2;
                                double deta2 = eta_gen_dau2 - eta2;
                                double dR2 = sqrt(dphi2*dphi2 + deta2*deta2);
                                double dpt2 = pt_gen_dau2 - pt2;
                                //dR Fill daughter2 hist
                                //dpt2 Fill daughter2 hist
                                if(dR2 < 0.5 && dpt2/pt_gen_dau2 < 0.5)
                                {
                                    theNewV0Cands->push_back(*v0cand);
                                    continue;
                                }
                            }
                        }
                        else
                        {
                            double dphi1 = phi_gen_dau1 - phi2;
                            double deta1 = eta_gen_dau1 - eta2;
                            double dR1 = sqrt(dphi1*dphi1 + deta1*deta1);
                            double dpt1 = pt_gen_dau1 - pt2;
                            //dR Fill daughter1 hist
                            //dpt1 Fill daughter1 hist
                            if(dR1 < 0.5 && dpt1/pt_gen_dau1 < 0.5)
                            {
                                double dphi2 = phi_gen_dau2 - phi1;
                                double deta2 = eta_gen_dau2 - eta1;
                                double dR2 = sqrt(dphi2*dphi2 + deta2*deta2);
                                double dpt2 = pt_gen_dau2 - pt1;
                                //dR Fill daughter2 hist
                                //dpt2 Fill daughter2 hist
                                if(dR2 < 0.5 && dpt2/pt_gen_dau2 < 0.5)
                                {
                                    theNewV0Cands->push_back(*v0cand);
                                    continue;
                                }
                            }
                        }
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
}


void MatchSelector::endJob() {
}
