// -*- C++ -*-
//
// Package:    GenSelector
// Class:      GenSelector
// 
/**\class GenSelector GenSelector.cc RiceHIG/V0Analysis/src/GenSelector.cc

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

#include "XiAnalyzer/XiAnalyzer/interface/GenSelector.h"

// Constructor
GenSelector::GenSelector(const edm::ParameterSet& iConfig)
{
    using std::string;

    v0IDName_     = iConfig.getParameter<string>("v0IDName");
    _gnCollection = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gnCollection"));
    // Trying this with Candidates instead of the simple reco::Vertex
    produces< reco::GenParticleCollection >(v0IDName_);

}

// (Empty) Destructor
GenSelector::~GenSelector() {

}

//
// Methods
//

// Producer Method
void GenSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

   using namespace edm;

   edm::Handle<reco::GenParticleCollection> gencand;
   iEvent.getByToken(_gnCollection, gencand);
   if(!gencand.isValid()) 
   {
       std::cout << "Actual GenCollection is invalid" << std::endl;
       return;
   }

   // Create auto_ptr for each collection to be stored in the Event
   std::auto_ptr< reco::GenParticleCollection >
     theGenCands( new reco::GenParticleCollection );


   for( reco::GenParticleCollection::const_iterator gncand = gencand->begin();
         gncand != gencand->end();
         gncand++)
   {
       int id = gncand->pdgId();
       //int st = gncand->status();
       int rap = gncand->rapidity();

       if(v0IDName_ == "Kshort")
       {
           if(fabs(id) == 310 && fabs(rap) < 1.0)
           {
               theGenCands->push_back(*gncand);
               std::cout << "Kshort fill" << std::endl;
           }
       }

       if(v0IDName_ == "Lambda")
       {
           if(fabs(id) == 3122 && fabs(rap) < 1.0)
           {
               int mid = 0;
               if(gncand->numberOfMothers()==1)
               {
                   const reco::Candidate * mom = gncand->mother();
                   mid = mom->pdgId();
                   if(mom->numberOfMothers()==1)
                   {
                       const reco::Candidate * mom1 = mom->mother();
                       mid = mom1->pdgId();
                   }
               }
               if(fabs(mid)==3322 || fabs(mid)==3312 || fabs(mid)==3324 || fabs(mid)==3314 || fabs(mid)==3334) continue;
               theGenCands->push_back(*gncand);
               std::cout << "Lambda fill" << std::endl;
           }
       }
   }

   // Write the collections to the Event
   iEvent.put( theGenCands, std::string(v0IDName_) );
}


void GenSelector::beginJob() {
}


void GenSelector::endJob() {
}
