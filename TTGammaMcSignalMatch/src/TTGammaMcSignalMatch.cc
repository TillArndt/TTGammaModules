// -*- C++ -*-
//
// Package:    TTGammaMcSignalMatch
// Class:      TTGammaMcSignalMatch
// 
/**\class TTGammaMcSignalMatch TTGammaMcSignalMatch.cc TTGammaModules/TTGammaMcSignalMatch/src/TTGammaMcSignalMatch.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Heiner Tholen
//         Created:  Wed Mar 13 11:56:59 CET 2013
// $Id: TTGammaMcSignalMatch.cc,v 1.1 2013/03/27 12:52:01 htholen Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <DataFormats/PatCandidates/interface/Photon.h>
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//

class TTGammaMcSignalMatch : public edm::EDProducer {
   public:
      explicit TTGammaMcSignalMatch(const edm::ParameterSet&);
      ~TTGammaMcSignalMatch();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
   edm::InputTag patPhotons_;
   edm::InputTag genSignal_;
};

//
// constructors and destructor
//
TTGammaMcSignalMatch::TTGammaMcSignalMatch(const edm::ParameterSet& iConfig):
    patPhotons_(iConfig.getParameter<edm::InputTag>("patPhotons")),
    genSignal_(iConfig.getParameter<edm::InputTag>("genSignal"))
{
    produces<std::vector<pat::Photon> >();
}


TTGammaMcSignalMatch::~TTGammaMcSignalMatch()
{
}


// ------------ method called to produce the data  ------------
void
TTGammaMcSignalMatch::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    using reco::GenParticle;

    edm::Handle<std::vector<pat::Photon> > patPhotons;
    iEvent.getByLabel(patPhotons_, patPhotons);

    edm::Handle<std::vector<GenParticle> > genSignal;
    iEvent.getByLabel(genSignal_, genSignal);

    // loop over photons and find daugthers of particles in all
    std::vector<pat::Photon>* signalPhotons = new std::vector<pat::Photon>();
    std::vector<pat::Photon>::const_iterator patPhoton = patPhotons->begin();
    for ( ; patPhoton != patPhotons->end(); ++patPhoton) {
        GenParticle* gp = (GenParticle*) patPhoton->genParticle();
        if (!gp) continue;
        // set gp to uppermost photon
        while(abs(gp->mother()->pdgId()) == 22) {
            gp = (GenParticle*) gp->mother();
        }
        std::vector<GenParticle>::const_iterator signalPho = genSignal->begin();
        for (; signalPho != genSignal->end(); ++signalPho) {
            if ((signalPho->p4() - gp->p4()).M() < 1e-43) {// consider same, when spatial momentum difference vanishes
                signalPhotons->push_back(*patPhoton);
                break;
            }
        }
    }

   std::auto_ptr<std::vector<pat::Photon> > pOut(signalPhotons);
   iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void 
TTGammaMcSignalMatch::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TTGammaMcSignalMatch::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
TTGammaMcSignalMatch::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TTGammaMcSignalMatch::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TTGammaMcSignalMatch::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TTGammaMcSignalMatch::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTGammaMcSignalMatch::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TTGammaMcSignalMatch);
