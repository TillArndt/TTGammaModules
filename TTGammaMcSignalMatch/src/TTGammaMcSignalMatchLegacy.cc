// -*- C++ -*-
//
// Package:    TTGammaMcSignalMatchLegacy
// Class:      TTGammaMcSignalMatchLegacy
// 
/**\class TTGammaMcSignalMatchLegacy TTGammaMcSignalMatchLegacy.cc TTGammaModules/TTGammaMcSignalMatch/src/TTGammaMcSignalMatchLegacy.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Heiner Tholen
//         Created:  Wed Mar 13 11:56:59 CET 2013
// $Id: TTGammaMcSignalMatchLegacy.cc,v 1.1 2013/03/27 12:52:01 htholen Exp $
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

class TTGammaMcSignalMatchLegacy : public edm::EDProducer {
   public:
      explicit TTGammaMcSignalMatchLegacy(const edm::ParameterSet&);
      ~TTGammaMcSignalMatchLegacy();

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
   edm::InputTag photons_;
   const bool is2to5_;
   const bool is2to7_;
};

//
// constructors and destructor
//
TTGammaMcSignalMatchLegacy::TTGammaMcSignalMatchLegacy(const edm::ParameterSet& iConfig):
    photons_(iConfig.getParameter<edm::InputTag>("src")),
    is2to5_(iConfig.getUntrackedParameter<bool>("is2to5", false)),
    is2to7_(iConfig.getUntrackedParameter<bool>("is2to7", false))
{
    assert( !(is2to5_ && is2to7_) );
    produces<std::vector<pat::Photon> >();
}


TTGammaMcSignalMatchLegacy::~TTGammaMcSignalMatchLegacy()
{
}


// ------------ method called to produce the data  ------------
void
TTGammaMcSignalMatchLegacy::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    using reco::GenParticle;

    Handle<TtGenEvent> ttGenEvent;
    iEvent.getByLabel(InputTag("genEvt"), ttGenEvent);

    edm::Handle<std::vector<pat::Photon> > photons;
    iEvent.getByLabel(photons_, photons);

    // get all particles where a photon could be radiated from 
    vector<const GenParticle*> all;
    const GenParticle* top    = ttGenEvent->top();
    const GenParticle* topBar = ttGenEvent->topBar();
    if (is2to7_) {
        // if not semimuonic, this is not simulated in ttgamma 2 to 7 ME
        if (!ttGenEvent->isTtBar()) return;
        if (!ttGenEvent->isSemiLeptonic(WDecay::kMuon)) return;

        const GenParticle* gp = ttGenEvent->lepton();
        if (!gp) gp = ttGenEvent->leptonBar();
        all.push_back(gp);
        gp = ttGenEvent->leptonicDecayB();
        all.push_back(gp);
        gp = ttGenEvent->hadronicDecayB();
        all.push_back(gp);
        gp = ttGenEvent->hadronicDecayQuark();
        all.push_back(gp);
        gp = ttGenEvent->hadronicDecayQuarkBar();
        all.push_back(gp);
        all.push_back(ttGenEvent->leptonicDecayW());
        all.push_back(ttGenEvent->hadronicDecayW());
    } else if (is2to5_) {
        for (unsigned i = 0; i < top->numberOfDaughters(); ++i) {
            const GenParticle* gp = (GenParticle*) top->daughter(i);
            if (abs(gp->pdgId()) < 6 || abs(gp->pdgId()) == 24) {
                all.push_back(gp);
            }
        }
        for (unsigned i = 0; i < topBar->numberOfDaughters(); ++i) {
            const GenParticle* gp = (GenParticle*) topBar->daughter(i);
            if (abs(gp->pdgId()) < 6 || abs(gp->pdgId()) == 24) {
                all.push_back(gp);
            }
        }
    }
    all.push_back(top);
    all.push_back(topBar);
    for (unsigned i = 0; i < top->numberOfMothers(); ++i)
        all.push_back((GenParticle*) top->mother(i));
    for (unsigned i = 0; i < topBar->numberOfMothers(); ++i)
        all.push_back((GenParticle*) topBar->mother(i));

    // get photons radiated from all
    std::vector<GenParticle*> genPhotons;
    for (unsigned i = 0; i < all.size(); ++i) {
        const GenParticle *gp = all[i];
        for (unsigned j = 0; j < gp->numberOfDaughters(); ++j) {
            GenParticle* daughter = (GenParticle*) gp->daughter(j);
            if (daughter->pdgId() == 22) 
                genPhotons.push_back(daughter);
        }
    }

    // loop over photons and find daugthers of particles in all
    std::vector<pat::Photon>* signalPhotons = new std::vector<pat::Photon>();
    std::vector<pat::Photon>::const_iterator photon = photons->begin();
    for ( ; photon != photons->end(); ++photon) {
        GenParticle* gp = (GenParticle*) photon->genParticle();
        if (!gp) continue;
        // set gp to uppermost photon
        while(abs(gp->mother()->pdgId()) == 22) {
            gp = (GenParticle*) gp->mother();
        }
        for (unsigned i = 0; i < genPhotons.size(); ++i) {
            if ((genPhotons[i]->p4() - gp->p4()).M() < 1e-43) {// consider same, when spatial momentum difference vanishes
                signalPhotons->push_back(*photon);
                break;
            }
        }
    }

   std::auto_ptr<std::vector<pat::Photon> > pOut(signalPhotons);
   iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void 
TTGammaMcSignalMatchLegacy::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TTGammaMcSignalMatchLegacy::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
TTGammaMcSignalMatchLegacy::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TTGammaMcSignalMatchLegacy::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TTGammaMcSignalMatchLegacy::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TTGammaMcSignalMatchLegacy::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTGammaMcSignalMatchLegacy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TTGammaMcSignalMatchLegacy);
