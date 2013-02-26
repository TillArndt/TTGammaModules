// -*- C++ -*-
//
// Package:    TTGammaMerger
// Class:      TTGammaMerger
// 
/**\class TTGammaMerger TTGammaMerger.cc MyPackage/TTGammaMerger/src/TTGammaMerger.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Heiner Tholen
//         Created:  Wed May 23 20:38:31 CEST 2012
// $Id: TTGammaMerger.cc,v 1.4 2012/11/26 10:10:30 htholen Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <iostream>
#include <TH1D.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//

class TTGammaMerger : public edm::EDFilter {
   public:
      explicit TTGammaMerger(const edm::ParameterSet&);
      ~TTGammaMerger();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      const double ptCut_;
      const double drCut_;
      const double legPtCut_;
      const bool is2to7_;
      TH1D *kickedPhotons_;
      TH1D *survivingPhotons_;
      TH1D *allPhotons_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TTGammaMerger::TTGammaMerger(const edm::ParameterSet& iConfig) :
    ptCut_(iConfig.getParameter<double>("ptCut")),
    drCut_(iConfig.getParameter<double>("drCut")),
    legPtCut_(iConfig.getUntrackedParameter<double>("legPtCut", 0.)),
    is2to7_(iConfig.getUntrackedParameter<bool>("is2to7", false))
{
    edm::Service<TFileService> fs;
    kickedPhotons_      = fs->make<TH1D>("kickedPhotons",    ";photon e_{T} / GeV;number of photons", 70, 0., 700.);
    survivingPhotons_   = fs->make<TH1D>("survivingPhotons", ";photon e_{T} / GeV;number of photons", 70, 0., 700.);
    allPhotons_         = fs->make<TH1D>("allPhotons",       ";photon e_{T} / GeV;number of photons", 70, 0., 700.);
}


TTGammaMerger::~TTGammaMerger()
{
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TTGammaMerger::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    using reco::GenParticle;
    using reco::deltaR;

    // Why all this?
    // I want to sort out events, that have been simulated with ttgamma matrix element.

    Handle<vector<reco::GenParticle> > genParticles;
    iEvent.getByLabel(InputTag("genParticles"), genParticles);

    Handle<TtGenEvent> ttGenEvent;
    iEvent.getByLabel(InputTag("genEvt"), ttGenEvent);

    // if not semimuonic, this is surely not simulated in ttgamma 2 to 7 ME
    if (!ttGenEvent->isTtBar()) return true;
    if (!ttGenEvent->isSemiLeptonic(WDecay::kMuon)) return true;


    // find legs and all relevant particles
    vector<const GenParticle*> legs;
    vector<const GenParticle*> all;
    const GenParticle* tlep = ttGenEvent->leptonicDecayTop();
    const GenParticle* thad = ttGenEvent->hadronicDecayTop();
    if (is2to7_) {
        const GenParticle* gp = ttGenEvent->lepton();
        if (!gp) gp = ttGenEvent->leptonBar();
        legs.push_back(gp);
        all.push_back(gp);
        gp = ttGenEvent->leptonicDecayB();
        legs.push_back(gp);
        all.push_back(gp);
        gp = ttGenEvent->hadronicDecayB();
        legs.push_back(gp);
        all.push_back(gp);
        gp = ttGenEvent->hadronicDecayQuark();
        legs.push_back(gp);
        all.push_back(gp);
        gp = ttGenEvent->hadronicDecayQuarkBar();
        legs.push_back(gp);
        all.push_back(gp);
        all.push_back(ttGenEvent->leptonicDecayW());
        all.push_back(ttGenEvent->hadronicDecayW());
    } else {
        legs.push_back(tlep);
        legs.push_back(thad);
    }
    all.push_back(tlep);
    all.push_back(thad);
    for (unsigned i = 0; i < tlep->numberOfMothers(); ++i)
        all.push_back((GenParticle*)tlep->mother(i));
    for (unsigned i = 0; i < thad->numberOfMothers(); ++i)
        all.push_back((GenParticle*) thad->mother(i));


    // check legs pt cut (which is 0. by default)
    if (legPtCut_ > 1e-43 && is2to7_) {
        for (unsigned i = 0; i < legs.size(); ++i) {
            if (legs.at(i)->pt() < legPtCut_)
                return true;
        }
    }

    // find relevant photons
    vector<const GenParticle*> photons;
    for (unsigned i = 0; i < all.size(); ++i) {
        for (unsigned j = 0; j < all.at(i)->numberOfDaughters(); ++j) {
            const GenParticle* daughter = (const GenParticle*) all.at(i)->daughter(j);
            if (daughter->pdgId()*daughter->pdgId() == 22*22) {
                 photons.push_back(daughter);
                 allPhotons_->Fill(daughter->et());
            }
        }
    }

    // sort out fails (must fulfill both cuts)
    bool foundNoDrUnderCut = true;
    for (unsigned i = 0; i < photons.size(); ++i) {
        const GenParticle* photon = photons.at(i);
        if (photon->pt() > ptCut_) {
            for (unsigned j = 0; j < legs.size(); ++j) {
                const GenParticle* leg = legs.at(j);
                if (deltaR(*photon, *leg) < drCut_) {
                    foundNoDrUnderCut = false;
                    cout << "<TTGammaMerger>: removing Event! "
                         << "Photon pt < ptCut: (" << photon->pt() << " > " << ptCut_
                         << ") and no deltaR to a leg smaller than " << drCut_ << endl;
                    break;
                }
            }
        }
    }

    if (photons.size() && foundNoDrUnderCut) {
        for (unsigned i = 0; i < photons.size(); ++i)
            kickedPhotons_->Fill(photons.at(i)->et());
        return false;
    } else {
        for (unsigned i = 0; i < photons.size(); ++i)
            survivingPhotons_->Fill(photons.at(i)->et());
    }
    return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
TTGammaMerger::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TTGammaMerger::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
TTGammaMerger::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
TTGammaMerger::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
TTGammaMerger::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
TTGammaMerger::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTGammaMerger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TTGammaMerger);
