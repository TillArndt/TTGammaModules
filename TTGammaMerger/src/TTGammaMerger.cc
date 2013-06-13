// -*- C++ -*-
//
// Package:    TTGammaMerger
// Class:      TTGammaMerger
// 
/**\class TTGammaMerger TTGammaMerger.cc MyPackage/TTGammaMerger/src/TTGammaMerger.cc

 Description: Removes Signal Overlap in TTbar sample

 Implementation:
// Sort out events, that have been simulated with ttgamma matrix element.

*/
//
// Original Author:  Heiner Tholen
//         Created:  Wed May 23 20:38:31 CEST 2012
// $Id: TTGammaMerger.cc,v 1.4 2013/06/13 12:05:02 htholen Exp $
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
      const bool is2to5_;
      const bool is2to7_;
      TH1D *etKickedPhotons_;
      TH1D *etSurvivingPhotons_;
      TH1D *etAllPhotons_;
      TH1D *etaKickedPhotons_;
      TH1D *etaSurvivingPhotons_;
      TH1D *etaAllPhotons_;
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
    is2to5_(iConfig.getUntrackedParameter<bool>("is2to5", false)),
    is2to7_(iConfig.getUntrackedParameter<bool>("is2to7", false))
{
    assert( !(is2to5_ && is2to7_) );
    edm::Service<TFileService> fs;
    etaKickedPhotons_     = fs->make<TH1D>("etaKickedPhotons",    ";photon E_{T} / GeV;number of photons", 80, -4., 4.);
    etaSurvivingPhotons_  = fs->make<TH1D>("etaSurvivingPhotons", ";photon E_{T} / GeV;number of photons", 80, -4., 4.);
    etaAllPhotons_        = fs->make<TH1D>("etaAllPhotons",       ";photon E_{T} / GeV;number of photons", 80, -4., 4.);
    etKickedPhotons_      = fs->make<TH1D>("etKickedPhotons",     ";photon E_{T} / GeV;number of photons", 70, 0., 700.);
    etSurvivingPhotons_   = fs->make<TH1D>("etSurvivingPhotons",  ";photon E_{T} / GeV;number of photons", 70, 0., 700.);
    etAllPhotons_         = fs->make<TH1D>("etAllPhotons",        ";photon E_{T} / GeV;number of photons", 70, 0., 700.);
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

    Handle<vector<reco::GenParticle> > genParticles;
    iEvent.getByLabel(InputTag("genParticles"), genParticles);

    Handle<TtGenEvent> ttGenEvent;
    iEvent.getByLabel(InputTag("genEvt"), ttGenEvent);

    // find legs and all relevant particles
    vector<const GenParticle*> legs;
    vector<const GenParticle*> all;
    const GenParticle* top    = ttGenEvent->top();
    const GenParticle* topBar = ttGenEvent->topBar();
    if (is2to7_) {
        // if not semimuonic, this is not simulated in ttgamma 2 to 7 ME
        if (!ttGenEvent->isTtBar()) return true;
        if (!ttGenEvent->isSemiLeptonic(WDecay::kMuon)) return true;

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
    } else if (is2to5_) {
        for (unsigned i = 0; i < top->numberOfDaughters(); ++i) {
            const GenParticle* gp = (GenParticle*) top->daughter(i);
            if (abs(gp->pdgId()) < 6 || abs(gp->pdgId()) == 24) {
                all.push_back(gp);
            }
            if (abs(gp->pdgId()) < 6) {
                legs.push_back(gp);
            }
        }
        for (unsigned i = 0; i < topBar->numberOfDaughters(); ++i) {
            const GenParticle* gp = (GenParticle*) topBar->daughter(i);
            if (abs(gp->pdgId()) < 6 || abs(gp->pdgId()) == 24) {
                all.push_back(gp);
            }
            if (abs(gp->pdgId()) < 6) {
                legs.push_back(gp);
            }
        }
    } else { // 2 to 3
        legs.push_back(top);
        legs.push_back(topBar);
    }
    all.push_back(top);
    all.push_back(topBar);
    for (unsigned i = 0; i < top->numberOfMothers(); ++i)
        all.push_back((GenParticle*) top->mother(i));
    for (unsigned i = 0; i < topBar->numberOfMothers(); ++i)
        all.push_back((GenParticle*) topBar->mother(i));

    // check legs pt cut (which is 0. by default)
    if (legPtCut_ > 1e-43 && (is2to7_ || is2to5_)) {
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
                 etAllPhotons_->Fill(daughter->et());
                 etaAllPhotons_->Fill(daughter->eta());
            }
        }
    }

    // sort out fails (must fulfill both cuts)
    bool foundNoSignalPhoton = true;
    for (unsigned i = 0; i < photons.size(); ++i) {
        const GenParticle* photon = photons.at(i);
        if (photon->pt() > ptCut_) {
            for (unsigned j = 0; j < legs.size(); ++j) {
                const GenParticle* leg = legs.at(j);
                if (deltaR(*photon, *leg) > drCut_) {
                    foundNoSignalPhoton = false;
                    cout << "<TTGammaMerger>: removing Event! "
                         << "Photon pt < ptCut: (" << photon->pt() << " > " << ptCut_
                         << ") and no deltaR to a leg smaller than " << drCut_ << endl;
                    break;
                }
            }
        }
    }

    if (foundNoSignalPhoton) {
        for (unsigned i = 0; i < photons.size(); ++i) {
            etSurvivingPhotons_->Fill(photons.at(i)->et());
            etaSurvivingPhotons_->Fill(photons.at(i)->eta());
        }
     } else {
        for (unsigned i = 0; i < photons.size(); ++i) {
            etKickedPhotons_->Fill(photons.at(i)->et());
            etaKickedPhotons_->Fill(photons.at(i)->eta());
        }
    }
    return foundNoSignalPhoton;
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
