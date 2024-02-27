// -*- C++ -*-
//
// Package:    PlayGround/EGammaNtuples
// Class:      EGammaNtuples
//
/**\class EGammaNtuples EGammaNtuples.cc PlayGround/EGammaNtuples/plugins/EGammaNtuples.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Matthewman
//         Created:  Tue, 27 Feb 2024 14:59:00 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/EgammaObject.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class EGammaNtuples : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit EGammaNtuples(const edm::ParameterSet&);
  ~EGammaNtuples() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleToken_;
  edm::EDGetTokenT<std::vector<trigger::EgammaObject>> eGammaObjectToken_;
  edm::EDGetTokenT<std::vector<trigger::EgammaObject>> eGammaObjectUnseededToken_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>> scBarrelL1SeededToken_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>> scHGCalL1SeededToken_;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
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
EGammaNtuples::EGammaNtuples(const edm::ParameterSet& iConfig)
  : genParticleToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles"))),
  eGammaObjectToken_(consumes<std::vector<trigger::EgammaObject>>(iConfig.getUntrackedParameter<edm::InputTag>("eGammaObjects"))),
  eGammaObjectUnseededToken_(consumes<std::vector<trigger::EgammaObject>>(iConfig.getUntrackedParameter<edm::InputTag>("eGammaObjectsUnSeeded"))),
  scBarrelL1SeededToken_(consumes<std::vector<reco::SuperCluster>>(iConfig.getUntrackedParameter<edm::InputTag>("scBarrelL1Seeded"))),
  scHGCalL1SeededToken_(consumes<std::vector<reco::SuperCluster>>(iConfig.getUntrackedParameter<edm::InputTag>("scHGCalL1Seeded")))

{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

EGammaNtuples::~EGammaNtuples() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//



// ------------ method called for each event  ------------
void EGammaNtuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Auxiliary Data
  std::cout << "----- Auxiliary Data ----- " << std::endl;
  std::cout << "Run Nr: " << iEvent.eventAuxiliary().id().run() << std::endl;
  std::cout << "LumiSec: " << iEvent.eventAuxiliary().id().luminosityBlock() << std::endl;
  std::cout << "Event Nr: " << iEvent.eventAuxiliary().id().event() << std::endl;

  // Gen Particle
  std::cout << "----- Gen Particles ----- " << std::endl;
  for (const auto& gp: iEvent.get(genParticleToken_)){
    std::cout << "Energy: " << gp.energy() << std::endl;
    std::cout << "PT: " << gp.pt() << std::endl;
    std::cout << "Eta: " << gp.eta() << std::endl;
    std::cout << "Phi: " << gp.phi() << std::endl;
    std::cout << "Vz: " << gp.vz() << std::endl;
  }

  // Eg Object
  std::cout << "----- Egamma Particles ----- " << std::endl;
  for (const auto& egobj: iEvent.get(eGammaObjectToken_)){
    std::cout << "Energy: " << egobj.energy() << std::endl;
    std::cout << "ET: "<< egobj.et() << std::endl;
    std::cout << "Eta: "<< egobj.eta() << std::endl;
    std::cout << "Phi: "<< egobj.phi() << std::endl;
  }

  // SuperCluster

  std::cout << "----- SuperClusters: Barrel, L1Seeded ----- " << std::endl;
  for (const auto& sc: iEvent.get(scBarrelL1SeededToken_)){
    std::cout << "rawEnergy: " << sc.rawEnergy() << std::endl;
    std::cout << "nrClus: " << sc.clusters().size() << std::endl;
    std::cout << "seedId: " << sc.seed()->seed().rawId() << std::endl;
    std::cout << "seedDet: " << sc.seed()->seed().det() << std::endl;
    std::cout << "clusterMaxDr: " << sc.clusterMaxDr() << std::endl;
    std::cout << "r9Frac: " << sc.r9Frac() << std::endl;
    std::cout << "isEb: " << sc.isEB() << std::endl;
    std::cout << "isEE: " << sc.isEE() << std::endl;
  }

  std::cout << "----- SuperClusters: HGCAL, Unseeded  ----- " << std::endl;
  for (const auto& sc: iEvent.get(scHGCalL1SeededToken_)){
    std::cout << "rawEnergy: " << sc.rawEnergy() << std::endl;
    std::cout << "nrClus: " << sc.clusters().size() << std::endl;
    std::cout << "seedId: " << sc.seed()->seed().rawId() << std::endl;
    std::cout << "seedDet: " << sc.seed()->seed().det() << std::endl;
    std::cout << "clusterMaxDr: " << sc.clusterMaxDr() << std::endl;
    std::cout << "r9Frac: " << sc.r9Frac() << std::endl;
    std::cout << "isEb: " << sc.isEB() << std::endl;
    std::cout << "isEE: " << sc.isEE() << std::endl;
  }



#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void EGammaNtuples::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void EGammaNtuples::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EGammaNtuples::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EGammaNtuples);
