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
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/EgammaObject.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class EGammaNtuples : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit EGammaNtuples(const edm::ParameterSet&);
  ~EGammaNtuples() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  bool isEE(const reco::SuperCluster&);
  bool isEB(const reco::SuperCluster&);
  float cal_cluster_maxdr(const reco::SuperCluster&);
  float cal_r9(const reco::SuperCluster&, const EcalRecHitCollection&, const CaloTopology&);

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleToken_;
  edm::EDGetTokenT<std::vector<trigger::EgammaObject>> eGammaObjectToken_;
  edm::EDGetTokenT<std::vector<trigger::EgammaObject>> eGammaObjectUnseededToken_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>> scBarrelL1SeededToken_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>> scHGCalL1SeededToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ebRecHitsToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> eeRecHitsToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopoToken_;

  TTree *tree = new TTree("tree","tree");

  std::vector<int> run_nr;
  std::vector<int> lumi_sec;
  std::vector<int> event_nr;

  std::vector<float> gen_energy;
  std::vector<float> gen_pt;
  std::vector<float> gen_eta;
  std::vector<float> gen_phi;
  std::vector<float> gen_vz;

  std::vector<float> eg_energy;
  std::vector<float> eg_et;
  std::vector<float> eg_eta;
  std::vector<float> eg_phi;

  std::vector<float> sc_rawEnergy;
  std::vector<int> sc_nrClus;
  std::vector<int> sc_seedId;
  std::vector<int> sc_seedDet;
  std::vector<float> sc_clusterMaxDr;
  std::vector<float> sc_r9;
  std::vector<int> sc_isEB;
  std::vector<int> sc_isEE; 

  const CaloTopology* ecalTopology_;
  const CaloGeometry* caloGeometry_;


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
  scHGCalL1SeededToken_(consumes<std::vector<reco::SuperCluster>>(iConfig.getUntrackedParameter<edm::InputTag>("scHGCalL1Seeded"))),
  ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getUntrackedParameter<edm::InputTag>("ebRecHits"))),
  eeRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getUntrackedParameter<edm::InputTag>("eeRecHits"))),
  caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
  caloTopoToken_(esConsumes<CaloTopology, CaloTopologyRecord>()) {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed

  usesResource("TFileService");
  edm::Service<TFileService> file;
  tree = file->make<TTree>("tree","tree");

  tree->Branch("run_nr", &run_nr);
  tree->Branch("lumi_sec", &lumi_sec);
  tree->Branch("event_nr", &event_nr);

  tree->Branch("gen_energy", &gen_energy);
  tree->Branch("gen_pt", &gen_pt);
  tree->Branch("gen_eta", &gen_eta);
  tree->Branch("gen_phi", &gen_phi);
  tree->Branch("gen_vz", &gen_vz);

  tree->Branch("eg_energy", &eg_energy);
  tree->Branch("eg_et", &eg_et);
  tree->Branch("eg_eta", &eg_eta);
  tree->Branch("eg_phi", &eg_phi);

  tree->Branch("sc_rawEnergy", &sc_rawEnergy);
  tree->Branch("sc_nrClus", &sc_nrClus);
  tree->Branch("sc_seedId", &sc_seedId);
  tree->Branch("sc_seedDet", &sc_seedDet);
  tree->Branch("sc_clusterMaxDr", &sc_clusterMaxDr);
  tree->Branch("sc_r9", &sc_r9);
  tree->Branch("sc_isEB", &sc_isEB);
  tree->Branch("sc_isEE", &sc_isEE);
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


float EGammaNtuples::cal_r9(const reco::SuperCluster& sc,
                                    const EcalRecHitCollection& recHits,
                                    const CaloTopology& topology){
  auto &cluster = sc.seed();
  float e3x3 = EcalClusterTools::e3x3(*cluster, &recHits, &topology);
  return e3x3/sc.rawEnergy();
}

/*
float EGammaNtuples::cal_r9(const reco::SuperCluster& sc, const edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>& hits){
  
  float seed_eta = sc.seed()->eta();
  float seed_id = sc.seed()->seed().rawId();
  if (seed_id.det()!=DetId::Ecal || sc.rawEnergy()==0){
    return 0;
  }
  float e3x3 = 0;

  seed_id = EBDetId(seed_id);
  
  for local_ieta
  }
}

float EGammaNtuples::cal_r28(const reco::SuperCluster& sc, bool)

*/


float EGammaNtuples::cal_cluster_maxdr(const reco::SuperCluster& sc){
  float max_dr2 = 0;
  float seed_eta = sc.seed()->eta();
  float seed_phi = sc.seed()->phi();
  for(auto& c : sc.clusters()){
    if (c == sc.seed()){
      continue;
    }
    float dr2 = reco::deltaR2(c->eta(),c->phi(),seed_eta,seed_phi);
    max_dr2 = std::max(max_dr2,dr2);
  }

  // ECAL takes 999. if no other cluster for maxDR2
  if (max_dr2==0 && (sc.seed()->seed().det()==DetId::Ecal || sc.seed()->seed().det()==DetId::HGCalEE)){
    return 999;
  } else {
    return std::sqrt(max_dr2);
  }
}

bool EGammaNtuples::isEB(const reco::SuperCluster& sc){
  bool eb = false;
  if (sc.seed()->hitsAndFractions()[0].first.subdetId() == 1) {
    eb = true;
  } 
  return eb;
}

bool EGammaNtuples::isEE(const reco::SuperCluster& sc){
  bool ee = false;
  if (sc.seed()->hitsAndFractions()[0].first.subdetId() == 0) {
    ee = true;
  } 
  return ee;
}

// ------------ method called for each event  ------------
void EGammaNtuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<EcalRecHitCollection> ebRecHitsHandle;  
  iEvent.getByToken(ebRecHitsToken_, ebRecHitsHandle);

  edm::Handle<EcalRecHitCollection> eeRecHitsHandle;
  iEvent.getByToken(eeRecHitsToken_, eeRecHitsHandle);

  const CaloGeometry &geom = iSetup.getData(caloGeomToken_);
  const CaloTopology &topo = iSetup.getData(caloTopoToken_);

  // Auxiliary Data
  std::cout << "----- Auxiliary Data ----- " << std::endl;
  std::cout << "Run Nr: " << iEvent.eventAuxiliary().id().run() << std::endl;
  std::cout << "LumiSec: " << iEvent.eventAuxiliary().id().luminosityBlock() << std::endl;
  std::cout << "Event Nr: " << iEvent.eventAuxiliary().id().event() << std::endl;

  run_nr.push_back(iEvent.eventAuxiliary().id().run());
  run_nr.push_back(iEvent.eventAuxiliary().id().luminosityBlock());
  run_nr.push_back(iEvent.eventAuxiliary().id().event());

  // Gen Particle
  std::cout << "----- Gen Particles ----- " << std::endl;
  for (const auto& gp: iEvent.get(genParticleToken_)){
    std::cout << "Energy: " << gp.energy() << std::endl;
    std::cout << "PT: " << gp.pt() << std::endl;
    std::cout << "Eta: " << gp.eta() << std::endl;
    std::cout << "Phi: " << gp.phi() << std::endl;
    std::cout << "Vz: " << gp.vz() << std::endl;
  
    gen_energy.push_back(gp.energy());
    gen_pt.push_back(gp.pt());
    gen_eta.push_back(gp.eta());
    gen_phi.push_back(gp.phi());
    gen_vz.push_back(gp.vz());
  }

  // Eg Object
  std::cout << "----- Egamma Particles ----- " << std::endl;
  for (const auto& egobj: iEvent.get(eGammaObjectToken_)){
    std::cout << "Energy: " << egobj.energy() << std::endl;
    std::cout << "ET: "<< egobj.et() << std::endl;
    std::cout << "Eta: "<< egobj.eta() << std::endl;
    std::cout << "Phi: "<< egobj.phi() << std::endl;

    eg_energy.push_back(egobj.energy());
    eg_et.push_back(egobj.et());
    eg_eta.push_back(egobj.eta());
    eg_phi.push_back(egobj.phi());
  }

  // SuperCluster

  std::cout << "----- SuperClusters: Barrel, L1Seeded ----- " << std::endl;
  for (const auto& sc: iEvent.get(scBarrelL1SeededToken_)){
    std::cout << "rawEnergy: " << sc.rawEnergy() << std::endl;
    std::cout << "nrClus: " << sc.clusters().size() << std::endl;
    std::cout << "seedId: " << sc.seed()->seed().rawId() << std::endl;
    std::cout << "seedDet: " << sc.seed()->seed().det() << std::endl;
    std::cout << "clusterMaxDr: " << cal_cluster_maxdr(sc) << std::endl;
    std::cout << "r9Frac: " << cal_r9(sc, *ebRecHitsHandle, topo) << std::endl;
    std::cout << "isEB: " << isEB(sc) << std::endl;
    std::cout << "isEE: " << isEE(sc) << std::endl;

    sc_rawEnergy.push_back(sc.rawEnergy());
    sc_nrClus.push_back(sc.clusters().size());
    sc_seedId.push_back(sc.seed()->seed().rawId());
    sc_seedDet.push_back(sc.seed()->seed().det());
    sc_clusterMaxDr.push_back(cal_cluster_maxdr(sc));
    sc_r9.push_back(cal_r9(sc, *ebRecHitsHandle, topo));
    sc_isEB.push_back(isEB(sc));
    sc_isEE.push_back(isEE(sc));
  }

  std::cout << "----- SuperClusters: HGCAL, Unseeded  ----- " << std::endl;
  for (const auto& sc: iEvent.get(scHGCalL1SeededToken_)){
    std::cout << "rawEnergy: " << sc.rawEnergy() << std::endl;
    std::cout << "nrClus: " << sc.clusters().size() << std::endl;
    std::cout << "seedId: " << sc.seed()->seed().rawId() << std::endl;
    std::cout << "seedDet: " << sc.seed()->seed().det() << std::endl;
    std::cout << "clusterMaxDr: " << cal_cluster_maxdr(sc) << std::endl;
    std::cout << "r9Frac: " << cal_r9(sc, *eeRecHitsHandle, topo) << std::endl;
    std::cout << "isEb: " << isEB(sc) << std::endl;
    std::cout << "isEE: " << isEE(sc) << std::endl;

    sc_rawEnergy.push_back(sc.rawEnergy());
    sc_nrClus.push_back(sc.clusters().size());
    sc_seedId.push_back(sc.seed()->seed().rawId());
    sc_seedDet.push_back(sc.seed()->seed().det());
    sc_clusterMaxDr.push_back(cal_cluster_maxdr(sc));
    sc_r9.push_back(cal_r9(sc, *ebRecHitsHandle, topo));
    sc_isEB.push_back(isEB(sc));
    sc_isEE.push_back(isEE(sc));
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
  tree->Write();
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
