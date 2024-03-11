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
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/Common/interface/getRef.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

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
  float cal_r9(const reco::SuperCluster&, const EcalRecHitCollection&, const CaloTopology&, const bool);
  float cal_r28(const reco::SuperCluster&, const bool);
  virtual void fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap, 
      const HGCRecHitCollection& rechitsHGCEE) const;
  std::vector<DetId> getNeighbors(const DetId, const float, const CaloSubdetectorGeometry*);

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleToken_;
  edm::EDGetTokenT<std::vector<trigger::EgammaObject>> eGammaObjectToken_;
  edm::EDGetTokenT<std::vector<trigger::EgammaObject>> eGammaObjectUnseededToken_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>> scBarrelL1SeededToken_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>> scHGCalL1SeededToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ebRecHitsToken_;
  edm::EDGetTokenT<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit>>> eeRecHitsToken_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> sigmaIEtaIEtaToken_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> sigmaIPhiIPhiToken_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> sigmaIEtaIEtaNoiseCleanedToken_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> sigmaIPhiIPhiNoiseCleanedToken_;
  edm::EDGetTokenT<int> nrHitsEB1GeVToken_;
  edm::EDGetTokenT<int> nrHitsEE1GeVToken_;
  edm::EDGetTokenT<std::vector<reco::RecoEcalCandidate>> eGammaCandidatesToken_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopoToken_;


  const CaloTopology* ecalTopology_;
  const CaloGeometry* caloGeometry_;
  std::map<DetId, const HGCRecHit*> hitMap;

  hgcal::RecHitTools recHitTools_;

  TFile *newfile = new TFile("output.root", "RECREATE","",207);
  TTree *tree = new TTree("egHLTRun3Tree","egHLTRun3Tree");

  int run_nr;
  int lumi_sec;
  int event_nr;

  int nrHitsEB1GeV;
  int nrHitsEE1GeV;

  float eg_sigmaIEtaIEta;
  float eg_sigmaIPhiIPhi;
  float eg_sigmaIEtaIEtaNoiseCleaned;
  float eg_sigmaIPhiIPhiNoiseCleaned;

  float gen_energy;
  float gen_pt;
  float gen_eta;
  float gen_phi;
  float gen_vz;

  float eg_energy;
  float eg_et;
  float eg_eta;
  float eg_phi;

  float sc_rawEnergy;
  int sc_nrClus;
  int sc_seedId;
  int sc_seedDet;
  float sc_clusterMaxDr;
  float sc_r9Frac;
  float sc_r9Full;
  int sc_isEB;
  int sc_isEE; 
  float sc_phiWidth; 

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
  eeRecHitsToken_(consumes<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit>>>(iConfig.getUntrackedParameter<edm::InputTag>("eeRecHits"))),
  sigmaIEtaIEtaToken_(consumes<reco::RecoEcalCandidateIsolationMap> (iConfig.getUntrackedParameter<edm::InputTag>("sigmaIEtaIEta"))),
  sigmaIPhiIPhiToken_(consumes<reco::RecoEcalCandidateIsolationMap> (iConfig.getUntrackedParameter<edm::InputTag>("sigmaIPhiIPhi"))),
  sigmaIEtaIEtaNoiseCleanedToken_(consumes<reco::RecoEcalCandidateIsolationMap> (iConfig.getUntrackedParameter<edm::InputTag>("sigmaIEtaIEtaNoiseCleaned"))),
  sigmaIPhiIPhiNoiseCleanedToken_(consumes<reco::RecoEcalCandidateIsolationMap> (iConfig.getUntrackedParameter<edm::InputTag>("sigmaIPhiIPhiNoiseCleaned"))),
  nrHitsEB1GeVToken_(consumes<int>(iConfig.getUntrackedParameter<edm::InputTag>("nrHitsEB1GeV"))),
  nrHitsEE1GeVToken_(consumes<int>(iConfig.getUntrackedParameter<edm::InputTag>("nrHitsEE1GeV"))),
  eGammaCandidatesToken_(consumes<std::vector<reco::RecoEcalCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("eGammaCandidates"))),
  caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
  caloTopoToken_(esConsumes<CaloTopology, CaloTopologyRecord>()) {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed


  TTree *newtree = new TTree("egHLTRun3Tree","egHLTRun3Tree");

  tree->Branch("runnr", &run_nr);
  tree->Branch("lumiSec", &lumi_sec);
  tree->Branch("eventnr", &event_nr);

  tree->Branch("nrHitsEB1GeV", &nrHitsEB1GeV);
  tree->Branch("nrHitsEE1GeV", &nrHitsEE1GeV);  
  
  tree->Branch("eg_et", &eg_et);
  tree->Branch("eg_energy", &eg_energy);
  tree->Branch("eg_rawEnergy", &sc_rawEnergy);
  tree->Branch("eg_eta", &eg_eta);
  tree->Branch("eg_phi", &eg_phi);
  tree->Branch("eg_phiWidth", &sc_phiWidth);
  tree->Branch("eg_nrClus", &sc_nrClus);
  tree->Branch("eg_seedId", &sc_seedId);
  tree->Branch("eg_seedDet", &sc_seedDet);
  tree->Branch("eg_sigmaIEtaIEta", &eg_sigmaIEtaIEta);
  tree->Branch("eg_sigmaIPhiIPhi", &eg_sigmaIPhiIPhi);
  tree->Branch("eg_sigmaIEtaIEtaNoise", &eg_sigmaIEtaIEtaNoiseCleaned);
  tree->Branch("eg_sigmaIPhiIPhiNoise", &eg_sigmaIPhiIPhiNoiseCleaned);
  tree->Branch("eg_clusterMaxDR", &sc_clusterMaxDr);
  tree->Branch("eg_r9Frac", &sc_r9Frac);
  tree->Branch("eg_r9Full", &sc_r9Full);
  tree->Branch("eg_isEB", &sc_isEB);
  tree->Branch("eg_isEE", &sc_isEE);

  tree->Branch("eg_gen_energy", &gen_energy);
  tree->Branch("eg_gen_pt", &gen_pt);
  tree->Branch("eg_gen_eta", &gen_eta);
  tree->Branch("eg_gen_phi", &gen_phi);
  tree->Branch("eg_gen_vz", &gen_vz);
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

/*
float EGammaNtuples::cal_r28(const reco::SuperCluster& sc,
                                    const CaloSubdetectorGeometry* subGeom)

  // Find seeds and collect RecHits in Cluster

  // Loop over clusters

    // get Neighbors from seed
  
    // add to e_28 and e

  // return frac
*/

float EGammaNtuples::cal_r9(const reco::SuperCluster& sc,
                                    const EcalRecHitCollection& recHits,
                                    const CaloTopology& topology,
                                    const bool frac){
  auto &cluster = sc.seed();
  float e3x3 = EcalClusterTools::e3x3(*cluster, &recHits, &topology);
  if (frac == false){
    return e3x3;
  }
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
*/

float EGammaNtuples::cal_r28(const reco::SuperCluster& sc,
                              const bool frac){
    float e28 = 0;
    float e = 0;
    float dR = 2.8;
    
    auto clusters = sc.clusters();
    for (auto& c : clusters){
      std::map<int, std::vector<DetId>> v_detid;
      std::vector<int> seedIds(47,0);
      std::vector<float> maxEnergy(47,0);

      // Find seeds of LCs
      for (auto& hf : c->hitsAndFractions()){
        auto h = hf.first;
        auto f = hf.second;
        int layer = recHitTools_.getLayerWithOffset(h) ;
        std::map<DetId,const HGCRecHit *>::const_iterator itcheck = hitMap.find(h); 
        if (itcheck != hitMap.end()){
          v_detid[layer].push_back(h);
          float energy = hitMap[h]->energy();
          e+=energy;
          if (maxEnergy[layer] < energy){
            seedIds[layer] = h.rawId();
            maxEnergy[layer] = energy;
          }
        }
      }

      // Collect cells in a given radius
      for (auto& seedId: seedIds){
        if (seedId==0) continue;
        int layer = recHitTools_.getLayerWithOffset(seedId);
        GlobalPoint p0 = recHitTools_.getPosition(seedId);
        for (auto& h: v_detid[layer]){
          GlobalPoint p =  recHitTools_.getPosition(h);
          float dist = std::sqrt((p0.x()-p.x())*(p0.x()-p.x()) + (p0.y()-p.y())*(p0.y()-p.y()) + (p0.z()-p.z())*(p0.z()-p.z()));
          if (dist<dR){
            e28 += hitMap[h]->energy();
          }
        }
      }
    }
  if (frac == false){
    return e28;
  }
  return e28/e;
}


std::vector<DetId> EGammaNtuples::getNeighbors(const DetId seedId, const float dR, const CaloSubdetectorGeometry* subGeom){

  const GlobalPoint pos = recHitTools_.getPosition(seedId);
  const int layer = recHitTools_.getLayerWithOffset(seedId);
  const int zside = recHitTools_.zside(seedId);
  const double dR2 = dR * dR;
  const double eta = pos.eta();
  const double phi = pos.phi();

  std::vector<DetId> neighbors;
  if (0.000001 < dR) {
    for (auto id : subGeom->getValidDetIds()) {
      const GlobalPoint& p = recHitTools_.getPosition(id);
      const int l = recHitTools_.getLayerWithOffset(id);
      if (l != layer && recHitTools_.zside(id)==zside) continue;
      const float eta0 = p.eta();
      if (fabs(eta - eta0) < dR) {
        const float phi0 = p.phi();
        float delp =fabs(phi - phi0);
        if (delp > M_PI)
          delp = 2 * M_PI - delp;
        if (delp < dR) {
          const float dist2 = reco::deltaR2(eta0, phi0, eta, phi);
          if (dist2 < dR2){
            neighbors.push_back(id);
            //dss.insert(m_validIds[i]);
          }
        }
      }
    }
  }
  return neighbors;
}

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

void EGammaNtuples::fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
                                const HGCRecHitCollection& rechitsHGCEE) const {
  hitMap.clear();
  for (const auto& hit : rechitsHGCEE) {
    hitMap.emplace(hit.detid(), &hit);
  }
} // end of EfficiencyStudies::fillHitMap

// ------------ method called for each event  ------------
void EGammaNtuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<EcalRecHitCollection> ebRecHitsHandle;  
  iEvent.getByToken(ebRecHitsToken_, ebRecHitsHandle);

  edm::Handle<HGCRecHitCollection> eeRecHitsHandle;
  iEvent.getByToken(eeRecHitsToken_, eeRecHitsHandle);

  const CaloGeometry &geom = iSetup.getData(caloGeomToken_);
  const CaloTopology &topo = iSetup.getData(caloTopoToken_);

  recHitTools_.setGeometry(geom);
  int hgcalEEId = DetId::HGCalEE;
  const CaloSubdetectorGeometry *subGeom = geom.getSubdetectorGeometry(DetId::Detector(hgcalEEId), ForwardSubdetector::ForwardEmpty);

  edm::Handle<reco::RecoEcalCandidateIsolationMap> sigmaIEtaIEtaHandle;
  iEvent.getByToken(sigmaIEtaIEtaToken_,sigmaIEtaIEtaHandle);

  edm::Handle<reco::RecoEcalCandidateIsolationMap> sigmaIPhiIPhiHandle;
  iEvent.getByToken(sigmaIPhiIPhiToken_,sigmaIPhiIPhiHandle);

  edm::Handle<reco::RecoEcalCandidateIsolationMap> sigmaIEtaIEtaNoiseCleanedHandle;
  iEvent.getByToken(sigmaIEtaIEtaNoiseCleanedToken_,sigmaIEtaIEtaNoiseCleanedHandle);

  edm::Handle<reco::RecoEcalCandidateIsolationMap> sigmaIPhiIPhiNoiseCleanedHandle;
  iEvent.getByToken(sigmaIPhiIPhiNoiseCleanedToken_,sigmaIPhiIPhiNoiseCleanedHandle);

  edm::Handle<std::vector<reco::RecoEcalCandidate>> eGammaCandidatesHandle;
  iEvent.getByToken(eGammaCandidatesToken_,eGammaCandidatesHandle);

  std::cout << "EgammaCandidates size:" << (*eGammaCandidatesHandle).size() << std::endl;
  std::cout << "ValueMap size:" << (*sigmaIEtaIEtaHandle).size() << std::endl;

  fillHitMap(hitMap, *eeRecHitsHandle);

  // Objects
  const auto& egObjects = iEvent.get(eGammaObjectToken_);
  const auto& gPs =  iEvent.get(genParticleToken_);
  const auto& bScs = iEvent.get(scBarrelL1SeededToken_);
  const auto& eScs = iEvent.get(scHGCalL1SeededToken_);

  std::cout << "Size bScs: " << bScs.size() << std::endl;
  std::cout << "Size eScs: " << eScs.size() << std::endl;

  std::vector<reco::SuperCluster> scs(bScs.size() + eScs.size()); 
  std::merge(bScs.begin(), bScs.end(), eScs.begin(), eScs.end(), scs.begin());

  std::cout << "egObjs: " << egObjects.size() << ", gPs: " << gPs.size() << ", bScs: " << bScs.size() << ", eScs: " << eScs.size() << std::endl;

  if (egObjects.size() == 2 && gPs.size() == 2 && scs.size()){
    for (int i = 0; i<static_cast<int>((*eGammaCandidatesHandle).size());i++){
      // Auxiliary Data
      std::cout << "----- Auxiliary Data ----- " << std::endl;
      std::cout << "Run Nr: " << iEvent.eventAuxiliary().id().run() << std::endl;
      std::cout << "LumiSec: " << iEvent.eventAuxiliary().id().luminosityBlock() << std::endl;
      std::cout << "Event Nr: " << iEvent.eventAuxiliary().id().event() << std::endl;

      run_nr = iEvent.eventAuxiliary().id().run();
      lumi_sec = iEvent.eventAuxiliary().id().luminosityBlock();
      event_nr = iEvent.eventAuxiliary().id().event();

      // NrHits
      std::cout << "------ Nr Hits -------------" << std::endl;
      std::cout << "nrHitsEB1GeV: " << iEvent.get(nrHitsEB1GeVToken_) << std::endl;
      std::cout << "nrHitsEE1GeV: " << iEvent.get(nrHitsEE1GeVToken_) << std::endl;

      // Gen Particle
      std::cout << "----- Gen Particles ----- " << std::endl;
      std::cout << "Energy: " << gPs[i].energy() << std::endl;
      std::cout << "PT: " << gPs[i].pt() << std::endl;
      std::cout << "Eta: " << gPs[i].eta() << std::endl;
      std::cout << "Phi: " << gPs[i].phi() << std::endl;
      std::cout << "Vz: " << gPs[i].vz() << std::endl;

      gen_energy = gPs[i].energy();
      gen_pt = gPs[i].pt();
      gen_eta = gPs[i].eta();
      gen_phi = gPs[i].phi();
      gen_vz = gPs[i].vz();
      
      // Eg Object
      std::cout << "----- Egamma Particles ----- " << std::endl;
      std::cout << "Energy: " << egObjects[i].energy() << std::endl;
      std::cout << "ET: "<< egObjects[i].et() << std::endl;
      std::cout << "Eta: "<< egObjects[i].eta() << std::endl;
      std::cout << "Phi: "<< egObjects[i].phi() << std::endl;

      eg_energy = egObjects[i].energy();
      eg_et = egObjects[i].et();
      eg_eta = egObjects[i].eta();
      eg_phi = egObjects[i].phi();

      // Sigma
      std::cout << "----- Sigma Variables -----" << std::endl;
      reco::RecoEcalCandidateRef candidateRef = getRef(eGammaCandidatesHandle, i);
      std::cout << (*sigmaIEtaIEtaHandle)[candidateRef] << std::endl;
      std::cout << (*sigmaIPhiIPhiHandle)[candidateRef] << std::endl;
      std::cout << (*sigmaIEtaIEtaNoiseCleanedHandle)[candidateRef] << std::endl;
      std::cout << (*sigmaIPhiIPhiNoiseCleanedHandle)[candidateRef] << std::endl;

      eg_sigmaIEtaIEta = (*sigmaIEtaIEtaHandle)[candidateRef];
      eg_sigmaIPhiIPhi = (*sigmaIPhiIPhiHandle)[candidateRef];
      eg_sigmaIEtaIEtaNoiseCleaned = (*sigmaIEtaIEtaNoiseCleanedHandle)[candidateRef];
      eg_sigmaIPhiIPhi = (*sigmaIPhiIPhiNoiseCleanedHandle)[candidateRef];
      //auto mapIt = valueMapHandles->find(candidateRef);

      // SuperCluster
      std::cout << "----- SuperClusters: Barrel, L1Seeded ----- " << std::endl;
      std::cout << "rawEnergy: " << scs[i].rawEnergy() << std::endl;
      std::cout << "nrClus: " << scs[i].clusters().size() << std::endl;
      std::cout << "seedId: " << scs[i].seed()->seed().rawId() << std::endl;
      std::cout << "seedDet: " << scs[i].seed()->seed().det() << std::endl;
      std::cout << "clusterMaxDr: " << cal_cluster_maxdr(scs[i]) << std::endl;
      std::cout << "r9Frac: " << cal_r9(scs[i], *ebRecHitsHandle, topo,true) << std::endl;
      std::cout << "r9Full: " << cal_r9(scs[i], *ebRecHitsHandle, topo,false) << std::endl;
      std::cout << "isEB: " << isEB(scs[i]) << std::endl;
      std::cout << "isEE: " << isEE(scs[i]) << std::endl;
      std::cout << "PhiWidth: " << scs[i].phiWidth() << std::endl;

      sc_rawEnergy = scs[i].rawEnergy();
      sc_nrClus = scs[i].clusters().size();
      sc_seedId = scs[i].seed()->seed().rawId();
      sc_seedDet = scs[i].seed()->seed().det();
      sc_clusterMaxDr = cal_cluster_maxdr(scs[i]);
      sc_r9Frac = cal_r9(scs[i], *ebRecHitsHandle, topo,true);
      sc_r9Full = cal_r9(scs[i], *ebRecHitsHandle, topo,false);
      sc_isEB = isEB(scs[i]);
      sc_isEE = isEE(scs[i]);
      sc_phiWidth = scs[i].phiWidth();
      std::cout << "Nr. of Rechits: " << (*ebRecHitsHandle).size() << std::endl;
      tree->Fill();
    }
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
