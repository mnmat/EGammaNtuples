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
#include "FWCore/Utilities/interface/isFinite.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"


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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "TTree.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "RecoEgamma/EgammaTools/interface/EgammaHGCALIDParamDefaults.h"
#include "RecoEgamma/EgammaTools/interface/HGCalShowerShapeHelper.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace reco;

namespace {
  //bool is if a valid dr was found, float is the dr
  std::pair<bool, float> getMaxDRNonSeedCluster(const reco::SuperCluster& sc) {
    float maxDR2 = 0.;
    const edm::Ptr<reco::CaloCluster>& seedClus = sc.seed();

    for (const auto& clus : sc.clusters()) {
      if (clus == seedClus) {
        continue;
      }

      // find cluster with max dR
      const double dr2 = reco::deltaR2(*clus, *seedClus);
      if (dr2 > maxDR2) {
        maxDR2 = dr2;
      }
    }
    return {sc.clustersSize() != 1, sc.clustersSize() != 1 ? std::sqrt(maxDR2) : 999.};
  }
  template <typename T>
  int countRecHits(const T& recHitHandle, float threshold) {
    int count = 0;
    if (recHitHandle.isValid()) {
      for (const auto& recHit : *recHitHandle) {
        if (recHit.energy() > threshold) {
          count++;
        }
      }
    }
    return count;
  }
}  // namespace

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
  void fillRegDataEGRun3Tree(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillRegDataEcalV1(const reco::SuperCluster& sc, const EcalRecHitCollection& recHits,const reco::GenParticle& gPs);
  void fillRegDataEcalHLTV1(const reco::SuperCluster& sc, const EcalRecHitCollection& recHits,const reco::GenParticle& gPs);
  void fillRegDataHGCALV1(const reco::SuperCluster& sc,const reco::GenParticle& gPs);
  void fillRegDataHGCALHLTV1(const reco::SuperCluster& sc,const reco::GenParticle& gPs);

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
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  edm::EDGetTokenT<reco::PFRecHitCollection> pfRecHitsHGCALToken_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopoToken_;


  const CaloTopology* caloTopology_;
  const CaloGeometry* caloGeometry_;
  std::map<DetId, const HGCRecHit*> hitMap;

  hgcal::RecHitTools recHitTools_;

  HGCalShowerShapeHelper hgcalShowerShapes_;
  float hgcalCylinderR_ = EgammaHGCALIDParamDefaults::kRCylinder;

  edm::Handle<reco::VertexCollection> vertices_;

  
  TFile *newfile = new TFile("output.root", "RECREATE","",207);
  TTree* egRun3Tree = new TTree("egHLTRun3Tree","egHLTRun3Tree");
  TTree* egRun4CompleteTree = new TTree("egOfflineRun4Tree","egOfflineRun4Tree")
  TTree* egRegDataEcalV1Tree = new TTree("egRegDataEcalV1","egRegDataEcalV1");
  TTree* egRegDataEcalHLTV1Tree = new TTree("egRegDataEcalHLTV1","egRegDataEcalHLTV1");
  TTree* egRegDataHGCALV1Tree = new TTree("egRegDataHGCALV1","egRegDataHGCALV1");
  TTree* egRegDataHGCALHLTV1Tree = new TTree("egRegDataHGCALHLTV1","egRegDataHGCALHLTV1");

  // egRun3Tree

  int run_nr;
  int lumi_sec;
  int event_nr;

  int nrHitsEB1GeV;
  int nrHitsEE1GeV;
  int hitsEnergyThreshold_ = 1;

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
  float sc_regressedEnergy;
  int sc_nrClus;
  int sc_seedId;
  int sc_seedDet;
  float sc_clusterMaxDR;
  float sc_r9Frac;
  float sc_r9Full;
  int sc_isEB;
  int sc_isEE; 
  float sc_phiWidth; 

  // egRegDataEcalV1 Features

  int ecalV1_run_nr;
  int ecalV1_lumi_sec;
  int ecalV1_event_nr;

  int ecalV1_nrHitsThreshold;

  float ecalV1_sigmaIEtaIEta;
  float ecalV1_sigmaIPhiIPhi;
  float ecalV1_sigmaIEtaIPhi;

  float ecalV1_rawEnergy;
  float ecalV1_regressedEnergy;
  int ecalV1_nrVert;
  float ecalV1_etaWidth;
  float ecalV1_phiWidth;
  float ecalV1_energySeedCluster;
  float ecalV1_energyFirstRecHit;
  float ecalV1_energySecondRecHit;
  float ecalV1_eLeftRightDiffSumRatio;
  float ecalV1_eTopBottomDiffSumRatio;
  int ecalV1_rvar;
  float ecalV1_numberOfSubClusters;
  float ecalV1_clusterMaxDR;
  float ecalV1_clusterMaxDRDPhi;
  float ecalV1_clusterMaxDRDEta;
  float ecalV1_clusterMaxDREnergyFraction;
  float ecalV1_subCluster1RawEnergyFraction;
  float ecalV1_subCluster2RawEnergyFraction;
  float ecalV1_subCluster3RawEnergyFraction;
  float ecalV1_subCluster1DPhiToSeed;
  float ecalV1_subCluster2DPhiToSeed;
  float ecalV1_subCluster3DPhiToSeed;
  float ecalV1_subCluster1DEtaToSeed;
  float ecalV1_subCluster2DEtaToSeed;
  float ecalV1_subCluster3DEtaToSeed;
  float ecalV1_seedClusteriEtaOrX;
  float ecalV1_seedClusteriPhiOrY;
  float ecalV1_seedClusterEta;
  float ecalV1_gen_energy;
  float ecalV1_gen_pt;
  float ecalV1_gen_eta;
  float ecalV1_gen_phi;
  float ecalV1_gen_vz;

  float ecalV1_eg_sigmaIEtaIEta;  
  float ecalV1_eg_sigmaIPhiIPhi;
  int ecalV1_sc_isEB;
  int ecalV1_sc_isEE; 

  // egRegDataEcalHLTV1 Features

  int ecalHLTV1_run_nr;
  int ecalHLTV1_lumi_sec;
  int ecalHLTV1_event_nr;

  float ecalHLTV1_rawEnergy;
  float ecalHLTV1_regressedEnergy;
  float ecalHLTV1_phiWidth; 
  float ecalHLTV1_eta;
  float ecalHLTV1_rvar;
  int ecalHLTV1_numberOfSubClusters;
  float ecalHLTV1_clusterMaxDR;
  int ecalHLTV1_nrHitsThreshold;

  float ecalHLTV1_sc_eta;
  float ecalHLTV1_gen_energy;
  float ecalHLTV1_gen_pt;
  float ecalHLTV1_gen_eta;
  float ecalHLTV1_gen_phi;
  float ecalHLTV1_gen_vz;

  float ecalHLTV1_eg_sigmaIEtaIEta;  
  float ecalHLTV1_eg_sigmaIPhiIPhi;
  int ecalHLTV1_sc_isEB;
  int ecalHLTV1_sc_isEE; 

  // egRegDataHGCALV1 Features

  int hgcalV1_run_nr;
  int hgcalV1_lumi_sec;
  int hgcalV1_event_nr; 

  float hgcalV1_rawEnergy;
  float hgcalV1_regressedEnergy;
  float hgcalV1_eta;
  float hgcalV1_etaWidth;
  float hgcalV1_phiWidth;
  int hgcalV1_numberOfSubClusters;
  int hgcalV1_nrRecHits;
  float hgcalV1_clusterMaxDR;
  float hgcalV1_dEtaSCToSeed;
  float hgcalV1_dPhiSCToSeed;
  float hgcalV1_firstRecHitEnergyFraction;
  float hgcalV1_secondRecHitEnergyFraction;
  float hgcalV1_sigma2uu;
  float hgcalV1_sigma2vv;
  float hgcalV1_sigma2ww;
  float hgcalV1_rvar;
  float hgcalV1_seedEnergyFraction;
  float hgcalV1_nrHitsThreshold;
  float hgcalV1_gen_energy;
  float hgcalV1_gen_pt;
  float hgcalV1_gen_eta;
  float hgcalV1_gen_phi;
  float hgcalV1_gen_vz;

  float hgcalV1_eg_sigmaIEtaIEta;  
  float hgcalV1_eg_sigmaIPhiIPhi;
  int hgcalV1_sc_isEB;
  int hgcalV1_sc_isEE; 

  // egRegDataHGCALHLTV1 Features

  int hgcalHLTV1_run_nr;
  int hgcalHLTV1_lumi_sec;
  int hgcalHLTV1_event_nr; 

  float hgcalHLTV1_rawEnergy;
  float hgcalHLTV1_regressedEnergy;
  float hgcalHLTV1_eta;
  float hgcalHLTV1_phiWidth;
  int hgcalHLTV1_numberOfSubClusters;
  float hgcalHLTV1_clusterMaxDR;
  float hgcalHLTV1_rvar;
  float hgcalHLTV1_nrHitsThreshold;

  float hgcalHLTV1_gen_energy;
  float hgcalHLTV1_gen_pt;
  float hgcalHLTV1_gen_eta;
  float hgcalHLTV1_gen_phi;
  float hgcalHLTV1_gen_vz;

  float hgcalHLTV1_eg_sigmaIEtaIEta;  
  float hgcalHLTV1_eg_sigmaIPhiIPhi;
  int hgcalHLTV1_sc_isEB;
  int hgcalHLTV1_sc_isEE; 
  

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
  : genParticleToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  eGammaObjectToken_(consumes<std::vector<trigger::EgammaObject>>(iConfig.getParameter<edm::InputTag>("eGammaObjects"))),
  eGammaObjectUnseededToken_(consumes<std::vector<trigger::EgammaObject>>(iConfig.getParameter<edm::InputTag>("eGammaObjectsUnSeeded"))),
  scBarrelL1SeededToken_(consumes<std::vector<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("scBarrelL1Seeded"))),
  scHGCalL1SeededToken_(consumes<std::vector<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("scHGCalL1Seeded"))),
  ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("ebRecHits"))),
  eeRecHitsToken_(consumes<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit>>>(iConfig.getParameter<edm::InputTag>("eeRecHits"))),
  sigmaIEtaIEtaToken_(consumes<reco::RecoEcalCandidateIsolationMap> (iConfig.getParameter<edm::InputTag>("sigmaIEtaIEta"))),
  sigmaIPhiIPhiToken_(consumes<reco::RecoEcalCandidateIsolationMap> (iConfig.getParameter<edm::InputTag>("sigmaIPhiIPhi"))),
  sigmaIEtaIEtaNoiseCleanedToken_(consumes<reco::RecoEcalCandidateIsolationMap> (iConfig.getParameter<edm::InputTag>("sigmaIEtaIEtaNoiseCleaned"))),
  sigmaIPhiIPhiNoiseCleanedToken_(consumes<reco::RecoEcalCandidateIsolationMap> (iConfig.getParameter<edm::InputTag>("sigmaIPhiIPhiNoiseCleaned"))),
  nrHitsEB1GeVToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("nrHitsEB1GeV"))),
  nrHitsEE1GeVToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("nrHitsEE1GeV"))),
  eGammaCandidatesToken_(consumes<std::vector<reco::RecoEcalCandidate>>(iConfig.getParameter<edm::InputTag>("eGammaCandidates"))),
  verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  pfRecHitsHGCALToken_(consumes<reco::PFRecHitCollection>(iConfig.getParameter<edm::InputTag>("pfHGCALRecHits"))),
  caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
  caloTopoToken_(esConsumes<CaloTopology, CaloTopologyRecord>()){
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed

  hgcalShowerShapes_.setTokens<edm::Transition::Event>(consumesCollector());
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
                                    const CaloSubdetectorGeometry* subDetectorGeometry)

  // Find seeds and collect RecHits in Cluster

  // Loop over clusters

    // get Neighbors from seed
  
    // add to e_28 and e

  // return frac
*/

void EGammaNtuples::fillRegDataEGRun3Tree(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  edm::Handle<EcalRecHitCollection> ebRecHitsHandle;  
  iEvent.getByToken(ebRecHitsToken_, ebRecHitsHandle);

  edm::Handle<HGCRecHitCollection> eeRecHitsHandle;
  iEvent.getByToken(eeRecHitsToken_, eeRecHitsHandle);

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
      std::cout << "clusterMaxDR: " << cal_cluster_maxdr(scs[i]) << std::endl;
      std::cout << "r9Frac: " << cal_r9(scs[i], *ebRecHitsHandle, *caloTopology_,true) << std::endl;
      std::cout << "r9Full: " << cal_r9(scs[i], *ebRecHitsHandle, *caloTopology_,false) << std::endl;
      std::cout << "isEB: " << isEB(scs[i]) << std::endl;
      std::cout << "isEE: " << isEE(scs[i]) << std::endl;
      std::cout << "PhiWidth: " << scs[i].phiWidth() << std::endl;

      sc_rawEnergy = scs[i].rawEnergy();
      sc_regressedEnergy = scs[i].energy();
      sc_nrClus = scs[i].clusters().size();
      sc_seedId = scs[i].seed()->seed().rawId();
      sc_seedDet = scs[i].seed()->seed().det();
      sc_clusterMaxDR = cal_cluster_maxdr(scs[i]);
      sc_r9Frac = cal_r9(scs[i], *ebRecHitsHandle, *caloTopology_,true);
      sc_r9Full = cal_r9(scs[i], *ebRecHitsHandle, *caloTopology_,false);
      sc_isEB = isEB(scs[i]);
      sc_isEE = isEE(scs[i]);
      sc_phiWidth = scs[i].phiWidth();
      std::cout << "Nr. of Rechits: " << (*ebRecHitsHandle).size() << std::endl;
      egRun3Tree->Fill();
    }
  }
}

void EGammaNtuples::fillRegDataEcalV1(const reco::SuperCluster& sc,  const EcalRecHitCollection& recHits, const reco::GenParticle& gPs) {

  const reco::CaloCluster& seedCluster = *(sc.seed());
  const bool iseb = seedCluster.hitsAndFractions()[0].first.subdetId() == EcalBarrel;
  //const EcalRecHitCollection* recHits = iseb ? recHitsEB_.product() : recHitsEE_.product();

  const double raw_energy = sc.rawEnergy();
  const double regressed_energy = sc.energy();
  const int numberOfClusters = sc.clusters().size();

  const auto& localCovariances = EcalClusterTools::localCovariances(seedCluster, &recHits, caloTopology_);

  const float eLeft = EcalClusterTools::eLeft(seedCluster, &recHits, caloTopology_);
  const float eRight = EcalClusterTools::eRight(seedCluster, &recHits, caloTopology_);
  const float eTop = EcalClusterTools::eTop(seedCluster, &recHits, caloTopology_);
  const float eBottom = EcalClusterTools::eBottom(seedCluster, &recHits, caloTopology_);

  float sigmaIetaIeta = sqrt(localCovariances[0]);
  float sigmaIetaIphi = std::numeric_limits<float>::max();
  float sigmaIphiIphi = std::numeric_limits<float>::max();

  if (!edm::isNotFinite(localCovariances[2]))
    sigmaIphiIphi = sqrt(localCovariances[2]);

  // extra shower shapes
  bool applySigmaIetaIphiBug_ = false;
  const float see_by_spp = sigmaIetaIeta * (applySigmaIetaIphiBug_ ? std::numeric_limits<float>::max() : sigmaIphiIphi);
  if (see_by_spp > 0) {
    sigmaIetaIphi = localCovariances[1] / see_by_spp;
  } else if (localCovariances[1] > 0) {
    sigmaIetaIphi = 1.f;
  } else {
    sigmaIetaIphi = -1.f;
  }

  // calculate sub-cluster variables
  std::vector<float> clusterRawEnergy;
  clusterRawEnergy.resize(std::max(3, numberOfClusters), 0);
  std::vector<float> clusterDEtaToSeed;
  clusterDEtaToSeed.resize(std::max(3, numberOfClusters), 0);
  std::vector<float> clusterDPhiToSeed;
  clusterDPhiToSeed.resize(std::max(3, numberOfClusters), 0);
  float clusterMaxDR = 999.;
  float clusterMaxDRDPhi = 999.;
  float clusterMaxDRDEta = 999.;
  float clusterMaxDRRawEnergy = 0.;

  size_t iclus = 0;
  float maxDR = 0;
  edm::Ptr<reco::CaloCluster> pclus;
  const edm::Ptr<reco::CaloCluster>& theseed = sc.seed();
  // loop over all clusters that aren't the seed
  auto clusend = sc.clustersEnd();
  for (auto clus = sc.clustersBegin(); clus != clusend; ++clus) {
    pclus = *clus;

    if (theseed == pclus)
      continue;
    clusterRawEnergy[iclus] = pclus->energy();
    clusterDPhiToSeed[iclus] = reco::deltaPhi(pclus->phi(), theseed->phi());
    clusterDEtaToSeed[iclus] = pclus->eta() - theseed->eta();

    // find cluster with max dR
    const auto the_dr = reco::deltaR(*pclus, *theseed);
    if (the_dr > maxDR) {
      maxDR = the_dr;
      clusterMaxDR = maxDR;
      clusterMaxDRDPhi = clusterDPhiToSeed[iclus];
      clusterMaxDRDEta = clusterDEtaToSeed[iclus];
      clusterMaxDRRawEnergy = clusterRawEnergy[iclus];
    }
    ++iclus;
  }

  ecalV1_nrVert = vertices_->size();
  ecalV1_rawEnergy = raw_energy;
  ecalV1_regressedEnergy = regressed_energy;
  ecalV1_etaWidth = sc.etaWidth();
  ecalV1_phiWidth = sc.phiWidth();
  ecalV1_rvar = EcalClusterTools::e3x3(seedCluster, &recHits, caloTopology_) / raw_energy;
  ecalV1_energySeedCluster = seedCluster.energy() / raw_energy;
  ecalV1_energyFirstRecHit = EcalClusterTools::eMax(seedCluster, &recHits) / raw_energy;
  ecalV1_energySecondRecHit = EcalClusterTools::e2nd(seedCluster, &recHits) / raw_energy;
  ecalV1_eLeftRightDiffSumRatio = (eLeft + eRight != 0.f ? (eLeft - eRight) / (eLeft + eRight) : 0.f);
  ecalV1_eTopBottomDiffSumRatio = (eTop + eBottom != 0.f ? (eTop - eBottom) / (eTop + eBottom) : 0.f);
  ecalV1_sigmaIEtaIEta = sigmaIetaIeta;
  ecalV1_sigmaIEtaIPhi = sigmaIetaIphi;
  ecalV1_sigmaIPhiIPhi = sigmaIphiIphi;
  ecalV1_numberOfSubClusters = std::max(0, numberOfClusters - 1);
  ecalV1_clusterMaxDR = clusterMaxDR;
  ecalV1_clusterMaxDRDPhi = clusterMaxDRDPhi;
  ecalV1_clusterMaxDRDEta = clusterMaxDRDEta;
  ecalV1_clusterMaxDREnergyFraction = clusterMaxDRRawEnergy / raw_energy;
  ecalV1_subCluster1RawEnergyFraction = clusterRawEnergy[0] / raw_energy;
  ecalV1_subCluster2RawEnergyFraction = clusterRawEnergy[1] / raw_energy;
  ecalV1_subCluster3RawEnergyFraction = clusterRawEnergy[2] / raw_energy;
  ecalV1_subCluster1DPhiToSeed = clusterDPhiToSeed[0];
  ecalV1_subCluster2DPhiToSeed = clusterDPhiToSeed[1];
  ecalV1_subCluster3DPhiToSeed = clusterDPhiToSeed[2];
  ecalV1_subCluster1DEtaToSeed = clusterDEtaToSeed[0];
  ecalV1_subCluster2DEtaToSeed = clusterDEtaToSeed[1];
  ecalV1_subCluster3DEtaToSeed = clusterDEtaToSeed[2];
  if (iseb) {
    EBDetId ebseedid(seedCluster.seed());
    ecalV1_seedClusteriEtaOrX = ebseedid.ieta();
    ecalV1_seedClusteriPhiOrY = ebseedid.iphi();
  } else {
    EEDetId eeseedid(seedCluster.seed());
    ecalV1_seedClusteriEtaOrX = eeseedid.ix();
    ecalV1_seedClusteriPhiOrY = eeseedid.iy();
    //seed cluster eta is only needed for the 106X Ultra Legacy regressions
    //and was not used in the 74X regression however as its just an extra varaible
    //at the end, its harmless to add for the 74X regression
    ecalV1_seedClusterEta = seedCluster.eta();
  }

  ecalV1_gen_energy = gPs.energy();
  ecalV1_gen_pt = gPs.pt();
  ecalV1_gen_eta = gPs.eta();
  ecalV1_gen_phi = gPs.phi();
  ecalV1_gen_vz = gPs.vz();


  ecalV1_eg_sigmaIEtaIEta = eg_sigmaIEtaIEta;  
  ecalV1_eg_sigmaIPhiIPhi = eg_sigmaIPhiIPhi;  
  ecalV1_isEB = true;
  ecalV1_isEE = false;

  egRegDataEcalV1Tree->Fill();
}

void EGammaNtuples::fillRegDataEcalHLTV1(const reco::SuperCluster& sc,const EcalRecHitCollection& recHits, const reco::GenParticle& gPs) {
  std::vector<float> eval(7, 0.);
  auto maxDRNonSeedClus = getMaxDRNonSeedCluster(sc);
  const float clusterMaxDR = maxDRNonSeedClus.first ? maxDRNonSeedClus.second : 999.;

  const reco::CaloCluster& seedCluster = *(sc.seed());
  const bool iseb = seedCluster.hitsAndFractions()[0].first.subdetId() == EcalBarrel;
  //const EcalRecHitCollection* recHits = iseb ? recHitsEB_.product() : recHitsEE_.product();

  const auto& localCovariances = EcalClusterTools::localCovariances(seedCluster, &recHits, caloTopology_);
  float sigmaIetaIeta = sqrt(localCovariances[0]);
  float sigmaIetaIphi = std::numeric_limits<float>::max();
  float sigmaIphiIphi = std::numeric_limits<float>::max();

  if (!edm::isNotFinite(localCovariances[2]))
    sigmaIphiIphi = sqrt(localCovariances[2]);


  ecalHLTV1_nrHitsThreshold = nrHitsEB1GeV + nrHitsEE1GeV;
  ecalHLTV1_eta = sc.eta();
  ecalHLTV1_phiWidth = sc.phiWidth();
  ecalHLTV1_rvar = EcalClusterTools::e3x3(seedCluster, &recHits, caloTopology_) / sc.rawEnergy();
  ecalHLTV1_numberOfSubClusters = std::max(0, static_cast<int>(sc.clusters().size()) - 1);
  ecalHLTV1_clusterMaxDR = clusterMaxDR;
  ecalHLTV1_rawEnergy = sc.rawEnergy();
  ecalHLTV1_regressedEnergy = sc.energy();

  ecalHLTV1_gen_energy = gPs.energy();
  ecalHLTV1_gen_pt = gPs.pt();
  ecalHLTV1_gen_eta = gPs.eta();
  ecalHLTV1_gen_phi = gPs.phi();
  ecalHLTV1_gen_vz = gPs.vz();

  ecalHLTV1_eg_sigmaIEtaIEta = sigmaIEtaIEta;  
  ecalHLTV1_eg_sigmaIPhiIPhi = sigmaIPhiIPhi;  
  ecalHLTV1_isEB = true;
  ecalHLTV1_isEE = false;

  egRegDataEcalHLTV1Tree->Fill();
}

void EGammaNtuples::fillRegDataHGCALV1(const reco::SuperCluster& sc, const reco::GenParticle& gPs) {
  std::vector<float> eval(17, 0.);

  std::cout << "I'm crashing here" << std::endl;
  std::cout << sc.rawEnergy() << std::endl;
  auto ssCalc = hgcalShowerShapes_.createCalc(sc);
  auto pcaWidths = ssCalc.getPCAWidths(hgcalCylinderR_);
  auto energyHighestHits = ssCalc.getEnergyHighestHits(2);

  auto maxDRNonSeedClus = getMaxDRNonSeedCluster(sc);
  const float clusterMaxDR = maxDRNonSeedClus.first ? maxDRNonSeedClus.second : 999.;

  hgcalV1_rawEnergy = sc.rawEnergy();
  hgcalV1_regressedEnergy = sc.energy();
  hgcalV1_eta = sc.eta();
  hgcalV1_etaWidth = sc.etaWidth();
  hgcalV1_phiWidth = sc.phiWidth();
  hgcalV1_numberOfSubClusters = sc.clusters().size();
  hgcalV1_nrRecHits = sc.hitsAndFractions().size();
  hgcalV1_clusterMaxDR = clusterMaxDR;
  hgcalV1_dEtaSCToSeed = sc.eta() - sc.seed()->eta();
  hgcalV1_dPhiSCToSeed = reco::deltaPhi(sc.phi(), sc.seed()->phi());
  hgcalV1_firstRecHitEnergyFraction = energyHighestHits[0] / sc.rawEnergy();
  hgcalV1_secondRecHitEnergyFraction = energyHighestHits[1] / sc.rawEnergy();
  hgcalV1_sigma2uu = std::sqrt(pcaWidths.sigma2uu);
  hgcalV1_sigma2vv = std::sqrt(pcaWidths.sigma2vv);
  hgcalV1_sigma2ww = std::sqrt(pcaWidths.sigma2ww);
  hgcalV1_rvar = ssCalc.getRvar(hgcalCylinderR_, sc.rawEnergy());
  hgcalV1_seedEnergyFraction = sc.seed()->energy() / sc.rawEnergy();
  hgcalV1_nrHitsThreshold = nrHitsEB1GeV + nrHitsEE1GeV;

  hgcalV1_gen_energy = gPs.energy();
  hgcalV1_gen_pt = gPs.pt();
  hgcalV1_gen_eta = gPs.eta();
  hgcalV1_gen_phi = gPs.phi();
  hgcalV1_gen_vz = gPs.vz();

  hgcalV1_eg_sigmaIEtaIEta = 0;  
  hgcalV1_eg_sigmaIPhiIPhi = 0;  
  hgcalV1_isEB = false;
  hgcalV1_isEE = true;

  egRegDataHGCALV1Tree->Fill();
}

void EGammaNtuples::fillRegDataHGCALHLTV1(const reco::SuperCluster& sc, const reco::GenParticle& gPs) {
  const float clusterMaxDR = getMaxDRNonSeedCluster(sc).second;
  auto ssCalc = hgcalShowerShapes_.createCalc(sc);

  hgcalHLTV1_nrHitsThreshold = nrHitsEB1GeV + nrHitsEE1GeV;
  hgcalHLTV1_eta = sc.eta();
  hgcalHLTV1_phiWidth = sc.phiWidth();
  hgcalHLTV1_numberOfSubClusters = std::max(0, static_cast<int>(sc.clusters().size()) - 1);
  hgcalHLTV1_rvar = ssCalc.getRvar(hgcalCylinderR_);
  hgcalHLTV1_clusterMaxDR = clusterMaxDR;
  hgcalHLTV1_rawEnergy = sc.rawEnergy();
  hgcalHLTV1_regressedEnergy = sc.energy();

  hgcalHLTV1_gen_energy = gPs.energy();
  hgcalHLTV1_gen_pt = gPs.pt();
  hgcalHLTV1_gen_eta = gPs.eta();
  hgcalHLTV1_gen_phi = gPs.phi();
  hgcalHLTV1_gen_vz = gPs.vz();

  hgcalHLTV1_eg_sigmaIEtaIEta = 0;  
  hgcalHLTV1_eg_sigmaIPhiIPhi = 0;  
  hgcalHLTV1_isEB = false;
  hgcalHLTV1_isEE = true;

  egRegDataHGCALHLTV1Tree->Fill();
}

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


std::vector<DetId> EGammaNtuples::getNeighbors(const DetId seedId, const float dR, const CaloSubdetectorGeometry* subDetectorGeometry){

  const GlobalPoint pos = recHitTools_.getPosition(seedId);
  const int layer = recHitTools_.getLayerWithOffset(seedId);
  const int zside = recHitTools_.zside(seedId);
  const double dR2 = dR * dR;
  const double eta = pos.eta();
  const double phi = pos.phi();

  std::vector<DetId> neighbors;
  if (0.000001 < dR) {
    for (auto id : subDetectorGeometry->getValidDetIds()) {
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

  // Ecal takes 999. if no other cluster for maxDR2
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

  caloGeometry_ = &iSetup.getData(caloGeomToken_);
  caloTopology_ = &iSetup.getData(caloTopoToken_);

  recHitTools_.setGeometry(*caloGeometry_);
  int hgcalEEId = DetId::HGCalEE;
  
  fillRegDataEGRun3Tree(iEvent, iSetup);
  std::cout << "Finished Building Run3Trees" << std::endl;

  const auto& bScs = iEvent.get(scBarrelL1SeededToken_);
  const auto& eScs = iEvent.get(scHGCalL1SeededToken_);
  const auto& gPs =  iEvent.get(genParticleToken_);

  edm::Handle<EcalRecHitCollection> ebRecHitsHandle;  
  iEvent.getByToken(ebRecHitsToken_, ebRecHitsHandle);
  edm::Handle<HGCRecHitCollection> eeRecHitsHandle;
  iEvent.getByToken(eeRecHitsToken_, eeRecHitsHandle);
  edm::Handle<reco::PFRecHitCollection> pfRecHitsHGCALHandle;
  iEvent.getByToken(pfRecHitsHGCALToken_, pfRecHitsHGCALHandle);

  std::vector<reco::SuperCluster> scs(bScs.size() + eScs.size()); 
  std::merge(bScs.begin(), bScs.end(), eScs.begin(), eScs.end(), scs.begin());

  std::cout << "Initialize Showershapes" << std::endl;
  hgcalShowerShapes_.initPerSetup(iSetup);
  hgcalShowerShapes_.initPerEvent(*pfRecHitsHGCALHandle);
  std::cout << (*pfRecHitsHGCALHandle).size() << std::endl;
  nrHitsEB1GeV = countRecHits(ebRecHitsHandle, hitsEnergyThreshold_);
  nrHitsEE1GeV = countRecHits(eeRecHitsHandle, hitsEnergyThreshold_);
  iEvent.getByToken(verticesToken_, vertices_);

  if (gPs.size() == 2 && scs.size()==2){
    int i = 0;
    for (auto& sc: bScs){
      fillRegDataEcalV1(sc,*ebRecHitsHandle,gPs[i]);
      fillRegDataEcalHLTV1(sc,*ebRecHitsHandle,gPs[i]);
      std::cout << sc.rawEnergy() << ", " << sc.energy() << ", " << sc.correctedEnergy() << ", " << sc.correctedEnergyUncertainty()<< std::endl;
      i++;
    }
    std::cout << "Done with ECAL trees" << std::endl;
    i = 0;
    for (auto& sc: eScs){
      fillRegDataHGCALV1(sc, gPs[i]); 
      fillRegDataHGCALHLTV1(sc,gPs[i]);
      std::cout << sc.rawEnergy() << ", " << sc.energy() << ", " << sc.correctedEnergy() << ", " << sc.correctedEnergyUncertainty()<< std::endl;
      i++;
    }

    std::cout << "Done with HGCAL trees" << std::endl;

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

  // egRun3Tree
  egRun3Tree->Branch("runnr", &run_nr);
  egRun3Tree->Branch("lumiSec", &lumi_sec);
  egRun3Tree->Branch("eventnr", &event_nr);

  egRun3Tree->Branch("nrHitsEB1GeV", &nrHitsEB1GeV);
  egRun3Tree->Branch("nrHitsEE1GeV", &nrHitsEE1GeV);  
  
  egRun3Tree->Branch("eg_et", &eg_et);
  egRun3Tree->Branch("eg_energy", &eg_energy);
  egRun3Tree->Branch("eg_regressedEnergy", &sc_regressedEnergy);
  egRun3Tree->Branch("eg_rawEnergy", &sc_rawEnergy);
  egRun3Tree->Branch("eg_eta", &eg_eta);
  egRun3Tree->Branch("eg_phi", &eg_phi);
  egRun3Tree->Branch("eg_phiWidth", &sc_phiWidth);
  egRun3Tree->Branch("eg_nrClus", &sc_nrClus);
  egRun3Tree->Branch("eg_seedId", &sc_seedId);
  egRun3Tree->Branch("eg_seedDet", &sc_seedDet);
  egRun3Tree->Branch("eg_sigmaIEtaIEta", &eg_sigmaIEtaIEta);
  egRun3Tree->Branch("eg_sigmaIPhiIPhi", &eg_sigmaIPhiIPhi);
  egRun3Tree->Branch("eg_sigmaIEtaIEtaNoise", &eg_sigmaIEtaIEtaNoiseCleaned);
  egRun3Tree->Branch("eg_sigmaIPhiIPhiNoise", &eg_sigmaIPhiIPhiNoiseCleaned);
  egRun3Tree->Branch("eg_clusterMaxDR", &sc_clusterMaxDR);
  egRun3Tree->Branch("eg_r9Frac", &sc_r9Frac);
  egRun3Tree->Branch("eg_r9Full", &sc_r9Full);
  egRun3Tree->Branch("eg_isEB", &sc_isEB);
  egRun3Tree->Branch("eg_isEE", &sc_isEE);

  egRun3Tree->Branch("eg_gen_energy", &gen_energy);
  egRun3Tree->Branch("eg_gen_pt", &gen_pt);
  egRun3Tree->Branch("eg_gen_eta", &gen_eta);
  egRun3Tree->Branch("eg_gen_phi", &gen_phi);
  egRun3Tree->Branch("eg_gen_vz", &gen_vz);

  // egRegDataEcalV1Tree
  egRegDataEcalV1Tree->Branch("nrVert", &ecalV1_nrVert);
  egRegDataEcalV1Tree->Branch("rawEnergy", &ecalV1_rawEnergy);
  egRegDataEcalV1Tree->Branch("regressedEnergy", &ecalV1_regressedEnergy);
  egRegDataEcalV1Tree->Branch("etaWidth", &ecalV1_etaWidth);
  egRegDataEcalV1Tree->Branch("phiWidth", &ecalV1_phiWidth);
  egRegDataEcalV1Tree->Branch("rvar", &ecalV1_rvar);
  egRegDataEcalV1Tree->Branch("energySeedCluster", &ecalV1_energySeedCluster);
  egRegDataEcalV1Tree->Branch("energyFirstRecHit", &ecalV1_energyFirstRecHit);
  egRegDataEcalV1Tree->Branch("energySecondRecHit", &ecalV1_energySecondRecHit);
  egRegDataEcalV1Tree->Branch("eLeftRightDiffSumRatio", &ecalV1_eLeftRightDiffSumRatio);
  egRegDataEcalV1Tree->Branch("eTopBottomDiffSumRatio", &ecalV1_eTopBottomDiffSumRatio);
  egRegDataEcalV1Tree->Branch("sigmaIEtaIEta", &ecalV1_sigmaIEtaIEta);
  egRegDataEcalV1Tree->Branch("sigmaIEtaIPhi", &ecalV1_sigmaIEtaIPhi);
  egRegDataEcalV1Tree->Branch("sigmaIPhiIPhi", &ecalV1_sigmaIPhiIPhi);
  egRegDataEcalV1Tree->Branch("numberOfSubClusters", &ecalV1_numberOfSubClusters);
  egRegDataEcalV1Tree->Branch("clusterMaxDR", &ecalV1_clusterMaxDR);
  egRegDataEcalV1Tree->Branch("clusterMaxDRDPhi", &ecalV1_clusterMaxDRDPhi);
  egRegDataEcalV1Tree->Branch("clusterMaxDRDEta", &ecalV1_clusterMaxDRDEta);
  egRegDataEcalV1Tree->Branch("clusterMaxDREnergyFraction", &ecalV1_clusterMaxDREnergyFraction);
  egRegDataEcalV1Tree->Branch("subCluster1RawEnergyFraction", &ecalV1_subCluster1RawEnergyFraction);
  egRegDataEcalV1Tree->Branch("subCluster2RawEnergyFraction", &ecalV1_subCluster2RawEnergyFraction);
  egRegDataEcalV1Tree->Branch("subCluster3RawEnergyFraction", &ecalV1_subCluster3RawEnergyFraction);
  egRegDataEcalV1Tree->Branch("subCluster1DPhiToSeed", &ecalV1_subCluster1DPhiToSeed);
  egRegDataEcalV1Tree->Branch("subCluster2DPhiToSeed", &ecalV1_subCluster2DPhiToSeed);
  egRegDataEcalV1Tree->Branch("subCluster3DPhiToSeed", &ecalV1_subCluster3DPhiToSeed);
  egRegDataEcalV1Tree->Branch("subCluster1DEtaToSeed", &ecalV1_subCluster1DEtaToSeed);
  egRegDataEcalV1Tree->Branch("subCluster2DEtaToSeed", &ecalV1_subCluster2DEtaToSeed);
  egRegDataEcalV1Tree->Branch("subCluster3DEtaToSeed", &ecalV1_subCluster3DEtaToSeed);
  egRegDataEcalV1Tree->Branch("seedClusteriEtaOrX", &ecalV1_seedClusteriEtaOrX);
  egRegDataEcalV1Tree->Branch("seedClusteriPhiOrY", &ecalV1_seedClusteriEtaOrX);
  egRegDataEcalV1Tree->Branch("seedClusterEta", &ecalV1_seedClusterEta);
  egRegDataEcalV1Tree->Branch("eg_gen_energy", &ecalV1_gen_energy);
  egRegDataEcalV1Tree->Branch("eg_gen_pt", &ecalV1_gen_pt);
  egRegDataEcalV1Tree->Branch("eg_gen_eta", &ecalV1_gen_eta);
  egRegDataEcalV1Tree->Branch("eg_gen_phi", &ecalV1_gen_phi);
  egRegDataEcalV1Tree->Branch("eg_gen_vz", &ecalV1_gen_vz);

  egRegDataEcalV1Tree->Branch("eg_isEB", &ecalV1_sc_isEB);
  egRegDataEcalV1Tree->Branch("eg_isEE", &ecalV1_sc_isEE);
  egRegDataEcalV1Tree->Branch("eg_sigmaIEtaIEta", &ecalV1_eg_sigmaIEtaIEta);
  egRegDataEcalV1Tree->Branch("eg_sigmaIPhiIPhi", &ecalV1_eg_sigmaIPhiIPhi);

  // egRegDataEcalHLTV1Tree
  egRegDataEcalHLTV1Tree->Branch("nrHitsThreshold", &ecalHLTV1_nrHitsThreshold);
  egRegDataEcalHLTV1Tree->Branch("eta", &ecalHLTV1_eta);
  egRegDataEcalHLTV1Tree->Branch("rawEnergy", &ecalHLTV1_rawEnergy);
  egRegDataEcalHLTV1Tree->Branch("regressedEnergy", &ecalHLTV1_regressedEnergy);
  egRegDataEcalHLTV1Tree->Branch("phiWidth", &ecalHLTV1_phiWidth);
  egRegDataEcalHLTV1Tree->Branch("rvar", &ecalHLTV1_rvar);
  egRegDataEcalHLTV1Tree->Branch("numberOfSubClusters", &ecalHLTV1_numberOfSubClusters);
  egRegDataEcalHLTV1Tree->Branch("clusterMaxDR", &ecalHLTV1_clusterMaxDR);
  egRegDataEcalHLTV1Tree->Branch("eg_gen_energy", &ecalHLTV1_gen_energy);
  egRegDataEcalHLTV1Tree->Branch("eg_gen_pt", &ecalHLTV1_gen_pt);
  egRegDataEcalHLTV1Tree->Branch("eg_gen_eta", &ecalHLTV1_gen_eta);
  egRegDataEcalHLTV1Tree->Branch("eg_gen_phi", &ecalHLTV1_gen_phi);
  egRegDataEcalHLTV1Tree->Branch("eg_gen_vz", &ecalHLTV1_gen_vz);

  egRegDataEcalHLTV1Tree->Branch("eg_isEB", &ecalHLTV1_sc_isEB);
  egRegDataEcalHLTV1Tree->Branch("eg_isEE", &ecalHLTV1_sc_isEE);
  egRegDataEcalHLTV1Tree->Branch("eg_sigmaIEtaIEta", &ecalHLTV1_eg_sigmaIEtaIEta);
  egRegDataEcalHLTV1Tree->Branch("eg_sigmaIPhiIPhi", &ecalHLTV1_eg_sigmaIPhiIPhi);

  // egRegDataHGCALV1Tree
  egRegDataHGCALV1Tree->Branch("rawEnergy", &hgcalV1_rawEnergy);
  egRegDataHGCALV1Tree->Branch("regressedEnergy", &hgcalV1_regressedEnergy);
  egRegDataHGCALV1Tree->Branch("eta", &hgcalV1_eta);
  egRegDataHGCALV1Tree->Branch("etaWidth", &hgcalV1_etaWidth);
  egRegDataHGCALV1Tree->Branch("phiWidth", &hgcalV1_phiWidth);
  egRegDataHGCALV1Tree->Branch("numberOfSubClusters", &hgcalV1_numberOfSubClusters);
  egRegDataHGCALV1Tree->Branch("nrRecHits", &hgcalV1_nrRecHits);
  egRegDataHGCALV1Tree->Branch("clusterMaxDR", &hgcalV1_clusterMaxDR);
  egRegDataHGCALV1Tree->Branch("dEtaSCToSeed", &hgcalV1_dEtaSCToSeed);
  egRegDataHGCALV1Tree->Branch("dPhiSCToSeed", &hgcalV1_dPhiSCToSeed);
  egRegDataHGCALV1Tree->Branch("firstRecHitEnergyFraction", &hgcalV1_firstRecHitEnergyFraction);
  egRegDataHGCALV1Tree->Branch("secondRecHitEnergyFraction", &hgcalV1_secondRecHitEnergyFraction);
  egRegDataHGCALV1Tree->Branch("sigma2uu", &hgcalV1_sigma2uu);
  egRegDataHGCALV1Tree->Branch("sigma2vv", &hgcalV1_sigma2vv);
  egRegDataHGCALV1Tree->Branch("sigma2ww", &hgcalV1_sigma2ww);
  egRegDataHGCALV1Tree->Branch("rvar", &hgcalV1_rvar);
  egRegDataHGCALV1Tree->Branch("seedEnergyFraction", &hgcalV1_seedEnergyFraction);
  egRegDataHGCALV1Tree->Branch("nrHitsThreshold", &hgcalV1_nrHitsThreshold);
  egRegDataHGCALV1Tree->Branch("eg_gen_energy", &hgcalV1_gen_energy);
  egRegDataHGCALV1Tree->Branch("eg_gen_pt", &hgcalV1_gen_pt);
  egRegDataHGCALV1Tree->Branch("eg_gen_eta", &hgcalV1_gen_eta);
  egRegDataHGCALV1Tree->Branch("eg_gen_phi", &hgcalV1_gen_phi);
  egRegDataHGCALV1Tree->Branch("eg_gen_vz", &hgcalV1_gen_vz);

  egRegDataHGCALV1Tree->Branch("eg_isEB", &hgcalV1_sc_isEB);
  egRegDataHGCALV1Tree->Branch("eg_isEE", &hgcalV1_sc_isEE);
  egRegDataHGCALV1Tree->Branch("eg_sigmaIEtaIEta", &hgcalV1_eg_sigmaIEtaIEta);
  egRegDataHGCALV1Tree->Branch("eg_sigmaIPhiIPhi", &hgcalV1_eg_sigmaIPhiIPhi);


 // egRegDataHGCALHLTV1Tree
  egRegDataHGCALHLTV1Tree->Branch("rawEnergy", &hgcalHLTV1_rawEnergy);
  egRegDataHGCALHLTV1Tree->Branch("regressedEnergy", &hgcalHLTV1_regressedEnergy);
  egRegDataHGCALHLTV1Tree->Branch("eta", &hgcalHLTV1_eta);
  egRegDataHGCALHLTV1Tree->Branch("phiWidth", &hgcalHLTV1_phiWidth);
  egRegDataHGCALHLTV1Tree->Branch("numberOfSubClusters", &hgcalHLTV1_numberOfSubClusters);
  egRegDataHGCALHLTV1Tree->Branch("clusterMaxDR", &hgcalHLTV1_clusterMaxDR);
  egRegDataHGCALHLTV1Tree->Branch("rvar", &hgcalHLTV1_rvar);
  egRegDataHGCALHLTV1Tree->Branch("nrHitsThreshold", &hgcalHLTV1_nrHitsThreshold);
  egRegDataHGCALHLTV1Tree->Branch("eg_gen_energy", &hgcalHLTV1_gen_energy);
  egRegDataHGCALHLTV1Tree->Branch("eg_gen_pt", &hgcalHLTV1_gen_pt);
  egRegDataHGCALHLTV1Tree->Branch("eg_gen_eta", &hgcalHLTV1_gen_eta);
  egRegDataHGCALHLTV1Tree->Branch("eg_gen_phi", &hgcalHLTV1_gen_phi);
  egRegDataHGCALHLTV1Tree->Branch("eg_gen_vz", &hgcalHLTV1_gen_vz);
  egRegDataHGCALV1Tree->Branch("eg_isEB",&hgcal_isEB);

  egRegDataHGCALHLTV1Tree->Branch("eg_isEB", &hgcalHLTV1_sc_isEB);
  egRegDataHGCALHLTV1Tree->Branch("eg_isEE", &hgcalHLTV1_sc_isEE);
  egRegDataHGCALHLTV1Tree->Branch("eg_sigmaIEtaIEta", &hgcalHLTV1_eg_sigmaIEtaIEta);
  egRegDataHGCALHLTV1Tree->Branch("eg_sigmaIPhiIPhi", &hgcalHLTV1_eg_sigmaIPhiIPhi);

  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void EGammaNtuples::endJob() {
  // please remove this method if not needed
  newfile->cd();
  egRun3Tree->Write();
  egRegDataEcalV1Tree->Write();
  egRegDataEcalHLTV1Tree->Write();
  egRegDataHGCALV1Tree->Write();
  egRegDataHGCALHLTV1Tree->Write();
  newfile->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EGammaNtuples::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genParticles",edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("eGammaObjects", edm::InputTag("hltEgammaHLTExtra",""));
  desc.add<edm::InputTag>("eGammaObjectsUnSeeded",edm::InputTag("hltEgammaHLTExtra","Unseeded"));
  desc.add<edm::InputTag>("scBarrelL1Seeded", edm::InputTag("hltParticleFlowSuperClusterECALL1Seeded","hltParticleFlowSuperClusterECALBarrel"));
  desc.add<edm::InputTag>("scHGCalL1Seeded", edm::InputTag("particleFlowSuperClusterHGCalFromTICLL1Seeded",""));
  desc.add<edm::InputTag>("ebRecHits", edm::InputTag("hltEgammaHLTExtra","EcalRecHitsEB"));
  desc.add<edm::InputTag>("eeRecHits", edm::InputTag("HGCalRecHit","HGCEERecHits"));
  desc.add<edm::InputTag>("sigmaIEtaIEta", edm::InputTag("hltEgammaClusterShapeL1Seeded","sigmaIEtaIEta5x5"));
  desc.add<edm::InputTag>("sigmaIPhiIPhi", edm::InputTag("hltEgammaClusterShapeL1Seeded","sigmaIPhiIPhi5x5"));
  desc.add<edm::InputTag>("sigmaIEtaIEtaNoiseCleaned", edm::InputTag("hltEgammaClusterShapeL1Seeded","sigmaIEtaIEta5x5NoiseCleaned"));
  desc.add<edm::InputTag>("sigmaIPhiIPhiNoiseCleaned", edm::InputTag("hltEgammaClusterShapeL1Seeded","sigmaIPhiIPhi5x5NoiseCleaned"));
  desc.add<edm::InputTag>("nrHitsEB1GeV", edm::InputTag("hltEgammaHLTExtra","countEcalRecHitsEcalRecHitsEBThres1GeV"));
  desc.add<edm::InputTag>("nrHitsEE1GeV", edm::InputTag("hltEgammaHLTExtra","countEcalRecHitsHGCalRecHitsThres1GeV"));
  desc.add<edm::InputTag>("eGammaCandidates", edm::InputTag("hltEgammaCandidatesL1Seeded"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("pfHGCALRecHits", edm::InputTag("particleFlowRecHitHGC"));
  descriptions.add("EGammaNtuples",desc);


  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EGammaNtuples);


