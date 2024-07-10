#ifndef RecoEgamma_EgammaTools_HGCalShowerShapeHelper_Dummy_h
#define RecoEgamma_EgammaTools_HGCalShowerShapeHelper_Dummy_h

// system include files
#include <map>
#include <memory>
#include <utility>
#include <vector>

// external include files
#include <Math/Vector3Dfwd.h>

// CMSSW include files
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/Transition.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

class HGCalShowerShapeHelper_Dummy {
  // Good to filter/compute/store this stuff beforehand as they are common to the shower shape variables.
  // No point in filtering, computing layer-wise centroids, etc. for each variable again and again.
  // Once intitialized, one can the calculate different variables one after another for a given object.
  // This is all handled by ShowerShapeCalc class which caches the layer-wise centroids and other
  // heavy variables for an object + set of cuts
  // It was changed to this approach so that we could use this in constant functions

  // In principle should consider the HGCalHSi and HGCalHSc hits (leakage) also.
  // Can have subdetector dependent thresholds and layer selection.
  // To be implemented.

public:
  static const double kLDWaferCellSize_;
  static const double kHDWaferCellSize_;

  struct ShowerWidths {
    double sigma2xx;
    double sigma2yy;
    double sigma2zz;

    double sigma2xy;
    double sigma2yz;
    double sigma2zx;

    double sigma2uu;
    double sigma2vv;
    double sigma2ww;

    ShowerWidths()
        : sigma2xx(0.0),
          sigma2yy(0.0),
          sigma2zz(0.0),
          sigma2xy(0.0),
          sigma2yz(0.0),
          sigma2zx(0.0),
          sigma2uu(0.0),
          sigma2vv(0.0),
          sigma2ww(0.0) {}
  };

  class ShowerShapeCalc {
  public:
    ShowerShapeCalc(std::shared_ptr<const hgcal::RecHitTools> recHitTools,
                    std::map<DetId, const HGCRecHit*>> hitMap,
                    const std::vector<std::pair<DetId, float> > &hitsAndFracs,
                    const double rawEnergy,
                    const double minHitE = 0,
                    const double minHitET = 0,
                    const int minLayer = 1,
                    const int maxLayer = -1,
                    DetId::Detector subDet = DetId::HGCalEE);

    double getCellSize(DetId detId) const;

    // Compute Rvar in a cylinder around the layer centroids
    double getRvar(double cylinderR, bool useFractions = true, bool useCellSize = true) const;

    // Compute PCA widths around the layer centroids
    ShowerWidths getPCAWidths(double cylinderR, bool useFractions = false) const;

    std::vector<double> getEnergyHighestHits(unsigned int nrHits, bool useFractions = true) const;

  private:
    void setFilteredHitsAndFractions(const std::vector<std::pair<DetId, float> > &hitsAndFracs);
    void setLayerWiseInfo();

    std::shared_ptr<const hgcal::RecHitTools> recHitTools_;
    //std::shared_ptr<const std::unordered_map<uint32_t, const reco::PFRecHit *> > pfRecHitPtrMap_;
    std::map<DetId, const HGCRecHit*>> hitMap;

    double rawEnergy_;

    double minHitE_;
    double minHitET_;
    double minHitET2_;
    int minLayer_;
    int maxLayer_;
    int nLayer_;
    DetId::Detector subDet_;

    std::vector<std::pair<DetId, float> > hitsAndFracs_;
    std::vector<double> hitEnergies_;
    std::vector<double> hitEnergiesWithFracs_;

    ROOT::Math::XYZVector centroid_;
    std::vector<double> layerEnergies_;
    std::vector<ROOT::Math::XYZVector> layerCentroids_;
  };

  HGCalShowerShapeHelper_Dummy();
  HGCalShowerShapeHelper_Dummy(edm::ConsumesCollector &&sumes);
  ~HGCalShowerShapeHelper_Dummy() = default;
  HGCalShowerShapeHelper_Dummy(const HGCalShowerShapeHelper_Dummy &rhs) = delete;
  HGCalShowerShapeHelper_Dummy(const HGCalShowerShapeHelper_Dummy &&rhs) = delete;
  HGCalShowerShapeHelper_Dummy &operator=(const HGCalShowerShapeHelper_Dummy &rhs) = delete;
  HGCalShowerShapeHelper_Dummy &operator=(const HGCalShowerShapeHelper_Dummy &&rhs) = delete;

  template <edm::Transition tr = edm::Transition::Event>
  void setTokens(edm::ConsumesCollector consumesCollector) {
    caloGeometryToken_ = consumesCollector.esConsumes<CaloGeometry, CaloGeometryRecord, tr>();
  }

  void initPerSetup(const edm::EventSetup &iSetup);
  void initPerEvent(const HGCRecHitCollection& rechitsEE, 
    const HGCRecHitCollection& rechitsFH,
    const HGCRecHitCollection& rechitsBH);
  void initPerEvent(const edm::EventSetup &iSetup, const std::vector<reco::PFRecHit> &recHits);

  HGCalShowerShapeHelper_Dummy::ShowerShapeCalc createCalc(const std::vector<std::pair<DetId, float> > &hitsAndFracs,
                                                     double rawEnergy,
                                                     double minHitE = 0,
                                                     double minHitET = 0,
                                                     int minLayer = 1,
                                                     int maxLayer = -1,
                                                     DetId::Detector subDet = DetId::HGCalEE) const;
  HGCalShowerShapeHelper_Dummy::ShowerShapeCalc createCalc(const reco::SuperCluster &sc,
                                                     double minHitE = 0,
                                                     double minHitET = 0,
                                                     int minLayer = 1,
                                                     int maxLayer = -1,
                                                     DetId::Detector subDet = DetId::HGCalEE) const {
    return createCalc(sc.hitsAndFractions(), sc.rawEnergy(), minHitE, minHitET, minLayer, maxLayer, subDet);
  }

private:
  // void setPFRecHitPtrMap(const std::vector<reco::PFRecHit> &recHits);
  
  void fillHitMap(std::map<DetId, std::pair<const HGCRecHit*,float>>& hitMap, 
      const HGCRecHitCollection& rechitsEE, 
      const HGCRecHitCollection& rechitsFH,
      const HGCRecHitCollection& rechitsBH,
      const std::vector<float>& lcmask) const;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  std::shared_ptr<hgcal::RecHitTools> recHitTools_;

  std::map<DetId, const HGCRecHit*>> hitMap;
  // std::shared_ptr<std::unordered_map<uint32_t, const reco::PFRecHit *>> pfRecHitPtrMap_;
};

#endif