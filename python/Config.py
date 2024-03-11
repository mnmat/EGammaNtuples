import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/work/m/mmatthew/private/egamma/CMSSW_13_1_0/src/Phase2_HLT.root'
            #'file:/afs/cern.ch/work/m/mmatthew/private/egamma/CMSSW_13_1_0/src/tutorialDir/DoubleElectron_FlatPt-1To100-gun/hlt_IDEAL.root'
            #'file:/afs/cern.ch/work/m/mmatthew/private/deleteMe/CMSSW_13_1_0/src/Phase2_HLT.root'    
                )
            )

process.demo = cms.EDAnalyzer('EGammaNtuples',
    genParticles = cms.untracked.InputTag('genParticles'),
    eGammaObjects = cms.untracked.InputTag('hltEgammaHLTExtra',""),
    eGammaObjectsUnSeeded = cms.untracked.InputTag('hltEgammaHLTExtra',"Unseeded"),
    scBarrelL1Seeded = cms.untracked.InputTag('hltParticleFlowSuperClusterECALL1Seeded','hltParticleFlowSuperClusterECALBarrel'),
    scHGCalL1Seeded = cms.untracked.InputTag('particleFlowSuperClusterHGCalFromTICLL1Seeded',''),
    ebRecHits = cms.untracked.InputTag('hltEgammaHLTExtra',"EcalRecHitsEB"),
    eeRecHits = cms.untracked.InputTag('HGCalRecHit',"HGCEERecHits"),
    sigmaIEtaIEta = cms.untracked.InputTag('hltEgammaClusterShapeL1Seeded','sigmaIEtaIEta5x5'),
    sigmaIPhiIPhi = cms.untracked.InputTag('hltEgammaClusterShapeL1Seeded','sigmaIPhiIPhi5x5'),
    sigmaIEtaIEtaNoiseCleaned = cms.untracked.InputTag('hltEgammaClusterShapeL1Seeded','sigmaIEtaIEta5x5NoiseCleaned'),
    sigmaIPhiIPhiNoiseCleaned = cms.untracked.InputTag('hltEgammaClusterShapeL1Seeded','sigmaIPhiIPhi5x5NoiseCleaned'),
    nrHitsEB1GeV = cms.untracked.InputTag('hltEgammaHLTExtra',"countEcalRecHitsEcalRecHitsEBThres1GeV"),
    nrHitsEE1GeV = cms.untracked.InputTag('hltEgammaHLTExtra',"countEcalRecHitsHGCalRecHitsThres1GeV"),
    eGammaCandidates = cms.untracked.InputTag('hltEgammaCandidatesL1Seeded')

                              )

process.p = cms.Path(process.demo)
