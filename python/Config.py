import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'root://xrootd-cms.infn.it//store/user/mmatthew/crab_projects/egamma/DoublePhoton_FlatPt-1To100-gun/crab_photon/240722_083350/0000/Phase2_IDEAL_HLT_14.root'
            #'file:/afs/cern.ch/work/m/mmatthew/private/egamma/CMSSW_13_1_0/src/tutorialDir/DoubleElectron_FlatPt-1To100-gun/hlt_IDEAL.root'
            #'file:/afs/cern.ch/work/m/mmatthew/private/deleteMe/CMSSW_13_1_0/src/Phase2_HLT.root'    
                )
            )


from PlayGround.EGammaNtuples.EGammaNtuples_cfi import EGammaNtuples
process.demo = EGammaNtuples.clone(pType = "pho")


process.p = cms.Path(process.demo)
