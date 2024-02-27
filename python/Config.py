import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/work/m/mmatthew/private/egamma/CMSSW_13_1_0/src/tutorialDir/DoubleElectron_FlatPt-1To100-gun/hlt_IDEAL.root'
            #'file:/afs/cern.ch/work/m/mmatthew/private/deleteMe/CMSSW_13_1_0/src/Phase2_HLT.root'    
                )
                            )

process.demo = cms.EDAnalyzer('EGammaNtuples',
                              )

process.p = cms.Path(process.demo)
