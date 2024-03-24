from PlayGround.EGammaNtuples.EGammaNtuples_cfi import EGammaNtuples
import FWCore.ParameterSet.Config as cms

def customiseHLTForEGammaNtuples(process):

    process.EGammaNtuples = EGammaNtuples.clone()

    process.TFileService = cms.Service("TFileService",
                                       fileName=cms.string("output.root")
                                       )
    process.FEVTDEBUGHLToutput_step = cms.EndPath(
        process.FEVTDEBUGHLToutput + process.EGammaNtuples)
    return process