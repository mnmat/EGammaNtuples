from PlayGround.EGammaNtuples.EGammaNtuples_cfi import EGammaNtuples
import FWCore.ParameterSet.Config as cms

def customiseHLTForEGammaNtuples(process):

    process.EGammaNtuples = EGammaNtuples.clone()

    process.FEVTDEBUGHLToutput_step = cms.EndPath(
        process.FEVTDEBUGHLToutput + process.EGammaNtuples)
    return process