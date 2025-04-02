from EgRegresNtuples.EGammaNtuples.EGammaNtuples_cfi import EGammaNtuples
import FWCore.ParameterSet.Config as cms

def customiseHLTForEGammaNtuples(process):

    process.EGammaNtuples = EGammaNtuples.clone(
        pType = "ele",
    )

    from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
    ticl_v5.toModify(process.EGammaNtuples, scHGCalL1Seeded = cms.InputTag('hltTiclEGammaSuperClusterProducerUnseeded'))

    process.FEVTDEBUGHLToutput_step = cms.EndPath(
        process.FEVTDEBUGHLToutput + process.EGammaNtuples)
    return process
