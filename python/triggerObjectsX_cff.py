import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.triggerObjects_cff import *

triggerObjectXTable = cms.EDProducer("TriggerObjectTableProducerX",
    name= cms.string("TrigObj"),
    src = cms.InputTag("unpackedPatTrigger"),
    selections = cms.VPSet(
        cms.PSet(
            name = cms.string("Muon"),
            id = cms.int32(13),
            sel = cms.string("type(83) && pt > 5 && coll('hltIterL3MuonCandidates')"), 
            qualityBits = cms.string("filter('hlt')"), qualityBitsDoc = cms.string("1 = Muon filters"),
        ),
    ),
)

triggerObjectXTables = cms.Sequence( unpackedPatTrigger + triggerObjectXTable )
