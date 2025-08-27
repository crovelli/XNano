import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

trgTable = cms.EDProducer( "TrgBitTableProducer",
                          hltresults = cms.InputTag("TriggerResults::HLT"),

                          # interesting paths for 2022
                          paths      = cms.vstring(
                                             "HLT_DoubleMu4_3_LowMass"
                                              ),
)

trgTables = cms.Sequence(trgTable)



