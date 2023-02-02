import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

trgTable = cms.EDProducer( "TrgBitTableProducer",
                          hltresults = cms.InputTag("TriggerResults::HLT"),

                          # interesting paths for 2017 and 2018
                          paths      = cms.vstring(
                                             "HLT_Dimuon25_Jpsi",                        # inclusive dimuon jpsi
                                             "HLT_DoubleMu4_JpsiTrk_Displaced"           # displaced jpistrk 
                                              ),

                          # interesting paths for 2016 
                          # paths      = cms.vstring(
                          #                   "HLT_Dimuon20_Jpsi",                        # inclusive dimuon jpsi. There is also Dimuon16, but only for part of the dataset
                          #                   "HLT_DoubleMu4_JpsiTrk_Displaced",          # displaced jpistrk 
                          #                    ),
)

trgTables = cms.Sequence(trgTable)



