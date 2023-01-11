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
                          #                   "HLT_Dimuon16_Jpsi",                        # inclusive dimuon jpsi
                          #                   "HLT_Dimuon10_Jpsi_Barrel",                 # inclusive dimuon jpsi
                          #                   "HLT_DoubleMu4_JpsiTrk_Displaced",          # displaced jpistrk 
                          #                   "HLT_Dimuon20_JpsiHLT_Dimuon20_Jpsi",       # 
                          #                   "HLT_DoubleMu4_3_Jpsi_Displaced",           # prescaled
                          #                   "HLT_DoubleMu4_3_Jpsi",                     # prescaled
                          #                   "HLT_DoubleMu4_Jpsi_Displaced",             # prescaled
                          #                   "HLT_Dimuon13_PsiPrime",                    # inclusive dimuon psi2s
                          #                   "HLT_Dimuon8_PsiPrime_Barrel",              # inclusive dimuon psi2s
                          #                   "HLT_DoubleMu4_PsiPrimeTrk_Displaced",      # displaced psi2s trk 
                          #                   "HLT_Dimuon0_Jpsi_Muon",                    # triple-mu (jpsi + muon)
                          #                   "HLT_Dimuon6_Jpsi_NoVertexing",             # jpsi no vertex very prescaled
                          #                   "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing",    # jpsi no vertexing paths
                          #                   "HLT_Dimuon0er16_Jpsi_NoVertexing"          # jpsi no vertexing paths
                          #                    ),
)

trgTables = cms.Sequence(trgTable)



