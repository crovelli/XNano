import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

trgTable = cms.EDProducer( "TrgBitTableProducer",
                          hltresults = cms.InputTag("TriggerResults::HLT"),

                          # interesting paths for 2017 and 2018
                          paths      = cms.vstring(
                                             "HLT_Dimuon25_Jpsi",                        # inclusive dimuon jpsi
                                             "HLT_Dimuon20_Jpsi_Barrel_Seagulls",        # inclusive dimuon jpsi
                                             "HLT_DoubleMu4_JpsiTrk_Displaced",          # displaced jpistrk or jpsi trktrk
                                             "HLT_DoubleMu4_JpsiTrkTrk_Displaced",       # displaced jpistrk or jpsi trktrk
                                             "HLT_DoubleMu4_3_Jpsi_Displaced",           # prescaled for 2017
                                             "HLT_DoubleMu4_3_Jpsi",                     # prescaled for 2018  
                                             "HLT_DoubleMu4_Jpsi_Displaced",             # prescaled
                                             "HLT_Dimuon18_PsiPrime",                    # inclusive dimuon psi2s
                                             "HLT_Dimuon10_PsiPrime_Barrel_Seagulls",    # inclusive dimuon psi2s
                                             "HLT_DoubleMu4_PsiPrimeTrk_Displaced",      # displaced psi2s trk 
                                             "HLT_Dimuon0_Jpsi3p5_Muon2",                # triple-mu (jpsi + muon)
                                             "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi",       # jpsi + 2 trkmu (phi->mumu)
                                             "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi",        # jpsi+2 trk(phi->KK) for 2017 
                                             "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05",    # jpsi+2 trk(phi->KK) for 2018
                                             "HLT_DoubleMu4_3_Bs"                        # Livia 
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



