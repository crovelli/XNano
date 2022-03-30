import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

Path2016=["HLT_Dimuon16_Jpsi","HLT_Dimuon10_Jpsi_Barrel","HLT_DoubleMu4_JpsiTrk_Displaced","HLT_Dimuon20_Jpsi","HLT_DoubleMu4_3_Jpsi_Displaced","HLT_DoubleMu4_3_Jpsi","HLT_DoubleMu4_Jpsi_Displaced","HLT_Dimuon13_PsiPrime","HLT_Dimuon8_PsiPrime_Barrel","HLT_DoubleMu4_PsiPrimeTrk_Displaced","HLT_Dimuon0_Jpsi_Muon","HLT_Dimuon6_Jpsi_NoVertexing","HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing","HLT_Dimuon0er16_Jpsi_NoVertexing"]

Path2017=["HLT_Dimuon25_Jpsi","HLT_Dimuon20_Jpsi_Barrel_Seagulls","HLT_DoubleMu4_JpsiTrk_Displaced","HLT_DoubleMu4_JpsiTrkTrk_Displaced","HLT_DoubleMu4_3_Jpsi_Displaced","HLT_DoubleMu4_3_Jpsi","HLT_DoubleMu4_Jpsi_Displaced","HLT_Dimuon18_PsiPrime","HLT_Dimuon10_PsiPrime_Barrel_Seagulls","HLT_DoubleMu4_PsiPrimeTrk_Displaced","HLT_Dimuon0_Jpsi3p5_Muon2","HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi","HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi","HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05"]

Path=Path2017

muonTrgSelector = cms.EDProducer("MuonTriggerSelector",
                                 muonCollection = cms.InputTag("slimmedMuons"), 
                                 bits = cms.InputTag("TriggerResults","","HLT"),
                                 objects = cms.InputTag("slimmedPatTrigger"),
                                 
                                 ## trigger match
                                 drForTriggerMatch = cms.double(0.1),        # to be tuned

                                 ## for the output selected collection 
                                 dzForCleaning_wrtTrgMuon = cms.double(-1.),     # chiara: to be used when trigger match is fixed
                                 ptMin = cms.double(0.5),                            
                                 absEtaMax = cms.double(2.4),

                                 HLTPaths=cms.vstring(Path)
                             )

# at least 1 triggering muon
countTrgMuons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),       
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("muonTrgSelector", "trgMuons")
)

# muons selection
muonXTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("muonTrgSelector:SelectedMuons"),
    cut = cms.string(""), 
    name = cms.string("Muon"),
    doc  = cms.string("slimmedMuons for X analysis after basic selection"),
    singleton = cms.bool(False),         
    extension = cms.bool(False),         
    variables = cms.PSet(CandVars,
        ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
        vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
        vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
        vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6),
        pfRelIso03_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
        pfRelIso04_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.4, total (deltaBeta corrections)"),
        isGlobal = Var("userInt('isGlobal')",bool,doc="muon is global muon"),
        softId = Var("passed('SoftCutBasedId')",bool,doc="soft cut-based ID"), 

        highQualityTrack = Var("userInt('highQualityTrack')",int,doc="muon track is high quality"),                 

        isTriggering = Var("userInt('isTriggering')", int,doc="flag the reco muon is also triggering"),

        fired_HLT_Dimuon25_Jpsi = Var("userInt('HLT_Dimuon25_Jpsi')",int,doc="reco muon fired this trigger"),
        fired_HLT_Dimuon20_Jpsi_Barrel_Seagulls = Var("userInt('HLT_Dimuon20_Jpsi_Barrel_Seagulls')",int,doc="reco muon fired this trigger"),
        fired_HLT_DoubleMu4_JpsiTrk_Displaced = Var("userInt('HLT_DoubleMu4_JpsiTrk_Displaced')",int,doc="reco muon fired this trigger"),
        fired_HLT_DoubleMu4_JpsiTrkTrk_Displaced = Var("userInt('HLT_DoubleMu4_JpsiTrkTrk_Displaced')",int,doc="reco muon fired this trigger"),
        fired_HLT_DoubleMu4_3_Jpsi_Displaced = Var("userInt('HLT_DoubleMu4_3_Jpsi_Displaced')",int,doc="reco muon fired this trigger"),
        fired_HLT_DoubleMu4_3_Jpsi = Var("userInt('HLT_DoubleMu4_3_Jpsi')",int,doc="reco muon fired this trigger"),
        fired_HLT_DoubleMu4_Jpsi_Displaced = Var("userInt('HLT_DoubleMu4_Jpsi_Displaced')",int,doc="reco muon fired this trigger"),
        fired_HLT_Dimuon18_PsiPrime = Var("userInt('HLT_Dimuon18_PsiPrime')",int,doc="reco muon fired this trigger"),
        fired_HLT_Dimuon10_PsiPrime_Barrel_Seagulls = Var("userInt('HLT_Dimuon10_PsiPrime_Barrel_Seagulls')",int,doc="reco muon fired this trigger"),
        fired_HLT_DoubleMu4_PsiPrimeTrk_Displaced = Var("userInt('HLT_DoubleMu4_PsiPrimeTrk_Displaced')",int,doc="reco muon fired this trigger"),
        fired_HLT_Dimuon0_Jpsi3p5_Muon2 = Var("userInt('HLT_Dimuon0_Jpsi3p5_Muon2')",int,doc="reco muon fired this trigger"),
        fired_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi = Var("userInt('HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi')",int,doc="reco muon fired this trigger"),
        fired_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi = Var("userInt('HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi')",int,doc="reco muon fired this trigger"),
        fired_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05 = Var("userInt('HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05')",int,doc="reco muon fired this trigger"),
    ),
)

muonsXMCMatchForTable = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = muonXTable.src,                         # final reco collection
    matched     = cms.InputTag("finalGenParticlesX"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),                             # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.03),                           # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
)

muonXMCTable = cms.EDProducer("CandMCMatchTableProducerX",
    src     = muonXTable.src,
    mcMap   = cms.InputTag("muonsXMCMatchForTable"),
    objName = muonXTable.name,
    objType = muonXTable.name, 
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 muons"),
)

selectedMuonsMCMatchEmbedded = cms.EDProducer(
    'MuonMatchEmbedder',
    src = cms.InputTag('muonTrgSelector:SelectedMuons'),
    matching = cms.InputTag('muonsXMCMatchForTable')
)

muonTriggerMatchedTable = muonXTable.clone(
    src = cms.InputTag("muonTrgSelector:trgMuons"),
    name = cms.string("TriggerMuon"),
    doc  = cms.string("HLT Muons matched with reco muons"), #reco muon matched to triggering muon"),
    variables = cms.PSet(CandVars,
        vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
        vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
        vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6)
   )
)


#muonXSequence = cms.Sequence(muonTrgSelector * countTrgMuons) # chiara
muonXSequence = cms.Sequence(muonTrgSelector)
muonXMC = cms.Sequence(muonXSequence + muonsXMCMatchForTable + selectedMuonsMCMatchEmbedded + muonXMCTable)
muonXTables = cms.Sequence(muonXTable)
muonTriggerMatchedTables = cms.Sequence(muonTriggerMatchedTable)   
