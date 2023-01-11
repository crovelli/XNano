import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

Path2016=["HLT_Dimuon16_Jpsi","HLT_Dimuon13_PsiPrime","HLT_DoubleMu4_JpsiTrk_Displaced","HLT_DoubleMu4_PsiPrimeTrk_Displaced"]

Path2017=["HLT_Dimuon25_Jpsi","HLT_DoubleMu4_JpsiTrk_Displaced"]

Path=Path2017

muonTrgSelector = cms.EDProducer("MuonTriggerSelector",
                                 muonCollection = cms.InputTag("slimmedMuons"), 
                                 beamSpot   = cms.InputTag("offlineBeamSpot"),
                                 bits = cms.InputTag("TriggerResults","","HLT"),
                                 objects = cms.InputTag("slimmedPatTrigger"),
                                 
                                 ## trigger match
                                 drForTriggerMatch = cms.double(0.1),        # to be tuned

                                 ## for the output selected collection 
                                 ptMin = cms.double(0.5),                            
                                 absEtaMax = cms.double(2.4),

                                 HLTPaths=cms.vstring(Path)
                             )


# we need at least 2 triggering muons
countTrgMuons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(2),
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
        variables = cms.PSet( 
        isGlobal = Var("userInt('isGlobal')",bool,doc="muon is global muon"),
        softId = Var("passed('SoftCutBasedId')",bool,doc="soft cut-based ID"), 
        looseId = Var("userInt('looseId')",bool,doc="loose cut-based ID"),
        charge = Var("userInt('charge')",bool,doc="charge"),
        trackQuality = Var("userInt('trackQuality')",int,doc="trackQuality"),
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


muonXSequence = cms.Sequence(muonTrgSelector * countTrgMuons)
muonXMC = cms.Sequence(muonXSequence + muonsXMCMatchForTable + selectedMuonsMCMatchEmbedded + muonXMCTable)
muonXTables = cms.Sequence(muonXTable)
