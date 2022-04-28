import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

Path2016=["HLT_DoubleMu4_JpsiTrk_Displaced","HLT_DoubleMu4_PsiPrimeTrk_Displaced"]

Path2017=["HLT_DoubleMu4_JpsiTrk_Displaced","HLT_DoubleMu4_PsiPrimeTrk_Displaced","HLT_DoubleMu4_JpsiTrkTrk_Displaced"]

Path=Path2017

tracksX = cms.EDProducer('TrackMerger',
                         beamSpot   = cms.InputTag("offlineBeamSpot"),
                         tracks     = cms.InputTag("packedPFCandidates"),
                         lostTracks = cms.InputTag("lostTracks"),
                         ## trigger match  
                         bits = cms.InputTag("TriggerResults","","HLT"), 
                         objects = cms.InputTag("slimmedPatTrigger"), 
                         drForTriggerMatch = cms.double(0.1),           # to be tuned
                         HLTPaths=cms.vstring(Path),
                         #
                         trkPtCut = cms.double(0.5),
                         trkEtaCut = cms.double(3.0),
                         muons      = cms.InputTag("slimmedMuons"),
                         vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                         trkNormChiMin = cms.int32(-1),
                         trkNormChiMax = cms.int32(-1)
                        )


trackXTable = cms.EDProducer(
    "SimpleCompositeCandidateFlatTableProducer",
    src = cms.InputTag("tracksX:SelectedTracks"),
    cut = cms.string(""),
    name = cms.string("ProbeTracks"),
    doc  = cms.string("track collection after basic selection"),
    singleton = cms.bool(False),
    extension = cms.bool(False), 
    variables = cms.PSet(
        #CandVars,
        #vx = Var("vx()", float, doc="x coordinate of vertex position, in cm", precision=10),
        #vy = Var("vy()", float, doc="y coordinate of vertex position, in cm", precision=10),
        #vz = Var("vz()", float, doc="z coordinate of vertex position, in cm", precision=10),
        isMatchedToMuon = Var("userInt('isMatchedToMuon')",bool,doc="track was used to build a muon", precision=10),
        isMatchedToLooseMuon = Var("userInt('isMatchedToLooseMuon')",bool,doc="track was used to build a muon passing LooseID", precision=10),
        isMatchedToSoftMuon = Var("userInt('isMatchedToSoftMuon')",bool,doc="track was used to build a muon passing softID", precision=10),
        nValidHits = Var("userInt('nValidHits')", int,doc="Number of valid hits on track", precision=10),
    ),
)

tracksXMCMatchForTable = cms.EDProducer("MCMatcher",   # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = trackXTable.src,                     # final reco collection
    matched     = cms.InputTag("finalGenParticlesX"),  # final mc-truth particle collection
    mcPdgId     = cms.vint32(321,211),                 # one or more PDG ID (321 = charged kaon, 211 = charged pion); absolute values (see below)
    checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.03),             # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
)

tracksXMCMatchEmbedded = cms.EDProducer(
    'CompositeCandidateMatchEmbedder',
    src = trackXTable.src,
    matching = cms.InputTag("tracksXMCMatchForTable")
)

tracksXMCTable = cms.EDProducer("CandMCMatchTableProducerX",
    src     = tracksXMCMatchForTable.src,
    mcMap   = cms.InputTag("tracksXMCMatchForTable"),
    objName = trackXTable.name,
    objType = trackXTable.name,
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 kaons or pions"),
)


tracksXSequence = cms.Sequence(tracksX)
tracksXTables = cms.Sequence(trackXTable)
tracksXMC = cms.Sequence(tracksXSequence + tracksXMCMatchForTable + tracksXMCMatchEmbedded + tracksXMCTable)


