import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

tracksX = cms.EDProducer('TrackMerger',
                         beamSpot   = cms.InputTag("offlineBeamSpot"),
                         trgMuon    = cms.InputTag("muonTrgSelector:trgMuons"),
                         tracks     = cms.InputTag("packedPFCandidates"),
                         lostTracks = cms.InputTag("lostTracks"),
                         trkPtCut = cms.double(0.4),
                         trkEtaCut = cms.double(3.0),
                         muons      = cms.InputTag("slimmedMuons"),
                         vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                         dzTrg_cleaning = cms.double(-1.),    # was 1
                         drTrg_Cleaning = cms.double(-1.),    # was 0.03
                         dcaSig = cms.double(-1.),
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
        CandVars,
        vx = Var("vx()", float, doc="x coordinate of vertex position, in cm", precision=10),
        vy = Var("vy()", float, doc="y coordinate of vertex position, in cm", precision=10),
        vz = Var("vz()", float, doc="z coordinate of vertex position, in cm", precision=10),
        isPacked = Var("userInt('isPacked')",int,doc="track from packedCandidate collection", precision=10),
        isLostTrk = Var("userInt('isLostTrk')",int,doc="track from lostTrack collection", precision=10),
        dz = Var("userFloat('dz')",float,doc="dz (with sign) wrt first PV, in cm", precision=10),
        dxy = Var("userFloat('dxy')",float,doc="dxy (with sign) wrt first PV, in cm", precision=10),
        dzS = Var("userFloat('dzS')", float, doc="dz/err (with sign) wrt first PV, in cm", precision=10),
        dxyS = Var("userFloat('dxyS')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),
        DCASig=Var("userFloat('DCASig')", float,doc="significance of xy-distance of closest approach wrt beamspot", precision=10),
        dzTrg = Var("userFloat('dzTrg')", float,doc="dz from the corresponding trigger muon, in cm", precision=10),
        isMatchedToMuon = Var("userInt('isMatchedToMuon')",bool,doc="track was used to build a muon", precision=10),
        isMatchedToLooseMuon = Var("userInt('isMatchedToLooseMuon')",bool,doc="track was used to build a muon passing LooseID", precision=10),
        isMatchedToSoftMuon = Var("userInt('isMatchedToSoftMuon')",bool,doc="track was used to build a muon passing softID", precision=10),
        isMatchedToMediumMuon = Var("userInt('isMatchedToMediumMuon')",bool,doc="track was used to build a muon passing mediumID", precision=10),
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


