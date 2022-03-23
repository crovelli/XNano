from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.XNano.trgbits_cff import *

## for gen and trigger muon
from PhysicsTools.XNano.genparticlesX_cff import *
from PhysicsTools.XNano.particlelevelX_cff import *
from PhysicsTools.XNano.triggerObjectsX_cff import *
from PhysicsTools.XNano.muonsX_cff import * 

## filtered input collections
from PhysicsTools.XNano.tracksX_cff import *

## B collections
from PhysicsTools.XNano.B0ToK0X_cff import *


nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectXTables + l1bits)

nanoSequence = cms.Sequence(nanoMetadata + 
                            vertexSequence +           
                            globalTables + vertexTables + 
                            triggerObjectXTables + l1bits)

nanoSequenceMC = cms.Sequence(particleLevelXSequence + genParticleXSequence + 
                              globalTablesMC + genWeightsTable + genParticleXTables + lheInfoTable) 


def nanoAOD_customizeMuonTriggerX(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + muonXSequence + muonXTables)
    return process

def nanoAOD_customizeTrackFilteredX(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + tracksXSequence + tracksXTables)
    return process

def nanoAOD_customizeTriggerBitsX(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + trgTables)
    return process

def nanoAOD_customizeB0ToK0X(process):
    process.nanoB0ToK0XSequence = cms.Sequence( JPsiMuMuSequence + JPsiToMuMuTableSequence + candPiPiSequence + candPiPiTableSequence + K0sToPiPiSequence + K0sToPiPiTableSequence + BToK0sMuMuPiPiSequence + BToK0sMuMuPiPiTableSequence )
    return process


from FWCore.ParameterSet.MassReplace import massSearchReplaceAnyInputTag
def nanoAOD_customizeMC(process):
    for name, path in process.paths.iteritems():
        # replace all the non-match embedded inputs with the matched ones
        massSearchReplaceAnyInputTag(path, 'muonTrgSelector:SelectedMuons', 'selectedMuonsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'tracksX:SelectedTracks', 'tracksXMCMatchEmbedded')

        # modify the path to include mc-specific info
        path.insert(0, nanoSequenceMC)
        path.replace(process.muonXSequence, process.muonXMC)
        path.replace(process.tracksXSequence, process.tracksXMC)
