from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

options.register('isMC', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('Era','2022postEE',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set data taking period : 2022preEE, 2022postEE, 2023preBPix, 2023postBPix"
)
options.register('physProcess','X3872',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set the phys process : X3872"
)
options.register('globalTag', 'NOTSET',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)
options.register('wantSummary', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('wantFullRECO', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('reportEvery', 10,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)
options.register('skip', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip first N events"
)

options.parseArguments()

era = options.Era
phys_process = options.physProcess
if not options.tag :
    tag = '_'.join([phys_process, era]) 
    options.setDefault('tag', tag)

# set number of events
options.setDefault('maxEvents', 100)

# set global tag:
global_tags_mc = {
    '2022preEE'     : '130X_mcRun3_2022_realistic_v5',
    '2022postEE'    : '130X_mcRun3_2022_realistic_postEE_v6',
    '2023preBPix'   : '130X_mcRun3_2023_realistic_v14',
    '2023postBPix'  : '130X_mcRun3_2023_realistic_postBPix_v2',
}
global_tags_data = {
    '2022preEE'     : '124X_dataRun3_PromptAnalysis_v1', #era CD use '124X_dataRun3_Prompt_v10' for era E
    '2022postEE'    : '130X_dataRun3_PromptAnalysis_v1', #era FG
    '2023preBPix'   : '130X_dataRun3_PromptAnalysis_v1', #era BC
    '2023postBPix'  : '130X_dataRun3_PromptAnalysis_v1', #era D
}
 
if options._beenSet['globalTag']:
    globaltag = options.globalTag
else:
    globaltag = global_tags_mc[era] if options.isMC else global_tags_data[era] 

extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['xNANO', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['xFullEvt', extension[options.isMC], options.tag])+'.root')

# test input files per process and era
test_mc_inFiles_process_era = {
     'X3872' : {
         '2022preEE'     : ['/store/mc/Run3Summer22MiniAODv4/BuToX3872Ks_X3872ToJPsiRho_JPsiToMuMu_RhoToPiPi_TuneCP5_13TeV_pythia8-evtgen/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/3db65477-4aea-4dd6-a194-fa111395ba8e.root'],   # NB: si chiama 13 TeV ma sembra essere 13.6, https://cms-pdmv-prod.web.cern.ch/mcm/requests?prepid=BPH-Run3Summer22GS-00020&shown=262271
        '2022postEE'    : ['/store/mc/Run3Summer22EEMiniAODv4/BuToX3872Ks_X3872ToJPsiRho_JPsiToMuMu_RhoToPiPi_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/60000/02bdd816-3a70-4b3d-bdf9-021756af1a04.root'],
        '2023preBPix'   : [''],
        '2023postBPix'  : [''], 
    },
}

test_data_inFiles_era = {
    '2022preEE'     : [], #era CDE
    '2022postEE'    : ['/store/data/Run2022F/ParkingDoubleMuonLowMass0/MINIAOD/22Sep2023-v1/60000/00807cd1-61e0-4a30-917a-f2957e4365b0.root',
                        '/store/data/Run2022F/ParkingDoubleMuonLowMass0/MINIAOD/22Sep2023-v1/60000/0163f2af-82a7-4540-8702-b111818a93d4.root'
                    ], #era FG
    '2023preBPix'   : [], #era BC
    '2023postBPix'  : [], #era D
}

# Define inputs
if not options.inputFiles :
    options.inputFiles = test_mc_inFiles_process_era[phys_process][era] if options.isMC else test_data_inFiles_era[era]
annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

#from Configuration.StandardSequences.Eras import eras
process = cms.Process('XNANO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.XNano.nanoX_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
)

process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileNANO,
    outputCommands = cms.untracked.vstring(
      'drop *',
      "keep nanoaodFlatTable_*Table_*_*",     # event data
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

from PhysicsTools.XNano.nanoX_cff import *
process = nanoAOD_customizeMuonTriggerX(process)
process = nanoAOD_customizeTrackFilteredX(process)
process = nanoAOD_customizeB0ToK0X(process)
process = nanoAOD_customizeTriggerBitsX(process)

# Path and EndPath definitions
process.nanoAOD_B0ToK0X_step = cms.Path(process.nanoSequence + process.nanoB0ToK0XSequence )

# customisation of the process.
if options.isMC:
    from PhysicsTools.XNano.nanoX_cff import nanoAOD_customizeMC
    nanoAOD_customizeMC(process)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
                                process.nanoAOD_B0ToK0X_step,
                                process.endjob_step,
                                process.NANOAODoutput_step
                               )

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
                                   'nanoAOD_B0ToK0X_step'
                                   )
)


process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
