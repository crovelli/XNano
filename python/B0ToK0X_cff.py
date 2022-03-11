import FWCore.ParameterSet.Config as cms
from PhysicsTools.XNano.common_cff import *


########## inputs preparation ################

# J/Psi -> mumu
muonPairsForJPsiMuMu = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    lep1Selection = cms.string('pt > 2.5'),
    lep2Selection = cms.string('pt > 2.5'),
    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
                                 '&& mass() > 1 && charge() == 0 && userFloat("lep_deltaR") > 0.02'),
    postVtxSelection =  cms.string('userFloat("sv_prob") > 1.e-3'
                                 '&& userFloat("fitted_mass") > 2.50 && userFloat("fitted_mass") < 4.00'  
                               ),
)


# K0s -> pi pi
K0sToPiPi = cms.EDProducer(
    'K0sBuilder',
    pfcands= cms.InputTag('tracksX', 'SelectedTracks'),
    transientTracks= cms.InputTag('tracksX', 'SelectedTransientTracks'),
    trk1Selection = cms.string('pt > 0.4 && abs(eta)<3.0'), 
    trk2Selection = cms.string('pt > 0.4 && abs(eta)<3.0'), 
    preVtxSelection = cms.string('abs(userCand("trk1").vz - userCand("trk2").vz)<=1.0' 
                                 ' && userFloat("trk_deltaR") > 0.02'
                                 ' &&  pt()>0.5 && (mass() < 0.750 && mass() > 0.250)'
                             ),
    postVtxSelection = cms.string('userFloat("sv_prob") > 1.e-2'
                                  ' && (  (userFloat("fitted_mass")<0.700 && userFloat("fitted_mass")>0.300))'
                              )
)


# pi+pi- candidate 
candPiPi = cms.EDProducer(
    'PiPiBuilder',
    pfcands= cms.InputTag('tracksX', 'SelectedTracks'),
    transientTracks= cms.InputTag('tracksX', 'SelectedTransientTracks'),
    trk1Selection = cms.string('pt > 0.4 && abs(eta)<3.0'),
    trk2Selection = cms.string('pt > 0.4 && abs(eta)<3.0'), 
    preVtxSelection = cms.string('abs(userCand("trk1").vz - userCand("trk2").vz) <= 1.0'
                                 '&& userFloat("trk_deltaR") > 0.02'),
    postVtxSelection = cms.string('userFloat("sv_prob") > 1.e-3'
                                  ' && userFloat("fitted_mass")>0.400'
                                  )
)


########################### B0 -> pi pi pi pi ll ##########################
BToK0sMuMuPiPi = cms.EDProducer(
    'BToK0sMuMuPiPiBuilder',

    dimuons = cms.InputTag('muonPairsForJPsiMuMu', 'SelectedDiMuons'),
    dimuonKinVtxs = cms.InputTag('muonPairsForJPsiMuMu', 'SelectedDiMuonKinVtxs'),
    muonTransientTracks = muonPairsForJPsiMuMu.transientTracksSrc,

    k0short = cms.InputTag('K0sToPiPi', 'SelectedK0s'),
    k0shortKinVtxs = cms.InputTag('K0sToPiPi', 'SelectedK0sKinVtxs'),  
    k0shortTransientTracks = cms.InputTag('tracksX', 'SelectedTransientTracks'),

    dipion = cms.InputTag('candPiPi', 'SelectedDiPions'),
    dipionKinVtxs = cms.InputTag('candPiPi', 'SelectedDiPionKinVtxs'), 
    pionTransientTracks = cms.InputTag('tracksX', 'SelectedTransientTracks'),

    beamSpot = cms.InputTag("offlineBeamSpot"),
    offlinePrimaryVertexSrc = cms.InputTag('offlineSlimmedPrimaryVertices'), 

    # to compute isolation variables
    #tracks = cms.InputTag("packedPFCandidates"),
    #lostTracks = cms.InputTag("lostTracks"),
    #isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.5'),
    
    preVtxSelection = cms.string(
        'pt > 0.5 && userFloat("min_dr") > 0.03'
        '&& (mass < 7. && mass > 4.) '
        ),

    postVtxSelection = cms.string(
        'userFloat("sv_prob") > 0.01 '
        '&& userFloat("sv_chi2") > 0. '
        #'&& userFloat("fitted_cosTheta2D_PV") >= 0'
        #'&& userFloat("lxySign_PV") >= 3'
        '&& (userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.)'
        '&& (userFloat("fitted_X_mass") > 3.0 && userFloat("fitted_X_mass") < 4.5)'
    )
)


################################### Tables #####################################

JPsiToMuMuTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag('muonPairsForJPsiMuMu', 'SelectedDiMuons'),
    cut = cms.string(""),
    name = cms.string("JPsiToMuMu"),
    doc = cms.string("JPsi Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
      CandVars,
      fitted_mass = ufloat('fitted_mass'),
      fitted_massErr = ufloat('fitted_massErr'),
      svprob = ufloat('sv_prob'),        
      svndof = ufloat('sv_ndof'),
      svchi2 = ufloat('sv_chi2'),         
      l1_idx = uint('l1_idx'),
      l2_idx = uint('l2_idx') #,
      #fitted_l1_pt = ufloat('fitted_l1_pt'),  
      #fitted_l2_pt = ufloat('fitted_l2_pt'),  
      #fitted_l1_eta = ufloat('fitted_l1_eta'),  
      #fitted_l2_eta = ufloat('fitted_l2_eta'),  
      #fitted_l1_phi = ufloat('fitted_l1_phi'),  
      #fitted_l2_phi = ufloat('fitted_l2_phi')
    )
)

K0sToPiPiTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("K0sToPiPi","SelectedK0s"),
    cut = cms.string(""),
    name = cms.string("K0s"),
    doc = cms.string("K0s Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
      CandVars,
      fitted_mass = ufloat('fitted_mass'),
      fitted_pt = ufloat('fitted_pt'),
      fitted_eta = ufloat('fitted_eta'),
      fitted_phi = ufloat('fitted_phi'),
      svprob = ufloat('sv_prob'),         
      svndof = ufloat('sv_ndof'),
      svchi2 = ufloat('sv_chi2'),         
      trk1_idx = uint('trk1_idx'),
      trk2_idx = uint('trk2_idx')  #,
      #fitted_trk1_pt = ufloat('fitted_trk1_pt'),  
      #fitted_trk2_pt = ufloat('fitted_trk2_pt'),  
      #fitted_trk1_eta = ufloat('fitted_trk1_eta'),  
      #fitted_trk2_eta = ufloat('fitted_trk2_eta'),  
      #fitted_trk1_phi = ufloat('fitted_trk1_phi'),  
      #fitted_trk2_phi = ufloat('fitted_trk2_phi')
    )
)

candPiPiTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("candPiPi","SelectedDiPions"),
    cut = cms.string(""),
    name = cms.string("pipi"),
    doc = cms.string("pi-pi Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
      CandVars,
      fitted_mass = ufloat('fitted_mass'),
      svprob = ufloat('sv_prob'),         
      svndof = ufloat('sv_ndof'),
      svchi2 = ufloat('sv_chi2'),         
      trk1_idx = uint('trk1_idx'),
      trk2_idx = uint('trk2_idx')    #,
      #fitted_trk1_pt = ufloat('fitted_trk1_pt'),  
      #fitted_trk2_pt = ufloat('fitted_trk2_pt'),  
      #fitted_trk1_eta = ufloat('fitted_trk1_eta'),  
      #fitted_trk2_eta = ufloat('fitted_trk2_eta'),  
      #fitted_trk1_phi = ufloat('fitted_trk1_phi'),  
      #fitted_trk2_phi = ufloat('fitted_trk2_phi')
    )
)

BToK0sMuMuPiPiTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToK0sMuMuPiPi"),
    cut = cms.string(""),
    name = cms.string("B0"),
    doc = cms.string("BToK0sMuMuPiPi Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        mu1_idx = uint('mu1_idx'),
        mu2_idx = uint('mu2_idx'),
        trk1Rho_idx = uint('trk1Rho_idx'),
        trk2Rho_idx = uint('trk2Rho_idx'),
        trk1k0s_idx = uint('trk1k0s_idx'),
        trk2k0s_idx = uint('trk2k0s_idx'),
        k0short_idx = uint('k0short_idx'),
        dipion_idx  = uint('dipion_idx'),
        dilepton_idx = uint('dilepton_idx'),
        pv_idx = uint('pv_idx'),
        min_dr = ufloat('min_dr'), 
        max_dr = ufloat('max_dr'),   
        # fit and vtx info
        svchi2 = ufloat('sv_chi2'),
        svndof = ufloat('sv_ndof'),
        svprob = ufloat('sv_prob'),
        lxyBS = ufloat('lxy_BS'),
        lxyPV = ufloat('lxy_PV'),
        lxyUncBS = ufloat('lxyUnc_BS'),
        lxyUncPV = ufloat('lxyUnc_PV'),
        # k0short/Rho/JPsi fitted in B0 vertex
        fit_k0short_mass = ufloat('fitted_k0short_mass'),
        fitted_Rho_mass  = ufloat('fitted_Rho_mass'),
        fitted_JPsi_mass = ufloat('fitted_JPsi_mass'),
        fitted_X_mass = ufloat('fitted_X_mass'),
        # Cos(theta)
        cos2D_BS = ufloat('cosTheta2D_BS'),
        cos2D_PV = ufloat('cosTheta2D_PV'),
        fit_cos2D_BS = ufloat('fitted_cosTheta2D_BS'),
        fit_cos2D_PV = ufloat('fitted_cosTheta2D_PV'),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),
        fit_pt = ufloat('fitted_pt'),
        fit_eta = ufloat('fitted_eta'),
        fit_phi = ufloat('fitted_phi'),
        # post-fit tracks/leptons
        #l1
        fit_mu1_pt  = ufloat('fitted_mu1_pt'),
        fit_mu1_eta = ufloat('fitted_mu1_eta'),
        fit_mu1_phi = ufloat('fitted_mu1_phi'),
        #l2
        fit_mu2_pt  = ufloat('fitted_mu2_pt'),
        fit_mu2_eta = ufloat('fitted_mu2_eta'),
        fit_mu2_phi = ufloat('fitted_mu2_phi'),
        #trk1
        fit_trk1Rho_pt  = ufloat('fitted_trk1Rho_pt'),
        fit_trk1Rho_eta = ufloat('fitted_trk1Rho_eta'),
        fit_trk1Rho_phi = ufloat('fitted_trk1Rho_phi'),
        #trk2
        fit_trk2Rho_pt  = ufloat('fitted_trk2Rho_pt'),
        fit_trk2Rho_eta = ufloat('fitted_trk2Rho_eta'),
        fit_trk2Rho_phi = ufloat('fitted_trk2Rho_phi'),
        #trk1k0s
        fit_trk1k0s_pt  = ufloat('fitted_trk1k0s_pt'),
        fit_trk1k0s_eta = ufloat('fitted_trk1k0s_eta'),
        fit_trk1k0s_phi = ufloat('fitted_trk1k0s_phi'),
        #trk2k0s
        fit_trk2k0s_pt  = ufloat('fitted_trk2k0s_pt'),
        fit_trk2k0s_eta = ufloat('fitted_trk2k0s_eta'),
        fit_trk2k0s_phi = ufloat('fitted_trk2k0s_phi'),
        #
        fit2_mass = ufloat('fitted2_mass'),
        # isolation 
        #l1_iso03 = ufloat('l1_iso03'),
        #l1_iso04 = ufloat('l1_iso04'),
        #l2_iso03 = ufloat('l2_iso03'),
        #l2_iso04 = ufloat('l2_iso04'),
        #tk1_iso03 = ufloat('tk1_iso03'),
        #tk1_iso04 = ufloat('tk1_iso04'),
        #tk2_iso03 = ufloat('tk2_iso03'),
        #tk2_iso04 = ufloat('tk2_iso04'),
        #b_iso03  = ufloat('b_iso03'),
        #b_iso04  = ufloat('b_iso04'),
    )
)

CountBToK0sMuMuPiPi = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToK0sMuMuPiPi")
)    

CountMuonPairs = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag('muonPairsForJPsiMuMu', 'SelectedDiMuons'),
)    

CountPiPi = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("candPiPi", "SelectedDiPions"),
)    

CountK0s = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("K0sToPiPi","SelectedK0s")
)    


########################### Sequencies  ############################

JPsiMuMuSequence = cms.Sequence(
    (muonPairsForJPsiMuMu * CountMuonPairs)
)

JPsiToMuMuTableSequence = cms.Sequence( JPsiToMuMuTable )

K0sToPiPiSequence = cms.Sequence(  
    ( K0sToPiPi * CountK0s )
)

K0sToPiPiTableSequence = cms.Sequence( K0sToPiPiTable )   

candPiPiSequence = cms.Sequence(
    ( candPiPi * CountPiPi )
)

candPiPiTableSequence = cms.Sequence( candPiPiTable )

BToK0sMuMuPiPiSequence = cms.Sequence( 
    ( BToK0sMuMuPiPi * CountBToK0sMuMuPiPi )
)

BToK0sMuMuPiPiTableSequence = cms.Sequence( BToK0sMuMuPiPiTable )

