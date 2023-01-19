import FWCore.ParameterSet.Config as cms
from PhysicsTools.XNano.common_cff import *


########## inputs preparation ################

# J/Psi -> mumu
muonPairsForJPsiMuMu = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    beamSpot   = cms.InputTag("offlineBeamSpot"),
    lep1Selection = cms.string('pt > 2.5'),
    lep2Selection = cms.string('pt > 2.5'),
    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
                                 '&& mass() > 1 && charge() == 0 && userFloat("lep_deltaR") > 0.02'),
    postVtxSelection =  cms.string('userFloat("sv_prob") > 0.01'
                                 '&& userFloat("fitted_mass") > 2.50 && userFloat("fitted_mass") < 4.00'  
                               ),
)

# K0s -> pi pi
K0sToPiPi = cms.EDProducer(
    'K0sBuilder',
    dipion = cms.InputTag('candPiPi', 'SelectedDiPions'),
    dimuons = cms.InputTag('muonPairsForJPsiMuMu', 'SelectedDiMuons'),
    svSrc = cms.InputTag("slimmedKshortVertices"),
    postVtxSelection = cms.string('userFloat("sv_prob") > 1.e-2'
                                  ' && (  (userFloat("fitted_mass_womc")<0.700 && userFloat("fitted_mass_womc")>0.300))'
                              )
)

# pi+pi- candidate 
candPiPi = cms.EDProducer(
    'PiPiBuilder',
    transientTracksSrc = cms.InputTag('tracksX', 'SelectedTransientTracks'),
    dimuons = cms.InputTag('muonPairsForJPsiMuMu', 'SelectedDiMuons'),
    pfcands= cms.InputTag('tracksX', 'SelectedTracks'),
    trkSelection = cms.string('pt > 0.4 && abs(eta)<3.0'),
    preVtxSelection = cms.string('abs(userCand("trk1").vz - userCand("trk2").vz) <= 1.0'
                                 '&& userFloat("trk_deltaR") > 0.02')
)


########################### B0 -> pi pi pi pi ll ##########################
BToK0sMuMuPiPi = cms.EDProducer(
    'BToK0sMuMuPiPiBuilder',

    pfcands= cms.InputTag('tracksX', 'SelectedTracks'),

    dimuons = cms.InputTag('muonPairsForJPsiMuMu', 'SelectedDiMuons'),
    muonTransientTracks = muonPairsForJPsiMuMu.transientTracksSrc,

    k0short = cms.InputTag('K0sToPiPi', 'SelectedK0s'),
    k0shortKinVtxsWMC = cms.InputTag('K0sToPiPi', 'SelectedK0sKinVtxsWMC'),  

    dipion = cms.InputTag('candPiPi', 'SelectedDiPions'),
    pionTransientTracks = cms.InputTag('tracksX', 'SelectedTransientTracks'),

    beamSpot = cms.InputTag("offlineBeamSpot"),
    offlinePrimaryVertexSrc = cms.InputTag('offlineSlimmedPrimaryVertices'), 

    postVtxSelection = cms.string(
        'userFloat("sv_prob") > 0.01 '
        '&& userFloat("sv_chi2") > 0. '
        '&& (userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.)'
        '&& (userFloat("finalFit_X_mass") > 3.0 && userFloat("finalFit_X_mass") < 4.5)'
    ),

    postVtxSelection2 = cms.string(
#        'userFloat("fit_cos3D_PV") >= 0'
#        '&& userFloat("lxySign_PV") >= 3'
        'userFloat("lxySign_PV") >= 3'
    ),
    
    drMatchTrack = cms.double(0.05)
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
        # CandVars,  
        # index
        mu1_idx = uint('mu1_idx'),
        mu2_idx = uint('mu2_idx'),
        pi1_idx = uint('trk1Rho_idx'),
        pi2_idx = uint('trk2Rho_idx'),
        k0short_idx = uint('k0short_idx'),
        dipion_idx  = uint('dipion_idx'),
        dimuon_idx = uint('dimuon_idx'),
        # fitted p4
        fitted_mass_womc = ufloat('fitted_mass_womc'),
        finalFit_mass = ufloat('fitted_mass'),
        finalFit_pt = ufloat('fitted_pt'),
        finalFit_eta = ufloat('fitted_eta'),
        finalFit_phi = ufloat('fitted_phi'),
        # fitted daughters
        finalFit_X_mass = ufloat('finalFit_X_mass'),        
        finalFit_Rho_mass  = ufloat('finalFit_Rho_mass'),
        finalFit_JPsi_mass = ufloat('finalFit_JPsi_mass'),
        finalFit_mu1_pt   = ufloat('finalFit_mu1_pt'),
        finalFit_mu1_eta  = ufloat('finalFit_mu1_eta'),
        finalFit_mu1_phi  = ufloat('finalFit_mu1_phi'),
        finalFit_mu2_pt   = ufloat('finalFit_mu2_pt'),
        finalFit_mu2_eta  = ufloat('finalFit_mu2_eta'),
        finalFit_mu2_phi  = ufloat('finalFit_mu2_phi'),
        finalFit_pi1_pt  = ufloat('finalFit_pi1Rho_pt'),
        finalFit_pi1_eta = ufloat('finalFit_pi1Rho_eta'),
        finalFit_pi1_phi = ufloat('finalFit_pi1Rho_phi'),
        finalFit_pi2_pt  = ufloat('finalFit_pi2Rho_pt'),
        finalFit_pi2_eta = ufloat('finalFit_pi2Rho_eta'),
        finalFit_pi2_phi = ufloat('finalFit_pi2Rho_phi'),
        finalFit_k0s_pt   = ufloat('finalFit_k0s_pt'),
        finalFit_k0s_eta  = ufloat('finalFit_k0s_eta'),
        finalFit_k0s_phi  = ufloat('finalFit_k0s_phi'),        
        # fit and vtx info
        svchi2 = ufloat('sv_chi2'),
        svprob = ufloat('sv_prob'),
        decayVtxX = ufloat('decayVtxX'),
        decayVtxY = ufloat('decayVtxY'),
        decayVtxZ = ufloat('decayVtxZ'),
        decayVtxXE = ufloat('decayVtxXE'),
        decayVtxYE = ufloat('decayVtxYE'),
        decayVtxZE = ufloat('decayVtxZE'),
        # Cos(theta) 
        lxySign_PV = ufloat('lxySign_PV'),
        lxySign_BS = ufloat('lxySign_BS'),
        cosAlpha3D_PV = ufloat('cosAlpha3D_PV'),
        #cosAlpha2D_PV = ufloat('cosAlpha2D_PV'),
        cosAlpha2D_BS = ufloat('cosAlpha2D_BS'),
        # vtx
        pv3D_idx = uint('pv3D_idx'),
        #pv2D_idx = uint('pv2D_idx'),
        PVx  = ufloat('PVx'),
        PVy  = ufloat('PVy'),
        PVz  = ufloat('PVz'),
        PVEx = ufloat('PVEx'),
        PVEy = ufloat('PVEy'),
        PVEz = ufloat('PVEz'),
        # mumu
        MuMu_sv_prob = ufloat('MuMu_sv_prob'),
        MuMu_fitted_mass = ufloat('MuMu_fitted_mass'),
        MuMu_fitted_pt  = ufloat('MuMu_fitted_pt'),
        MuMu_fitted_eta = ufloat('MuMu_fitted_eta'),
        MuMu_fitted_phi = ufloat('MuMu_fitted_phi'),
        MuMu_fitted_vtxX = ufloat('MuMu_fitted_vtxX'),
        MuMu_fitted_vtxY = ufloat('MuMu_fitted_vtxY'),
        MuMu_fitted_vtxZ = ufloat('MuMu_fitted_vtxZ'),
        MuMu_fitted_vtxXE = ufloat('MuMu_fitted_vtxXE'),
        MuMu_fitted_vtxYE = ufloat('MuMu_fitted_vtxYE'),
        MuMu_fitted_vtxZE = ufloat('MuMu_fitted_vtxZE'),
        MuMu_prefit_mu1_pt  = ufloat('MuMu_prefit_mu1_pt'),
        MuMu_prefit_mu1_eta = ufloat('MuMu_prefit_mu1_eta'),
        MuMu_prefit_mu1_phi = ufloat('MuMu_prefit_mu1_phi'),
        MuMu_prefit_mu2_pt  = ufloat('MuMu_prefit_mu2_pt'),
        MuMu_prefit_mu2_eta = ufloat('MuMu_prefit_mu2_eta'),
        MuMu_prefit_mu2_phi = ufloat('MuMu_prefit_mu2_phi'),
        MuMu_mu1_dr = ufloat('MuMu_mu1_dr'),
        MuMu_mu2_dr = ufloat('MuMu_mu2_dr'),
        MuMu_DCA = ufloat('MuMu_DCA'),
        MuMu_LxySign = ufloat('MuMu_LxySign'),
        MuMu_cosAlpha = ufloat('MuMu_cosAlpha'),
        MuMu_mu1_dzsign  = ufloat('MuMu_mu1_dzsign'),
        MuMu_mu2_dzsign  = ufloat('MuMu_mu2_dzsign'),
        MuMu_mu1_dxysign = ufloat('MuMu_mu1_dxysign'),
        MuMu_mu2_dxysign = ufloat('MuMu_mu2_dxysign'),
        MuMu_mu1_trackQuality = uint('MuMu_mu1_trackQuality'),
        MuMu_mu2_trackQuality = uint('MuMu_mu2_trackQuality'),
        MuMu_mu1_fired_Dimuon25_Jpsi                = uint('MuMu_mu1_fired_Dimuon25_Jpsi'),
        MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced  = uint('MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced'),
        MuMu_mu2_fired_Dimuon25_Jpsi                = uint('MuMu_mu2_fired_Dimuon25_Jpsi'),
        MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced  = uint('MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced'),
        MuMu_mu1_dr_Dimuon25_Jpsi                   = ufloat('MuMu_mu1_dr_Dimuon25_Jpsi'),
        MuMu_mu1_dr_DoubleMu4_JpsiTrk_Displaced     = ufloat('MuMu_mu1_dr_DoubleMu4_JpsiTrk_Displaced'),
        MuMu_mu2_dr_Dimuon25_Jpsi                   = ufloat('MuMu_mu2_dr_Dimuon25_Jpsi'),
        MuMu_mu2_dr_DoubleMu4_JpsiTrk_Displaced     = ufloat('MuMu_mu2_dr_DoubleMu4_JpsiTrk_Displaced'),

        # pipi   
        PiPi_prefit_pi1_pt  = ufloat('PiPi_prefit_pi1_pt'),
        PiPi_prefit_pi1_eta = ufloat('PiPi_prefit_pi1_eta'),
        PiPi_prefit_pi1_phi = ufloat('PiPi_prefit_pi1_phi'),
        PiPi_prefit_pi2_pt  = ufloat('PiPi_prefit_pi2_pt'),
        PiPi_prefit_pi2_eta = ufloat('PiPi_prefit_pi2_eta'),
        PiPi_prefit_pi2_phi = ufloat('PiPi_prefit_pi2_phi'),
        PiPi_prefit_pi1_vx  = ufloat('PiPi_prefit_pi1_vx'),
        PiPi_prefit_pi1_vy  = ufloat('PiPi_prefit_pi1_vy'),
        PiPi_prefit_pi1_vz  = ufloat('PiPi_prefit_pi1_vz'),
        PiPi_prefit_pi2_vx  = ufloat('PiPi_prefit_pi2_vx'),
        PiPi_prefit_pi2_vy  = ufloat('PiPi_prefit_pi2_vy'),
        PiPi_prefit_pi2_vz  = ufloat('PiPi_prefit_pi2_vz'),
        PiPi_pi1_d0sig      = ufloat('PiPi_pi1_d0sig'),
        PiPi_pi2_d0sig      = ufloat('PiPi_pi2_d0sig'),
        PiPi_pi1_maxd0PV    = ufloat('PiPi_pi1_maxd0PV'),
        PiPi_pi2_maxd0PV    = ufloat('PiPi_pi2_maxd0PV'),
        PiPi_pi1_mind0PV    = ufloat('PiPi_pi1_mind0PV'),
        PiPi_pi2_mind0PV    = ufloat('PiPi_pi2_mind0PV'),
        PiPi_pi1_dzsign     = ufloat('PiPi_pi1_dzsign'),
        PiPi_pi2_dzsign     = ufloat('PiPi_pi2_dzsign'),
        PiPi_pi1_dxysign    = ufloat('PiPi_pi1_dxysign'),
        PiPi_pi2_dxysign    = ufloat('PiPi_pi2_dxysign'),
        PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced = uint('PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced'),
        PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced = uint('PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced'),
        PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced    = ufloat('PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced'),
        PiPi_p2_dr_DoubleMu4_JpsiTrk_Displaced    = ufloat('PiPi_p2_dr_DoubleMu4_JpsiTrk_Displaced'),
        PiPi_sv_prob = ufloat('PiPi_sv_prob'),

        # K0s
        K0s_prefit_mass = ufloat('K0s_prefit_mass'),
        K0s_nmcFitted_mass = ufloat('K0s_nmcFitted_mass'),
        K0s_nmcFitted_pi1pt  = ufloat('K0s_nmcFitted_pi1pt'),
        K0s_nmcFitted_pi1eta = ufloat('K0s_nmcFitted_pi1eta'),
        K0s_nmcFitted_pi1phi = ufloat('K0s_nmcFitted_pi1phi'),
        K0s_nmcFitted_pi2pt  = ufloat('K0s_nmcFitted_pi2pt'),
        K0s_nmcFitted_pi2eta = ufloat('K0s_nmcFitted_pi2eta'),
        K0s_nmcFitted_pi2phi = ufloat('K0s_nmcFitted_pi2phi'),
        K0s_mcFitted_svprob  = ufloat('K0s_mcFitted_svprob'),
        K0s_mcFitted_mass = ufloat('K0s_mcFitted_mass'),
        K0s_mcFitted_pt   = ufloat('K0s_mcFitted_pt'),
        K0s_mcFitted_eta  = ufloat('K0s_mcFitted_eta'),
        K0s_mcFitted_phi  = ufloat('K0s_mcFitted_phi'),
        K0s_mcFitted_vtxX = ufloat('K0s_mcFitted_vtxX'),
        K0s_mcFitted_vtxY = ufloat('K0s_mcFitted_vtxY'),
        K0s_mcFitted_vtxZ = ufloat('K0s_mcFitted_vtxZ'),
        K0s_mcFitted_vtxXE = ufloat('K0s_mcFitted_vtxXE'),
        K0s_mcFitted_vtxYE = ufloat('K0s_mcFitted_vtxYE'),
        K0s_mcFitted_vtxZE = ufloat('K0s_mcFitted_vtxZE'),
        K0_lxySign_wrtBvtx = ufloat('K0_lxySign_wrtBvtx'),
        K0_cosAlpha3D = ufloat('K0_cosAlpha3D'),
        #K0_cosAlpha2D = ufloat('K0_cosAlpha2D'),
        K0s_matchTrack1_D0sign = ufloat('K0s_matchTrack1_D0sign'),
        K0s_matchTrack2_D0sign = ufloat('K0s_matchTrack2_D0sign'),
        K0s_matchTrack1_maxD0Pv = ufloat('K0s_matchTrack1_maxD0Pv'),
        K0s_matchTrack2_maxD0Pv = ufloat('K0s_matchTrack2_maxD0Pv'),
        K0s_matchTrack1_minD0Pv = ufloat('K0s_matchTrack1_minD0Pv'),
        K0s_matchTrack2_minD0Pv = ufloat('K0s_matchTrack2_minD0Pv'),
        K0s_matchTrack1_dR  = ufloat('K0s_matchTrack1_dR'),
        K0s_matchTrack2_dR  = ufloat('K0s_matchTrack2_dR'),
        K0s_matchTrack1_pt  = ufloat('K0s_matchTrack1_pt'),
        K0s_matchTrack2_pt  = ufloat('K0s_matchTrack2_pt'),
        K0s_matchTrack1_eta = ufloat('K0s_matchTrack1_eta'),
        K0s_matchTrack2_eta = ufloat('K0s_matchTrack2_eta'),
        K0s_matchTrack1_phi = ufloat('K0s_matchTrack1_phi'),
        K0s_matchTrack2_phi = ufloat('K0s_matchTrack2_phi'),
        K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced = uint('K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced'),
        K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced = uint('K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced'),
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

candPiPiSequence = cms.Sequence(
    ( candPiPi * CountPiPi )
)

candPiPiTableSequence = cms.Sequence( candPiPiTable )

K0sToPiPiSequence = cms.Sequence(  
    ( K0sToPiPi * CountK0s )
)

K0sToPiPiTableSequence = cms.Sequence( K0sToPiPiTable )   

BToK0sMuMuPiPiSequence = cms.Sequence( 
    ( BToK0sMuMuPiPi * CountBToK0sMuMuPiPi )
)

BToK0sMuMuPiPiTableSequence = cms.Sequence( BToK0sMuMuPiPiTable )

