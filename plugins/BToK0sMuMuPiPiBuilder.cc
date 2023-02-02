////////////////////// Code to produce B0 -> LL Pi Pi K0s candidates /////////////////////////

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
//
#include "KinVtxFitterWithMassConstraint.h"
#include "KinVtxFitter.h"
#include "helper.h"
#include <limits>
#include <algorithm>
//
#include <vector>
#include <memory>
#include <map>
#include <string>
#include "TLorentzVector.h"

class BToK0sMuMuPiPiBuilder : public edm::global::EDProducer<> {

public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  
  explicit BToK0sMuMuPiPiBuilder(const edm::ParameterSet &cfg):
    
    // selections
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    post_vtx_selection2_{cfg.getParameter<std::string>("postVtxSelection2")},

    // inputs
    pfcands_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pfcands") )},
    
    dimuons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dimuons") )},
    muons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("muonTransientTracks") )},
    
    k0short_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("k0short") )},
    k0short_kinVtxsWMC_{consumes<std::vector<KinVtxFitterWithMassConstraint> >( cfg.getParameter<edm::InputTag>("k0shortKinVtxsWMC") )},
    
    dipions_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dipion") )},
    pions_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("pionTransientTracks") )},
    
    // wrt PV / BS
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )},
    vertex_src_{consumes<reco::VertexCollection>( cfg.getParameter<edm::InputTag>("offlinePrimaryVertexSrc") )}  
  {
    //output
    produces<pat::CompositeCandidateCollection>();
  }
  
  ~BToK0sMuMuPiPiBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:

  // selections
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection2_; 

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pfcands_; 

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuons_;
  const edm::EDGetTokenT<TransientTrackCollection> muons_ttracks_;

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> k0short_;
  const edm::EDGetTokenT<std::vector<KinVtxFitterWithMassConstraint> > k0short_kinVtxsWMC_;

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dipions_;
  const edm::EDGetTokenT<TransientTrackCollection> pions_ttracks_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
  const edm::EDGetTokenT<reco::VertexCollection> vertex_src_;
};

void BToK0sMuMuPiPiBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  // -------------------------------------------------------------
  // Masses - NB : should be coherent with those in helper!!
  float the_MUON_MASS = 0.10565837;
  float the_MUON_SIGMA = 0.0000001;
  float the_PI_MASS = 0.139571;
  float the_PI_SIGMA = 0.000016;
  float the_JPSI_MASS = 3.096900;
  // -------------------------------------------------------------
  
  // inputs
  edm::Handle<pat::CompositeCandidateCollection> pfcands;
  evt.getByToken(pfcands_, pfcands);  

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  evt.getByToken(dimuons_, dimuons);  
  edm::Handle<TransientTrackCollection> muons_ttracks;
  evt.getByToken(muons_ttracks_, muons_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> k0shorts;
  evt.getByToken(k0short_, k0shorts);  
  edm::Handle<std::vector<KinVtxFitterWithMassConstraint> > k0short_kinVtxsWMC;
  evt.getByToken(k0short_kinVtxsWMC_, k0short_kinVtxsWMC);

  edm::Handle<pat::CompositeCandidateCollection> dipions;
  evt.getByToken(dipions_, dipions);  
  edm::Handle<TransientTrackCollection> pions_ttracks;
  evt.getByToken(pions_ttracks_, pions_ttracks);

  // Others
  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  edm::Handle<reco::VertexCollection> pvtxs;
  evt.getByToken(vertex_src_, pvtxs);

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());


  // For kin fits
  KinematicParticleFactoryFromTransientTrack pFactory;
  
  // All k0s, muon pair and pipi pairs already passed cuts; no need for more preselection
  // Retrieve infos for mumu and pipi using the index associated to the K0s

  // Loop over K0s candidates
  for (size_t k0short_idx = 0; k0short_idx < k0shorts->size(); ++k0short_idx) {
    
    // this is the K0s
    edm::Ptr<pat::CompositeCandidate> k0short_ptr(k0shorts, k0short_idx);

    // this is the virtual particle after the K0s kin vertex fit with mass constraint
    KinVtxFitterWithMassConstraint thek0short_kinVtxsWMC = k0short_kinVtxsWMC->at(k0short_idx);

    // take the dimuon pair associated to this K0s
    int mumu_idx = k0short_ptr->userInt("mumu_idx");
    edm::Ptr<pat::CompositeCandidate> ll_ptr(dimuons, mumu_idx);
    int l1_idx = ll_ptr->userInt("l1_idx");
    int l2_idx = ll_ptr->userInt("l2_idx");

    // take the pi-pi pair associated to this K0s
    int pipi_idx = k0short_ptr->userInt("pipi_idx");
    edm::Ptr<pat::CompositeCandidate> pipi_ptr(dipions, pipi_idx);
    int trk1_idx = pipi_ptr->userInt("trk1_idx");
    int trk2_idx = pipi_ptr->userInt("trk2_idx"); 
    
   
    // B0 candidate
    pat::CompositeCandidate cand;	
    cand.setCharge( 0 );    // neutral - all pairs are neutral, already imposed when building pairs

    // save indices
    cand.addUserInt("mu1_idx", l1_idx);
    cand.addUserInt("mu2_idx", l2_idx);
    cand.addUserInt("trk1Rho_idx", trk1_idx);
    cand.addUserInt("trk2Rho_idx", trk2_idx);
    cand.addUserInt("k0short_idx", k0short_idx);
    cand.addUserInt("dipion_idx",  pipi_idx);
    cand.addUserInt("dimuon_idx",  mumu_idx);
    
    // Now make the B
    float chi = 0.;
    float ndf = 0.;
    std::vector<RefCountedKinematicParticle> B_candidate_init;
    B_candidate_init.push_back(pFactory.particle(muons_ttracks->at(l1_idx), the_MUON_MASS, chi, ndf, the_MUON_SIGMA));
    B_candidate_init.push_back(pFactory.particle(muons_ttracks->at(l2_idx), the_MUON_MASS, chi, ndf, the_MUON_SIGMA));
    B_candidate_init.push_back(pFactory.particle(pions_ttracks->at(trk1_idx), the_PI_MASS, chi, ndf, the_PI_SIGMA));
    B_candidate_init.push_back(pFactory.particle(pions_ttracks->at(trk2_idx), the_PI_MASS, chi, ndf, the_PI_SIGMA));
    B_candidate_init.push_back(thek0short_kinVtxsWMC.fitted_particle()); 


    // First fit, without J/Psi mass constraint
    std::vector<RefCountedKinematicParticle> B_candidate = B_candidate_init;
    KinematicParticleVertexFitter theFitter; 
    RefCountedKinematicTree vtx_tree = theFitter.fit(B_candidate);   
    if (vtx_tree->isEmpty() || !vtx_tree->isValid() || !vtx_tree->isConsistent()) continue;
    
    vtx_tree->movePointerToTheTop();
    float mass_womc = vtx_tree->currentParticle()->currentState().mass();
    cand.addUserFloat("fitted_mass_womc", mass_womc);   

    // Now refit with J/Psi mass constraint                   
    B_candidate = B_candidate_init;
    ParticleMass thePM_JPSI_MASS = the_JPSI_MASS;
    MultiTrackKinematicConstraint *ConstraintJpsiMass = new TwoTrackMassKinematicConstraint(thePM_JPSI_MASS);
    KinematicConstrainedVertexFitter theMCFitter; 
    RefCountedKinematicTree vertexFitTree = theMCFitter.fit(B_candidate, ConstraintJpsiMass);
    if (vertexFitTree->isEmpty() || !vertexFitTree->isValid() || !vertexFitTree->isConsistent()) continue;

    
    // Infos post final fit
    vertexFitTree->movePointerToTheTop();
    RefCountedKinematicParticle fitted_particle = vertexFitTree->currentParticle();     
    KinematicState fitted_candidate = fitted_particle->currentState();

    TLorentzVector fit_p4;
    fit_p4.SetPtEtaPhiM(fitted_candidate.globalMomentum().perp(), 
			fitted_candidate.globalMomentum().eta(),
			fitted_candidate.globalMomentum().phi(),
			mass_womc);
    cand.addUserFloat("fitted_pt"  , fit_p4.Perp()); 
    cand.addUserFloat("fitted_eta" , fit_p4.Eta());
    cand.addUserFloat("fitted_phi" , fit_p4.Phi());
    cand.addUserFloat("fitted_mass", fitted_particle->currentState().mass());    

    RefCountedKinematicVertex fitted_vtx = vertexFitTree->currentDecayVertex();     
    cand.setVertex( reco::Candidate::Point(fitted_vtx->position().x(), fitted_vtx->position().y(), fitted_vtx->position().z()) );
    cand.addUserFloat("sv_chi2", fitted_vtx->chiSquared());
    cand.addUserFloat("sv_prob", ChiSquaredProbability(fitted_vtx->chiSquared(),fitted_vtx->degreesOfFreedom()));

    // Now get children from final B fit                                                    
    vertexFitTree->movePointerToTheFirstChild();
    RefCountedKinematicParticle mu1cand = vertexFitTree->currentParticle();
    vertexFitTree->movePointerToTheNextChild();
    RefCountedKinematicParticle mu2cand = vertexFitTree->currentParticle();
    vertexFitTree->movePointerToTheNextChild();
    RefCountedKinematicParticle pi1cand = vertexFitTree->currentParticle();
    vertexFitTree->movePointerToTheNextChild();
    RefCountedKinematicParticle pi2cand = vertexFitTree->currentParticle();
    vertexFitTree->movePointerToTheNextChild();
    RefCountedKinematicParticle kscand = vertexFitTree->currentParticle();
    vertexFitTree->movePointerToTheTop();
    //                                                                                        
    GlobalVector mu1cand_p = mu1cand->currentState().kinematicParameters().momentum();
    GlobalVector mu2cand_p = mu2cand->currentState().kinematicParameters().momentum();
    GlobalVector pi1cand_p = pi1cand->currentState().kinematicParameters().momentum();
    GlobalVector pi2cand_p = pi2cand->currentState().kinematicParameters().momentum();
    GlobalVector k0scand_p = kscand->currentState().kinematicParameters().momentum();
    // 
    TLorentzVector p4fit_mu1, p4fit_mu2, p4fit_pi1, p4fit_pi2, p4fit_k0s;
    p4fit_mu1.SetXYZM(mu1cand_p.x(), mu1cand_p.y(), mu1cand_p.z(), MUON_MASS);
    p4fit_mu2.SetXYZM(mu2cand_p.x(), mu2cand_p.y(), mu2cand_p.z(), MUON_MASS);
    p4fit_pi1.SetXYZM(pi1cand_p.x(), pi1cand_p.y(), pi1cand_p.z(), PI_MASS);
    p4fit_pi2.SetXYZM(pi2cand_p.x(), pi2cand_p.y(), pi2cand_p.z(), PI_MASS);
    p4fit_k0s.SetXYZM(k0scand_p.x(), k0scand_p.y(), k0scand_p.z(), KSHORT_MASS);

    // Variables to be added to the final tree
    cand.addUserFloat("finalFit_X_mass",    ((p4fit_mu1 + p4fit_mu2 + p4fit_pi1 + p4fit_pi2).M()) );
    cand.addUserFloat("finalFit_Rho_mass",  ((p4fit_pi1 + p4fit_pi2).M()) );
    cand.addUserFloat("finalFit_JPsi_mass", ((p4fit_mu1 + p4fit_mu2).M()) );

    // post fit selection 
    if( !post_vtx_selection_(cand) ) continue;        



    // Refitted daughters - to be added to the final tree  
    cand.addUserFloat("finalFit_mu1_pt",  p4fit_mu1.Perp());
    cand.addUserFloat("finalFit_mu1_eta", p4fit_mu1.Eta());
    cand.addUserFloat("finalFit_mu1_phi", p4fit_mu1.Phi());
    cand.addUserFloat("finalFit_mu2_pt",  p4fit_mu2.Perp());
    cand.addUserFloat("finalFit_mu2_eta", p4fit_mu2.Eta());
    cand.addUserFloat("finalFit_mu2_phi", p4fit_mu2.Phi());
    //
    cand.addUserFloat("finalFit_pi1Rho_pt",  p4fit_pi1.Perp());
    cand.addUserFloat("finalFit_pi1Rho_eta", p4fit_pi1.Eta());
    cand.addUserFloat("finalFit_pi1Rho_phi", p4fit_pi1.Phi());
    cand.addUserFloat("finalFit_pi2Rho_pt",  p4fit_pi2.Perp());
    cand.addUserFloat("finalFit_pi2Rho_eta", p4fit_pi2.Eta());
    cand.addUserFloat("finalFit_pi2Rho_phi", p4fit_pi2.Phi());
    // 
    cand.addUserFloat("finalFit_k0s_pt",  p4fit_k0s.Perp()); 
    cand.addUserFloat("finalFit_k0s_eta", p4fit_k0s.Eta()); 
    cand.addUserFloat("finalFit_k0s_phi", p4fit_k0s.Phi()); 
    

    // Select the PV with best cosAlphaXYZ wrt the fitted B candidate
    int pv3D_idx = -1;
    Double_t lip3D = -1000.;
    math::XYZVector myBpt3D(fit_p4.Px(), fit_p4.Py(), fit_p4.Pz());
    //
    for (size_t vtx_idx = 0; vtx_idx < pvtxs->size(); ++vtx_idx) {
	  
      edm::Ptr<reco::Vertex> thisPV(pvtxs, vtx_idx);
      Double_t dx = fitted_vtx->position().x() - thisPV->x();
      Double_t dy = fitted_vtx->position().y() - thisPV->y();
      Double_t dz = fitted_vtx->position().z() - thisPV->z();

      math::XYZVector myDelta3D(dx, dy, dz);
      double myDen3D = myDelta3D.R()*myBpt3D.R();
      double cosAlphaXYZb = (myDelta3D.Dot(myBpt3D))/myDen3D;

      if (cosAlphaXYZb>lip3D) {
	lip3D = cosAlphaXYZb;
	pv3D_idx = vtx_idx;
      }

    } // Loop over PVs

    // chosen PV
    cand.addUserInt("pv3D_idx", pv3D_idx);	
    edm::Ptr<reco::Vertex> chosenPV(pvtxs, pv3D_idx);
    
    // CosAlpha 2Dim wrt beamspot without z-correction
    math::XYZVector myBpt2D(fit_p4.Px(), fit_p4.Py(), 0.);
    Double_t dxBS = fitted_vtx->position().x() - beamspot->x0();
    Double_t dyBS = fitted_vtx->position().y() - beamspot->y0();
    math::XYZVector myDeltaBS_2D(dxBS, dyBS, 0.);
    double myDenBS_2D = myDeltaBS_2D.R()*myBpt2D.R();
    double cosAlphaBS_XYb = (myDeltaBS_2D.Dot(myBpt2D))/myDenBS_2D;

    // CosAlpha 2Dim wrt beamspot with z-correction
    GlobalPoint the_fitted_vtx = fitted_vtx->position();
    math::XYZVector myDeltaBS_2Dwithz(the_fitted_vtx.x() - (beamspot->position(the_fitted_vtx.z())).x(), the_fitted_vtx.y() - (beamspot->position(the_fitted_vtx.z())).y(), 0.);  
    double myDenBS_2Dwithz = myDeltaBS_2Dwithz.R()*myBpt2D.R();
    double cosAlphaBS_XYbWithZ = (myDeltaBS_2Dwithz.Dot(myBpt2D))/myDenBS_2Dwithz;

    // Save cosAlpha with different algos
    cand.addUserFloat("cosAlpha3D_PV", lip3D);
    cand.addUserFloat("cosAlpha2D_BS", cosAlphaBS_XYb);
    cand.addUserFloat("cosAlpha2D_BSwithZ", cosAlphaBS_XYbWithZ);

    // Save BS coordinates
    cand.addUserFloat("BSxRaw", beamspot->x0());
    cand.addUserFloat("BSyRaw", beamspot->y0());
    cand.addUserFloat("BSxWithZ", (beamspot->position(the_fitted_vtx.z())).x());
    cand.addUserFloat("BSyWithZ", (beamspot->position(the_fitted_vtx.z())).y());

    // B vertex displacement significance wrt PV and BS
    float BSx   = beamspot->x0();
    float BSy   = beamspot->y0();
    float BSz   = 0.;
    float BSxE  = beamspot->covariance(0,0);
    float BSyE  = beamspot->covariance(1,1);
    float BSzE  = 0.;
    float BSxyE = beamspot->covariance(0,1);    
    //
    float chosenPVx   = chosenPV->x();
    float chosenPVy   = chosenPV->y();
    float chosenPVz   = 0.;
    float chosenPVxE  = chosenPV->covariance(0,0);
    float chosenPVyE  = chosenPV->covariance(1,1);
    float chosenPVzE  = 0.;
    float chosenPVxyE = chosenPV->covariance(0,1);
    // 
    float fittedVx    = fitted_vtx->position().x();
    float fittedVy    = fitted_vtx->position().y();
    float fittedVz    = 0.;
    float fittedVxE   = fitted_vtx->error().cxx();
    float fittedVyE   = fitted_vtx->error().cyy();
    float fittedVzE   = 0.;
    float fittedVxyE  = fitted_vtx->error().matrix()(0,1);

    // wrt PV chosen above
    float LxyB_PV = sqrt( (fittedVx-chosenPVx)*(fittedVx-chosenPVx) + (fittedVy-chosenPVy)*(fittedVy-chosenPVy) + (fittedVz-chosenPVz)*(fittedVz-chosenPVz) );
    float LxyBErr_PV = -1;
    if (LxyB_PV > 0.)
      LxyBErr_PV = sqrt((fittedVx-chosenPVx) * (fittedVx-chosenPVx) * fittedVxE +
			(fittedVy-chosenPVy) * (fittedVy-chosenPVy) * fittedVyE +
			(fittedVz-chosenPVz) * (fittedVz-chosenPVz) * fittedVzE +
		     
			(fittedVx-chosenPVx) * (fittedVy-chosenPVy) * 2.*fittedVxyE +
			(fittedVx-chosenPVx) * (fittedVz-chosenPVz) * 2.*0. +
			(fittedVy-chosenPVy) * (fittedVz-chosenPVz) * 2.*0. +
		     
			(fittedVx-chosenPVx) * (fittedVx-chosenPVx) * chosenPVxE +
			(fittedVy-chosenPVy) * (fittedVy-chosenPVy) * chosenPVyE +
			(fittedVz-chosenPVz) * (fittedVz-chosenPVz) * chosenPVzE +
		     
			(fittedVx-chosenPVx) * (fittedVy-chosenPVy) * 2.*chosenPVxyE +
			(fittedVx-chosenPVx) * (fittedVz-chosenPVz) * 2.*0. +
			(fittedVy-chosenPVy) * (fittedVz-chosenPVz) * 2.*0.) / LxyB_PV;

    // wrt BS without z correction
    float LxyB_BS = sqrt( (fittedVx-BSx)*(fittedVx-BSx) + (fittedVy-BSy)*(fittedVy-BSy) + (fittedVz-BSz)*(fittedVz-BSz) );
    float LxyBErr_BS = -1;
    if (LxyB_BS > 0.)
      LxyBErr_BS = sqrt((fittedVx-BSx) * (fittedVx-BSx) * fittedVxE +
			(fittedVy-BSy) * (fittedVy-BSy) * fittedVyE +
			(fittedVz-BSz) * (fittedVz-BSz) * fittedVzE +
		     
			(fittedVx-BSx) * (fittedVy-BSy) * 2.*fittedVxyE +
			(fittedVx-BSx) * (fittedVz-BSz) * 2.*0. +
			(fittedVy-BSy) * (fittedVz-BSz) * 2.*0. +
		     
			(fittedVx-BSx) * (fittedVx-BSx) * BSxE +
			(fittedVy-BSy) * (fittedVy-BSy) * BSyE +
			(fittedVz-BSz) * (fittedVz-BSz) * BSzE +
		     
			(fittedVx-BSx) * (fittedVy-BSy) * 2.*BSxyE +
			(fittedVx-BSx) * (fittedVz-BSz) * 2.*0. +
			(fittedVy-BSy) * (fittedVz-BSz) * 2.*0.) / LxyB_BS;
    
    // wrt BS with z correction (done as at HLT)
    GlobalError fitted_vtx_err( fitted_vtx->error().cxx(), fitted_vtx->error().cyx(), fitted_vtx->error().cyy(), fitted_vtx->error().czx(), fitted_vtx->error().czy(), fitted_vtx->error().czz() );
    GlobalPoint dispFromBS( -1*( (beamspot->x0() - fitted_vtx->position().x()) + (fitted_vtx->position().z() - beamspot->z0()) * beamspot->dxdz()), -1*((beamspot->y0() - fitted_vtx->position().y()) + (fitted_vtx->position().z() - beamspot->z0()) * beamspot->dydz()), 0);
    float LxyB_BSwithZ = dispFromBS.perp();
    float LxyBErr_BSwithZ = sqrt(fitted_vtx_err.rerr(dispFromBS));

    cand.addUserFloat("lxySign_PV", LxyB_PV/LxyBErr_PV);
    cand.addUserFloat("lxySign_BS", LxyB_BS/LxyBErr_BS); 
    cand.addUserFloat("lxySign_BSwithZ", LxyB_BSwithZ/LxyBErr_BSwithZ); 

    // post fit selection
    if( !post_vtx_selection2_(cand) ) continue;        


    // impact parameters wrt this vertex: muons
    GlobalPoint chosenPVpoint(chosenPV->position().x(), chosenPV->position().y(), chosenPV->position().z());
    
    TrajectoryStateClosestToPoint trajm1 = (muons_ttracks->at(l1_idx)).trajectoryStateClosestToPoint(chosenPVpoint);
    TrajectoryStateClosestToPoint trajm2 = (muons_ttracks->at(l2_idx)).trajectoryStateClosestToPoint(chosenPVpoint);
    float dzSignMu1  = trajm1.perigeeParameters().longitudinalImpactParameter()/trajm1.perigeeError().longitudinalImpactParameterError();
    float dzSignMu2  = trajm2.perigeeParameters().longitudinalImpactParameter()/trajm2.perigeeError().longitudinalImpactParameterError();
    float dxySignMu1 = trajm1.perigeeParameters().transverseImpactParameter()/trajm1.perigeeError().transverseImpactParameterError();
    float dxySignMu2 = trajm2.perigeeParameters().transverseImpactParameter()/trajm2.perigeeError().transverseImpactParameterError();
    
    // impact parameters wrt this vertex: pions
    TrajectoryStateClosestToPoint trajt1 = (pions_ttracks->at(trk1_idx)).trajectoryStateClosestToPoint(chosenPVpoint);
    TrajectoryStateClosestToPoint trajt2 = (pions_ttracks->at(trk2_idx)).trajectoryStateClosestToPoint(chosenPVpoint);
    float dzSignTr1  = trajt1.perigeeParameters().longitudinalImpactParameter()/trajt1.perigeeError().longitudinalImpactParameterError();
    float dzSignTr2  = trajt2.perigeeParameters().longitudinalImpactParameter()/trajt2.perigeeError().longitudinalImpactParameterError();
    float dxySignTr1 = trajt1.perigeeParameters().transverseImpactParameter()/trajt1.perigeeError().transverseImpactParameterError();
    float dxySignTr2 = trajt2.perigeeParameters().transverseImpactParameter()/trajt2.perigeeError().transverseImpactParameterError();


    // Save all wanted infos: B0-related
    cand.addUserFloat("decayVtxX", fitted_vtx->position().x());
    cand.addUserFloat("decayVtxY", fitted_vtx->position().y());
    cand.addUserFloat("decayVtxZ", fitted_vtx->position().z());
    cand.addUserFloat("decayVtxXE", fitted_vtx->error().cxx());
    cand.addUserFloat("decayVtxYE", fitted_vtx->error().cyy()); 
    cand.addUserFloat("decayVtxZE", fitted_vtx->error().czz()); 

    // Save all wanted infos: MuMu-related
    cand.addUserFloat("MuMu_sv_prob", ll_ptr->userFloat("sv_prob"));
    cand.addUserFloat("MuMu_fitted_mass", ll_ptr->userFloat("fitted_mass"));
    cand.addUserFloat("MuMu_fitted_pt",  ll_ptr->userFloat("fitted_pt"));
    cand.addUserFloat("MuMu_fitted_eta", ll_ptr->userFloat("fitted_eta"));
    cand.addUserFloat("MuMu_fitted_phi", ll_ptr->userFloat("fitted_phi"));
    cand.addUserFloat("MuMu_fitted_vtxX", ll_ptr->userFloat("fitted_vtxX"));
    cand.addUserFloat("MuMu_fitted_vtxY", ll_ptr->userFloat("fitted_vtxY"));
    cand.addUserFloat("MuMu_fitted_vtxZ", ll_ptr->userFloat("fitted_vtxZ"));
    cand.addUserFloat("MuMu_fitted_vtxXE", ll_ptr->userFloat("fitted_vtxEx"));
    cand.addUserFloat("MuMu_fitted_vtxYE", ll_ptr->userFloat("fitted_vtxEy"));
    cand.addUserFloat("MuMu_fitted_vtxZE", ll_ptr->userFloat("fitted_vtxEz"));
    cand.addUserFloat("MuMu_prefit_mu1_pt",  ll_ptr->userFloat("mu1_pt"));
    cand.addUserFloat("MuMu_prefit_mu1_eta", ll_ptr->userFloat("mu1_eta"));
    cand.addUserFloat("MuMu_prefit_mu1_phi", ll_ptr->userFloat("mu1_phi"));
    cand.addUserFloat("MuMu_prefit_mu2_pt",  ll_ptr->userFloat("mu2_pt"));
    cand.addUserFloat("MuMu_prefit_mu2_eta", ll_ptr->userFloat("mu2_eta"));
    cand.addUserFloat("MuMu_prefit_mu2_phi", ll_ptr->userFloat("mu2_phi"));
    cand.addUserFloat("MuMu_mu1_dr", ll_ptr->userFloat("mu1_dr")); 
    cand.addUserFloat("MuMu_mu2_dr", ll_ptr->userFloat("mu2_dr")); 
    cand.addUserFloat("MuMu_DCA", ll_ptr->userFloat("DCA")); 
    cand.addUserFloat("MuMu_LxySign", ll_ptr->userFloat("LxySign"));
    cand.addUserFloat("MuMu_cosAlpha", ll_ptr->userFloat("cosAlpha"));
    cand.addUserFloat("MuMu_mu1_dzsign",  dzSignMu1);
    cand.addUserFloat("MuMu_mu2_dzsign",  dzSignMu2);
    cand.addUserFloat("MuMu_mu1_dxysign", dxySignMu1);
    cand.addUserFloat("MuMu_mu2_dxysign", dxySignMu2);
    cand.addUserInt("MuMu_mu1_trackQuality", ll_ptr->userInt("mu1_trackQuality")); 
    cand.addUserInt("MuMu_mu2_trackQuality", ll_ptr->userInt("mu2_trackQuality")); 
    cand.addUserInt("MuMu_mu1_fired_Dimuon20_Jpsi",               ll_ptr->userInt("mu1_fired_Dimuon20_Jpsi"));
    cand.addUserInt("MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced", ll_ptr->userInt("mu1_fired_DoubleMu4_JpsiTrk_Displaced"));
    cand.addUserInt("MuMu_mu2_fired_Dimuon20_Jpsi",               ll_ptr->userInt("mu2_fired_Dimuon20_Jpsi"));
    cand.addUserInt("MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced", ll_ptr->userInt("mu2_fired_DoubleMu4_JpsiTrk_Displaced"));
    cand.addUserFloat("MuMu_mu1_dr_Dimuon20_Jpsi",                ll_ptr->userFloat("mu1_dr_Dimuon20_Jpsi"));
    cand.addUserFloat("MuMu_mu1_dr_DoubleMu4_JpsiTrk_Displaced",  ll_ptr->userFloat("mu1_dr_DoubleMu4_JpsiTrk_Displaced"));
    cand.addUserFloat("MuMu_mu2_dr_Dimuon20_Jpsi",                ll_ptr->userFloat("mu2_dr_Dimuon20_Jpsi"));
    cand.addUserFloat("MuMu_mu2_dr_DoubleMu4_JpsiTrk_Displaced",  ll_ptr->userFloat("mu2_dr_DoubleMu4_JpsiTrk_Displaced"));

    // Save all wanted infos: PiPi (from Rho)-related
    cand.addUserFloat("PiPi_prefit_pi1_pt",  pipi_ptr->userFloat("pi1_pt"));
    cand.addUserFloat("PiPi_prefit_pi1_eta", pipi_ptr->userFloat("pi1_eta"));
    cand.addUserFloat("PiPi_prefit_pi1_phi", pipi_ptr->userFloat("pi1_phi"));
    cand.addUserFloat("PiPi_prefit_pi1_vx",  pipi_ptr->userFloat("pi1_vx"));
    cand.addUserFloat("PiPi_prefit_pi1_vy",  pipi_ptr->userFloat("pi1_vy"));
    cand.addUserFloat("PiPi_prefit_pi1_vz",  pipi_ptr->userFloat("pi1_vz"));
    cand.addUserFloat("PiPi_pi1_d0sig",      pipi_ptr->userFloat("pi1_d0sig"));
    cand.addUserFloat("PiPi_pi1_maxd0PV",    pipi_ptr->userFloat("pi1_maxd0PV")); 
    cand.addUserFloat("PiPi_pi1_mind0PV",    pipi_ptr->userFloat("pi1_mind0PV")); 
    cand.addUserFloat("PiPi_prefit_pi2_pt",  pipi_ptr->userFloat("pi2_pt"));
    cand.addUserFloat("PiPi_prefit_pi2_eta", pipi_ptr->userFloat("pi2_eta"));
    cand.addUserFloat("PiPi_prefit_pi2_phi", pipi_ptr->userFloat("pi2_phi"));
    cand.addUserFloat("PiPi_prefit_pi2_vx",  pipi_ptr->userFloat("pi2_vx"));
    cand.addUserFloat("PiPi_prefit_pi2_vy",  pipi_ptr->userFloat("pi2_vy"));
    cand.addUserFloat("PiPi_prefit_pi2_vz",  pipi_ptr->userFloat("pi2_vz"));
    cand.addUserFloat("PiPi_pi2_d0sig",      pipi_ptr->userFloat("pi2_d0sig"));
    cand.addUserFloat("PiPi_pi2_maxd0PV",    pipi_ptr->userFloat("pi2_maxd0PV")); 
    cand.addUserFloat("PiPi_pi2_mind0PV",    pipi_ptr->userFloat("pi2_mind0PV")); 
    cand.addUserFloat("PiPi_pi1_dzsign",  dzSignTr1);
    cand.addUserFloat("PiPi_pi2_dzsign",  dzSignTr2);
    cand.addUserFloat("PiPi_pi1_dxysign", dxySignTr1);
    cand.addUserFloat("PiPi_pi2_dxysign", dxySignTr2);
    cand.addUserInt("PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced", pipi_ptr->userInt("p1_fired_DoubleMu4_JpsiTrk_Displaced"));
    cand.addUserInt("PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced", pipi_ptr->userInt("p2_fired_DoubleMu4_JpsiTrk_Displaced"));
    cand.addUserFloat("PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced",  pipi_ptr->userFloat("p1_dr_DoubleMu4_JpsiTrk_Displaced"));
    cand.addUserFloat("PiPi_p2_dr_DoubleMu4_JpsiTrk_Displaced",  pipi_ptr->userFloat("p2_dr_DoubleMu4_JpsiTrk_Displaced"));
    cand.addUserFloat("PiPi_sv_prob", pipi_ptr->userFloat("sv_prob"));    

    // Save all wanted infos: K0s-related
    cand.addUserFloat("K0s_prefit_mass",      k0short_ptr->userFloat("prefit_mass"));
    cand.addUserFloat("K0s_nmcFitted_mass",   k0short_ptr->userFloat("fitted_nmc_mass"));
    cand.addUserFloat("K0s_nmcFitted_pi1pt",  k0short_ptr->userFloat("fitted_nmc_pi1pt"));
    cand.addUserFloat("K0s_nmcFitted_pi1eta", k0short_ptr->userFloat("fitted_nmc_pi1eta"));
    cand.addUserFloat("K0s_nmcFitted_pi1phi", k0short_ptr->userFloat("fitted_nmc_pi1phi"));
    cand.addUserFloat("K0s_nmcFitted_pi2pt",  k0short_ptr->userFloat("fitted_nmc_pi2pt"));
    cand.addUserFloat("K0s_nmcFitted_pi2eta", k0short_ptr->userFloat("fitted_nmc_pi2eta"));
    cand.addUserFloat("K0s_nmcFitted_pi2phi", k0short_ptr->userFloat("fitted_nmc_pi2phi"));
    //
    cand.addUserFloat("K0s_mcFitted_svprob", k0short_ptr->userFloat("sv_prob"));
    cand.addUserFloat("K0s_mcFitted_mass",   k0short_ptr->userFloat("fitted_mass"));
    cand.addUserFloat("K0s_mcFitted_pt",     k0short_ptr->userFloat("fitted_pt"));
    cand.addUserFloat("K0s_mcFitted_eta",    k0short_ptr->userFloat("fitted_eta"));
    cand.addUserFloat("K0s_mcFitted_phi",    k0short_ptr->userFloat("fitted_phi"));
    cand.addUserFloat("K0s_mcFitted_vtxX",   k0short_ptr->userFloat("fitted_vtxX"));
    cand.addUserFloat("K0s_mcFitted_vtxY",   k0short_ptr->userFloat("fitted_vtxY"));
    cand.addUserFloat("K0s_mcFitted_vtxZ",   k0short_ptr->userFloat("fitted_vtxZ"));
    cand.addUserFloat("K0s_mcFitted_vtxXE",  k0short_ptr->userFloat("fitted_vtxEx"));
    cand.addUserFloat("K0s_mcFitted_vtxYE",  k0short_ptr->userFloat("fitted_vtxEy"));
    cand.addUserFloat("K0s_mcFitted_vtxZE",  k0short_ptr->userFloat("fitted_vtxEz"));

    
    // Compute the K0s displacement wrt the B vertex
    float k0sPVx   = k0short_ptr->userFloat("fitted_vtxX");
    float k0sPVy   = k0short_ptr->userFloat("fitted_vtxY");
    float k0sPVz   = 0.;
    float k0sPVxE  = k0short_ptr->userFloat("fitted_vtxEx");
    float k0sPVyE  = k0short_ptr->userFloat("fitted_vtxEy");
    float k0sPVzE  = 0.;
    float k0sPVxyE = k0short_ptr->userFloat("fitted_vtxExy");
    
    float LxyKs = sqrt( (fittedVx-k0sPVx)*(fittedVx-k0sPVx) + (fittedVy-k0sPVy)*(fittedVy-k0sPVy) + (fittedVz-k0sPVz)*(fittedVz-k0sPVz) );
    float LxyKsErr = -1;
    if (LxyKs > 0.)
      LxyKsErr = sqrt((fittedVx-k0sPVx) * (fittedVx-k0sPVx) * fittedVxE +
		      (fittedVy-k0sPVy) * (fittedVy-k0sPVy) * fittedVyE +
		      (fittedVz-k0sPVz) * (fittedVz-k0sPVz) * fittedVzE +
		      
		      (fittedVx-k0sPVx) * (fittedVy-k0sPVy) * 2.*fittedVxyE +
		      (fittedVx-k0sPVx) * (fittedVz-k0sPVz) * 2.*0. +
		      (fittedVy-k0sPVy) * (fittedVz-k0sPVz) * 2.*0. +
		      
		      (fittedVx-k0sPVx) * (fittedVx-k0sPVx) * k0sPVxE +
		      (fittedVy-k0sPVy) * (fittedVy-k0sPVy) * k0sPVyE +
		      (fittedVz-k0sPVz) * (fittedVz-k0sPVz) * k0sPVzE +
		      
		      (fittedVx-k0sPVx) * (fittedVy-k0sPVy) * 2.*k0sPVxyE +
		      (fittedVx-k0sPVx) * (fittedVz-k0sPVz) * 2.*0. +
		      (fittedVy-k0sPVy) * (fittedVz-k0sPVz) * 2.*0.) / LxyKs;
    
    cand.addUserFloat("K0_lxySign_wrtBvtx", (LxyKs/LxyKsErr));


    // Compute cosAlpha between K0 vertex, B0 vertex and K0 momentum
    TLorentzVector p4K0s;
    p4K0s.SetPtEtaPhiM(cand.userFloat("K0s_mcFitted_pt"), cand.userFloat("K0s_mcFitted_eta"), cand.userFloat("K0s_mcFitted_phi"), KSHORT_MASS);
    math::XYZVector myK0pt3D(p4K0s.Px(), p4K0s.Py(), p4K0s.Pz());
    //
    Double_t dxK0s = fitted_vtx->position().x() - k0sPVx;
    Double_t dyK0s = fitted_vtx->position().y() - k0sPVy;
    Double_t dzK0s = fitted_vtx->position().z() - k0sPVz;
    //
    math::XYZVector myK0BDelta3D(dxK0s, dyK0s, dzK0s);
    double myK0BDen3D = myK0BDelta3D.R()*myK0pt3D.R();
    double cosAlphaK0BXYZ = (myK0BDelta3D.Dot(myK0pt3D))/myK0BDen3D;
    //
    cand.addUserFloat("K0_cosAlpha3D", cosAlphaK0BXYZ);     


    // To emulate the trigger: look for tracks in the collection 
    // matching pi1 and pi2 from the selected k0s
    float minDr1 = 1000.;
    float minDr2 = 1000.;
    float matchD0sign1  = -1.;
    float matchD0sign2  = -1.;
    float matchMaxD0Pv1 = -1.;
    float matchMaxD0Pv2 = -1.;
    float matchMinD0Pv1 = -1.;
    float matchMinD0Pv2 = -1.;
    float matchPt1  = -1.;
    float matchPt2  = -1.;
    float matchEta1 = -1.;
    float matchEta2 = -1.;
    float matchPhi1 = -1.;
    float matchPhi2 = -1.;
    int fired_DoubleMu4_JpsiTrk_Displaced_1     = -1;
    int fired_DoubleMu4_JpsiTrk_Displaced_2     = -1;
    
    TLorentzVector p4_K0sP1, p4_K0sP2;
    p4_K0sP1.SetPtEtaPhiM(cand.userFloat("K0s_nmcFitted_pi1pt"), cand.userFloat("K0s_nmcFitted_pi1eta"), cand.userFloat("K0s_nmcFitted_pi1phi"), PI_MASS);
    p4_K0sP2.SetPtEtaPhiM(cand.userFloat("K0s_nmcFitted_pi2pt"), cand.userFloat("K0s_nmcFitted_pi2eta"), cand.userFloat("K0s_nmcFitted_pi2phi"), PI_MASS);

    for(size_t trk_idx = 0; trk_idx < pfcands->size(); ++trk_idx ){
      edm::Ptr<pat::CompositeCandidate> trk_ptr( pfcands, trk_idx );
      TLorentzVector thetrk;
      thetrk.SetPtEtaPhiM(trk_ptr->pt(), trk_ptr->eta(), trk_ptr->phi(), PI_MASS);
      float dr1 = thetrk.DeltaR(p4_K0sP1);
      float dr2 = thetrk.DeltaR(p4_K0sP2);
      if (dr1<minDr1) {
	minDr1 = dr1;
	matchD0sign1  = trk_ptr->userFloat("d0sig");
	matchMaxD0Pv1 = trk_ptr->userFloat("maxd0PV");   
	matchMinD0Pv1 = trk_ptr->userFloat("mind0PV");   
	matchPt1  = trk_ptr->pt();
	matchEta1 = trk_ptr->eta();
	matchPhi1 = trk_ptr->phi();
	fired_DoubleMu4_JpsiTrk_Displaced_1     = trk_ptr->userInt("HLT_DoubleMu4_JpsiTrk_Displaced");
      }
      if (dr2<minDr2) {
	minDr2 = dr2;
	matchD0sign2  = trk_ptr->userFloat("d0sig");
	matchMaxD0Pv2 = trk_ptr->userFloat("maxd0PV");   
	matchMinD0Pv2 = trk_ptr->userFloat("mind0PV");   
	matchPt2  = trk_ptr->pt();
	matchEta2 = trk_ptr->eta();
	matchPhi2 = trk_ptr->phi();
	fired_DoubleMu4_JpsiTrk_Displaced_2     = trk_ptr->userInt("HLT_DoubleMu4_JpsiTrk_Displaced");
      }
    }
    cand.addUserFloat("K0s_matchTrack1_D0sign",  matchD0sign1);
    cand.addUserFloat("K0s_matchTrack2_D0sign",  matchD0sign2);
    cand.addUserFloat("K0s_matchTrack1_maxD0Pv", matchMaxD0Pv1);
    cand.addUserFloat("K0s_matchTrack2_maxD0Pv", matchMaxD0Pv2);
    cand.addUserFloat("K0s_matchTrack1_minD0Pv", matchMinD0Pv1);
    cand.addUserFloat("K0s_matchTrack2_minD0Pv", matchMinD0Pv2);
    cand.addUserFloat("K0s_matchTrack1_dR",      minDr1);
    cand.addUserFloat("K0s_matchTrack2_dR",      minDr2);  
    cand.addUserFloat("K0s_matchTrack1_pt",      matchPt1);
    cand.addUserFloat("K0s_matchTrack2_pt",      matchPt2);  
    cand.addUserFloat("K0s_matchTrack1_eta",     matchEta1);
    cand.addUserFloat("K0s_matchTrack2_eta",     matchEta2);
    cand.addUserFloat("K0s_matchTrack1_phi",     matchPhi1);
    cand.addUserFloat("K0s_matchTrack2_phi",     matchPhi2);
    cand.addUserInt("K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced", fired_DoubleMu4_JpsiTrk_Displaced_1);
    cand.addUserInt("K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced", fired_DoubleMu4_JpsiTrk_Displaced_2);

    // Save all wanted infos: PV
    cand.addUserFloat("PVx", chosenPVx);
    cand.addUserFloat("PVy", chosenPVy);
    cand.addUserFloat("PVz", chosenPVz);
    cand.addUserFloat("PVEx", chosenPVxE);
    cand.addUserFloat("PVEy", chosenPVyE);
    cand.addUserFloat("PVEz", chosenPVzE);

    // Save output 
    ret_val->push_back(cand);
    
  }  // K0s Loop

  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToK0sMuMuPiPiBuilder);
