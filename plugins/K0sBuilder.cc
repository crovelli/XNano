///////////////////////// Code to produce K0s -> pi pi candidates ////////////////////////

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"

#include <vector>
#include <string>
#include "TLorentzVector.h"

#include "helper.h"
#include "KinVtxFitter.h"
#include "KinVtxFitterWithMassConstraint.h"

#include "TVectorD.h"    // for fixing tracks
#include "TMatrixDSym.h" // for fixing tracks

class K0sBuilder : public edm::global::EDProducer<> {

public:

  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit K0sBuilder(const edm::ParameterSet &cfg):

    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},

    transientTrackRecordToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
    dimuons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dimuons") )},    
    dipions_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dipion") )},
    svs_{consumes<reco::VertexCompositePtrCandidateCollection>( cfg.getParameter<edm::InputTag>("svSrc") )} {
      
      // output
      produces<pat::CompositeCandidateCollection>("SelectedK0s");
      produces<std::vector<KinVtxFitterWithMassConstraint> >("SelectedK0sKinVtxsWMC");
      produces<std::vector<KinVtxFitter> >("SelectedK0sKinVtxsNMC");
    }

  ~K0sBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

private:

  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; 

  reco::Track fix_track(const reco::Track *tk, double delta) const;
  
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackRecordToken_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuons_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dipions_;
  const edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svs_;
};


void K0sBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &iSetup ) const {

  // inputs  
  auto const& theB = iSetup.getHandle(transientTrackRecordToken_);
  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  evt.getByToken(dimuons_, dimuons);  
  edm::Handle<pat::CompositeCandidateCollection> dipions;
  evt.getByToken(dipions_, dipions);  
  edm::Handle<reco::VertexCompositePtrCandidateCollection> v0Coll;
  evt.getByToken(svs_, v0Coll);

  // output 
  std::unique_ptr<pat::CompositeCandidateCollection> k0s_out(new pat::CompositeCandidateCollection());
  std::unique_ptr<std::vector<KinVtxFitterWithMassConstraint> > kinVtx_out( new std::vector<KinVtxFitterWithMassConstraint> );  


  // main loop
  for (reco::VertexCompositePtrCandidateCollection::const_iterator vee = v0Coll->begin(); vee != v0Coll->end(); vee++) {

    // Order pions in decreasing pT
    reco::TransientTrack pi1TT; 
    reco::TransientTrack pi2TT; 
    /*
    if (vee->daughter(0)->momentum().mag2() < vee->daughter(1)->momentum().mag2()) {
      pi1TT = theB->build( vee->daughter(0)->bestTrack() );  
      pi2TT = theB->build( vee->daughter(1)->bestTrack() );
    } else {
      pi1TT = theB->build( vee->daughter(1)->bestTrack() );
      pi2TT = theB->build( vee->daughter(0)->bestTrack() );
    }
    */
    if (vee->daughter(0)->momentum().mag2() < vee->daughter(1)->momentum().mag2()) {
      pi1TT = theB->build( fix_track( vee->daughter(0)->bestTrack(), 1e-8 ) );  
      pi2TT = theB->build( fix_track( vee->daughter(1)->bestTrack(), 1e-8 ) );
    } else {
      pi1TT = theB->build( fix_track( (vee->daughter(1)->bestTrack()), 1e-8) );
      pi2TT = theB->build( fix_track( (vee->daughter(0)->bestTrack()), 1e-8) );
    }
    if(!pi1TT.isValid()) continue;
    if(!pi2TT.isValid()) continue;

    // K0s candidate
    pat::CompositeCandidate k0s_tempcand;
    float prefit_mass = vee->mass();
    

    // Vertex fit without mass constraint
    KinVtxFitter fitterNMC(
			   {pi1TT, pi2TT},
			   {PI_MASS, PI_MASS },
			   {PI_SIGMA, PI_SIGMA} 
			   );
    if ( !fitterNMC.success() ) continue;           
    float fitted_mass_womc = fitterNMC.fitted_candidate().mass();
    k0s_tempcand.addUserFloat("fitted_mass_womc", fitted_mass_womc);

    // Vertex fit refinement, adding mass constraint
    KinVtxFitterWithMassConstraint fitter(KSHORT_MASS, KSHORT_SIGMA, fitterNMC.vtx_tree());
    if ( !fitter.success() ) continue;   
    
    // Save quantities after fit
    TLorentzVector p4k0s;
    p4k0s.SetPtEtaPhiM(fitter.fitted_candidate().globalMomentum().perp(), 
		       fitter.fitted_candidate().globalMomentum().eta(),
		       fitter.fitted_candidate().globalMomentum().phi(),
		       KSHORT_MASS);
    RefCountedKinematicVertex fitted_vtx = fitter.fitted_refvtx();

    // Selection after fit 
    k0s_tempcand.addUserFloat("sv_prob", fitter.prob());    
    k0s_tempcand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass() );
    if( !post_vtx_selection_(k0s_tempcand) ) continue;

    
    
    // Now check compatibility with selected dimuons and dipions
    // Loop over pipi candidates
    for (size_t pipi_idx = 0; pipi_idx < dipions->size(); ++pipi_idx) {
      edm::Ptr<pat::CompositeCandidate> pipi_ptr(dipions, pipi_idx);
      TLorentzVector p4tr1, p4tr2;
      p4tr1.SetPtEtaPhiM(pipi_ptr->userFloat("pi1_pt"), pipi_ptr->userFloat("pi1_eta"), pipi_ptr->userFloat("pi1_phi"), PI_MASS);
      p4tr2.SetPtEtaPhiM(pipi_ptr->userFloat("pi2_pt"), pipi_ptr->userFloat("pi2_eta"), pipi_ptr->userFloat("pi2_phi"), PI_MASS);
      
      // Now take the mumu index corresponding to this pair
      int mumu_idx = pipi_ptr->userInt("mumu_idx");
      edm::Ptr<pat::CompositeCandidate> ll_ptr(dimuons, mumu_idx);
      TLorentzVector p4mu1, p4mu2;
      p4mu1.SetPtEtaPhiM(ll_ptr->userFloat("mu1_pt"), ll_ptr->userFloat("mu1_eta"), ll_ptr->userFloat("mu1_phi"), MUON_MASS);
      p4mu2.SetPtEtaPhiM(ll_ptr->userFloat("mu2_pt"), ll_ptr->userFloat("mu2_eta"), ll_ptr->userFloat("mu2_phi"), MUON_MASS);

      // Hardcoded - compatibility checks
      if ((p4mu1 + p4mu2 + p4tr1 + p4tr2 + p4k0s).M() - (p4mu1 + p4mu2).M() + JPSI_MASS > 5.9) continue;
      if ((p4mu1 + p4mu2 + p4tr1 + p4tr2 + p4k0s).M() - (p4mu1 + p4mu2).M() + JPSI_MASS < 4.9) continue;

      // For selected K0s, create a candidate to be saved in the event
      pat::CompositeCandidate k0s_cand;
      k0s_cand.setP4(vee->polarP4());
      if (vee->daughter(0)->momentum().mag2() < vee->daughter(1)->momentum().mag2()) {
	k0s_cand.addUserCand("trk1", vee->daughterPtr(0)); 
	k0s_cand.addUserCand("trk2", vee->daughterPtr(1));
      } else {
	k0s_cand.addUserCand("trk1", vee->daughterPtr(1)); 
	k0s_cand.addUserCand("trk2", vee->daughterPtr(0));
      }

     
      // Infos after first fit:
      k0s_cand.addUserFloat("fitted_nmc_mass", fitterNMC.fitted_candidate().mass() );   
      math::PtEtaPhiMLorentzVector pi1TLV = fitterNMC.daughter_p4(0);                   
      math::PtEtaPhiMLorentzVector pi2TLV = fitterNMC.daughter_p4(1);                   
      k0s_cand.addUserFloat("fitted_nmc_pi1pt",  pi1TLV.Pt());
      k0s_cand.addUserFloat("fitted_nmc_pi1eta", pi1TLV.Eta());
      k0s_cand.addUserFloat("fitted_nmc_pi1phi", pi1TLV.Phi());
      k0s_cand.addUserFloat("fitted_nmc_pi2pt",  pi2TLV.Pt());
      k0s_cand.addUserFloat("fitted_nmc_pi2eta", pi2TLV.Eta());
      k0s_cand.addUserFloat("fitted_nmc_pi2phi", pi2TLV.Phi());
      
      // Infos after final fit
      k0s_cand.addUserFloat("sv_prob", fitter.prob());    
      k0s_cand.addUserFloat("prefit_mass", prefit_mass);  
      k0s_cand.addUserFloat("fitted_mass_womc", fitted_mass_womc);
      k0s_cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass() );
      k0s_cand.addUserFloat("fitted_pt",  p4k0s.Pt());
      k0s_cand.addUserFloat("fitted_eta", p4k0s.Eta());     
      k0s_cand.addUserFloat("fitted_phi", p4k0s.Phi()); 
      k0s_cand.addUserFloat("fitted_vtxX",  fitted_vtx->position().x());
      k0s_cand.addUserFloat("fitted_vtxY",  fitted_vtx->position().y());
      k0s_cand.addUserFloat("fitted_vtxZ",  fitted_vtx->position().z());
      k0s_cand.addUserFloat("fitted_vtxEx", fitted_vtx->error().cxx());
      k0s_cand.addUserFloat("fitted_vtxEy", fitted_vtx->error().cyy());
      k0s_cand.addUserFloat("fitted_vtxEz", fitted_vtx->error().czz());
      k0s_cand.addUserFloat("fitted_vtxExy", fitted_vtx->error().matrix()(0,1));

      // Index
      k0s_cand.addUserInt("mumu_idx", mumu_idx );
      k0s_cand.addUserInt("pipi_idx", pipi_idx );
	
      // Put in the event
      k0s_out->push_back(k0s_cand);  
      kinVtx_out->push_back(fitter);  
      
    }   // Loop over rhos
  }     // Loop over Koshort
  
  evt.put(std::move(k0s_out), "SelectedK0s");
  evt.put(std::move(kinVtx_out), "SelectedK0sKinVtxsWMC");   
}

// O. Cerri's code to deal with not positive definite covariance matrices
// https://github.com/ocerri/BPH_RDntuplizer/blob/master/plugins/VtxUtils.cc
/* Check for a not positive definite covariance matrix. If the covariance matrix is not positive definite, we force it to be positive definite by
 * adding the minimum eigenvalue to the diagonal of the covariance matrix plus `delta`.
 * See https://nhigham.com/2020/12/22/what-is-a-modified-cholesky-factorization/ */

reco::Track K0sBuilder::fix_track(const reco::Track *tk, double delta) const {

  unsigned int i, j;
  double min_eig = 1;

  // Get the original covariance matrix. 
  reco::TrackBase::CovarianceMatrix cov = tk->covariance();

  // Convert it from an SMatrix to a TMatrixD so we can get the eigenvalues. 
  TMatrixDSym new_cov(cov.kRows);
  for (i = 0; i < cov.kRows; i++) {
    for (j = 0; j < cov.kRows; j++) {
    // Need to check for nan or inf, because for some reason these
    // cause a segfault when calling Eigenvectors().
    //
    // No idea what to do here or why this happens. 
    if (std::isnan(cov(i,j)) || std::isinf(cov(i,j)))
	cov(i,j) = 1e-6;
      new_cov(i,j) = cov(i,j);
    }
  }

  // Get the eigenvalues. 
  TVectorD eig(cov.kRows);
  new_cov.EigenVectors(eig);
  for (i = 0; i < cov.kRows; i++)
    if (eig(i) < min_eig)
      min_eig = eig(i);

  // If the minimum eigenvalue is less than zero, then subtract it from the
  // diagonal and add `delta`. 
  if (min_eig < 0) {
    for (i = 0; i < cov.kRows; i++)
      cov(i,i) -= min_eig - delta;
  }

  return reco::Track(tk->chi2(), tk->ndof(), tk->referencePoint(), tk->momentum(), tk->charge(), cov, tk->algo(), (reco::TrackBase::TrackQuality) tk->qualityMask());
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(K0sBuilder);

