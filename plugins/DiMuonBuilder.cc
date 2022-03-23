#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"

#include <vector>
#include <string>

#include "helper.h"
#include "KinVtxFitter.h"

class DiMuonBuilder : public edm::global::EDProducer<> {

public:

  typedef std::vector<pat::Muon> MuonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit DiMuonBuilder(const edm::ParameterSet &cfg):
    l1_selection_{cfg.getParameter<std::string>("lep1Selection")},
    l2_selection_{cfg.getParameter<std::string>("lep2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    src_{consumes<MuonCollection>( cfg.getParameter<edm::InputTag>("src") )},
    ttracks_src_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracksSrc") )} {
      produces<pat::CompositeCandidateCollection>("SelectedDiMuons");
    }

  ~DiMuonBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::Muon> l1_selection_;    
  const StringCutObjectSelector<pat::Muon> l2_selection_;    
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_;
  const edm::EDGetTokenT<MuonCollection> src_;
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_src_;
};

void DiMuonBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {
  
  // input
  edm::Handle<MuonCollection> muons;
  evt.getByToken(src_, muons);
  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_src_, ttracks);

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  
  for(size_t l1_idx = 0; l1_idx < muons->size(); ++l1_idx) {
    edm::Ptr<pat::Muon> l1_ptr(muons, l1_idx);
    if(!l1_selection_(*l1_ptr)) continue; 
    
    for(size_t l2_idx = l1_idx + 1; l2_idx < muons->size(); ++l2_idx) {
      edm::Ptr<pat::Muon> l2_ptr(muons, l2_idx);
      if(!l2_selection_(*l2_ptr)) continue;


      // Muons must be different
      if (l1_idx==l2_idx) continue;

      pat::CompositeCandidate muon_pair;
      muon_pair.setP4(l1_ptr->p4() + l2_ptr->p4());
      muon_pair.setCharge(l1_ptr->charge() + l2_ptr->charge());
      muon_pair.addUserFloat("lep_deltaR", reco::deltaR(*l1_ptr, *l2_ptr));
      
      // Put the lepton passing the corresponding selection
      muon_pair.addUserInt("l1_idx", l1_idx );
      muon_pair.addUserInt("l2_idx", l2_idx );

      // Use UserCands as they should not use memory but keep the Ptr itself
      muon_pair.addUserCand("l1", l1_ptr );
      muon_pair.addUserCand("l2", l2_ptr );

      if( !pre_vtx_selection_(muon_pair) ) continue; // before making the SV, cut on the info we have

      KinVtxFitter fitter(
			  {ttracks->at(l1_idx), ttracks->at(l2_idx)},
			  {l1_ptr->mass(), l2_ptr->mass()},
			  {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
			  );
      if ( !fitter.success() ) continue;

      // save quantities after fit needed for selection and to be saved in the final ntuples
      KinematicState fitted_cand = fitter.fitted_candidate();

      muon_pair.addUserFloat("sv_prob", fitter.prob());
      muon_pair.addUserFloat("fitted_mass", fitter.success() ? fitted_cand.mass() : -1);

      // cut on the SV info
      if( !post_vtx_selection_(muon_pair) ) continue;


      // save further quantities, to be saved in the final ntuples: JPsi infos after fit
      TVector3 B_J(fitted_cand.globalMomentum().x(),
		   fitted_cand.globalMomentum().y(),
		   fitted_cand.globalMomentum().z());
      muon_pair.addUserFloat("fitted_pt",  B_J.Pt());
      muon_pair.addUserFloat("fitted_eta", B_J.Eta());
      muon_pair.addUserFloat("fitted_phi", B_J.Phi());

      // save further quantities, to be saved in the final ntuples: JPsi vertex after fit
      RefCountedKinematicVertex fitted_vtx = fitter.fitted_refvtx();
      muon_pair.addUserFloat("fitted_vtxX",  fitted_vtx->position().x());
      muon_pair.addUserFloat("fitted_vtxY",  fitted_vtx->position().y());
      muon_pair.addUserFloat("fitted_vtxZ",  fitted_vtx->position().z());
      muon_pair.addUserFloat("fitted_vtxEx", fitted_vtx->error().cxx());
      muon_pair.addUserFloat("fitted_vtxEy", fitted_vtx->error().cyy());
      muon_pair.addUserFloat("fitted_vtxEz", fitted_vtx->error().czz());

      // save further quantities, to be saved in the final ntuples: muons before fit
      // Muons post fit are saved only after the very final B fit
      muon_pair.addUserFloat("mu1_pt",  l1_ptr->pt());
      muon_pair.addUserFloat("mu1_eta", l1_ptr->eta());
      muon_pair.addUserFloat("mu1_phi", l1_ptr->phi());
      muon_pair.addUserFloat("mu2_pt",  l2_ptr->pt());
      muon_pair.addUserFloat("mu2_eta", l2_ptr->eta());
      muon_pair.addUserFloat("mu2_phi", l2_ptr->phi());

      
      // push in the event
      ret_value->push_back(muon_pair);
    }
  }

  evt.put(std::move(ret_value),  "SelectedDiMuons");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DiMuonBuilder);
