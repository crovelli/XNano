///////////////////////// Code to produce pi pi candidates ////////////////////////

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

class PiPiBuilder : public edm::global::EDProducer<> {

public:

  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  
  explicit PiPiBuilder(const edm::ParameterSet &cfg):
    trk1_selection_{cfg.getParameter<std::string>("trk1Selection")},
    trk2_selection_{cfg.getParameter<std::string>("trk2Selection")}, 
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    pfcands_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pfcands") )},
    ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracks") )} {
      produces<pat::CompositeCandidateCollection>("SelectedDiPions");
      produces<std::vector<KinVtxFitter> >("SelectedDiPionKinVtxs");  
    }

  ~PiPiBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> trk1_selection_;     
  const StringCutObjectSelector<pat::CompositeCandidate> trk2_selection_;     
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_;
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pfcands_; 
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_; 
};


void PiPiBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  // inputs  
  edm::Handle<pat::CompositeCandidateCollection> pfcands;
  evt.getByToken(pfcands_, pfcands);  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_, ttracks);
  
  // output
  std::unique_ptr<pat::CompositeCandidateCollection> pipi_out(new pat::CompositeCandidateCollection());
  std::unique_ptr<std::vector<KinVtxFitter> > kinVtx_out( new std::vector<KinVtxFitter> );   
  
  // main loop
  for(size_t trk1_idx = 0; trk1_idx < pfcands->size(); ++trk1_idx ){
    edm::Ptr<pat::CompositeCandidate> trk1_ptr( pfcands, trk1_idx );
    if(!trk1_selection_(*trk1_ptr)) continue; 
    
    for(size_t trk2_idx = trk1_idx + 1; trk2_idx < pfcands->size(); ++trk2_idx) {
      edm::Ptr<pat::CompositeCandidate> trk2_ptr( pfcands, trk2_idx );
      if(!trk2_selection_(*trk2_ptr)) continue;

      // opposite charge
      if (trk1_ptr->charge() == trk2_ptr->charge()) continue; 

      // create a pi+ / pi- pair
      pat::CompositeCandidate pipi_cand;
      auto trk1_p4=trk1_ptr->polarP4();
      auto trk2_p4=trk2_ptr->polarP4();
      trk1_p4.SetM(PI_MASS);
      trk2_p4.SetM(PI_MASS);

      // adding stuff for pre fit selection
      pipi_cand.setP4(trk1_p4 + trk2_p4);
      pipi_cand.addUserFloat("trk_deltaR", reco::deltaR(*trk1_ptr, *trk2_ptr));

      // Save candidates
      pipi_cand.addUserCand("trk1", trk1_ptr );
      pipi_cand.addUserCand("trk2", trk2_ptr );

      // save indices
      pipi_cand.addUserInt("trk1_idx", trk1_idx );
      pipi_cand.addUserInt("trk2_idx", trk2_idx );
      
      // selection before fit
      if( !pre_vtx_selection_(pipi_cand) ) continue;

      KinVtxFitter fitter(
			  {ttracks->at(trk1_idx), ttracks->at(trk2_idx)},
			  {PI_MASS, PI_MASS },
			  {PI_SIGMA, PI_SIGMA} 
			  );
      if ( !fitter.success() ) continue;           
      
      // save quantities after fit
      pipi_cand.addUserFloat("sv_chi2", fitter.chi2());
      pipi_cand.addUserFloat("sv_ndof", fitter.dof()); 
      pipi_cand.addUserFloat("sv_prob", fitter.prob());    
      pipi_cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass());

      // refitted daughters (tracks)     
      // std::vector<std::string> dnames{ "trk1", "trk2" };
      // for (size_t idaughter=0; idaughter<dnames.size(); idaughter++){
      // pipi_cand.addUserFloat("fitted_" + dnames[idaughter] + "_pt" ,fitter.daughter_p4(idaughter).pt() );
      // pipi_cand.addUserFloat("fitted_" + dnames[idaughter] + "_eta",fitter.daughter_p4(idaughter).eta() );
      // pipi_cand.addUserFloat("fitted_" + dnames[idaughter] + "_phi",fitter.daughter_p4(idaughter).phi() );
      // }

      // after fit selection 
      if( !post_vtx_selection_(pipi_cand) ) continue;

      pipi_out->push_back(pipi_cand);
      kinVtx_out->push_back(fitter);   
    }
  }
  
  evt.put(std::move(pipi_out), "SelectedDiPions");
  evt.put(std::move(kinVtx_out), "SelectedDiPionKinVtxs");   
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PiPiBuilder);

