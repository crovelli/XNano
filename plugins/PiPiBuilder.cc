///////////////////////// Code to produce pi pi candidates ////////////////////////

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include <vector>
#include <string>
#include "TLorentzVector.h"

#include "helper.h"

class PiPiBuilder : public edm::global::EDProducer<> {

public:

  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  
  explicit PiPiBuilder(const edm::ParameterSet &cfg):
    trk_selection_{cfg.getParameter<std::string>("trkSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    dimuons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dimuons") )},
    pfcands_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pfcands") )} {
      produces<pat::CompositeCandidateCollection>("SelectedDiPions");
    }

  ~PiPiBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> trk_selection_;     
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuons_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pfcands_; 
};

void PiPiBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  // inputs  
  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  evt.getByToken(dimuons_, dimuons);  
  edm::Handle<pat::CompositeCandidateCollection> pfcands;
  evt.getByToken(pfcands_, pfcands);  
  
  // output
  std::unique_ptr<pat::CompositeCandidateCollection> pipi_out(new pat::CompositeCandidateCollection());
  

  // Start from muons, to reduce the combinatorics
  for(size_t ll_idx = 0; ll_idx < dimuons->size(); ++ll_idx) {   
    
    edm::Ptr<pat::CompositeCandidate> ll_ptr(dimuons, ll_idx);
    
    TLorentzVector p4mu1, p4mu2;
    p4mu1.SetPtEtaPhiM(ll_ptr->userFloat("mu1_pt"), ll_ptr->userFloat("mu1_eta"), ll_ptr->userFloat("mu1_phi"), MUON_MASS);
    p4mu2.SetPtEtaPhiM(ll_ptr->userFloat("mu2_pt"), ll_ptr->userFloat("mu2_eta"), ll_ptr->userFloat("mu2_phi"), MUON_MASS);

    // Loop over pions, I
    for(size_t trk1_idx = 0; trk1_idx < pfcands->size(); ++trk1_idx ){
      edm::Ptr<pat::CompositeCandidate> trk1_ptr( pfcands, trk1_idx );
      if(!trk_selection_(*trk1_ptr)) continue; 
      
      TLorentzVector trk1;
      trk1.SetPtEtaPhiM(trk1_ptr->pt(), trk1_ptr->eta(), trk1_ptr->phi(), PI_MASS);

      // hardcoded...
      if (trk1.DeltaR(p4mu1) < 0.01) continue;
      if (trk1.DeltaR(p4mu2) < 0.01) continue;
      if ((p4mu1 + p4mu2 + trk1).M() - (p4mu1 + p4mu2).M() + JPSI_MASS > 4.0) continue; 


      // Loop over pions, II
      for(size_t trk2_idx = trk1_idx + 1; trk2_idx < pfcands->size(); ++trk2_idx) {
	edm::Ptr<pat::CompositeCandidate> trk2_ptr( pfcands, trk2_idx );
	if(!trk_selection_(*trk2_ptr)) continue;
	
	// opposite charge
	if (trk1_ptr->charge() == trk2_ptr->charge()) continue; 

	TLorentzVector trk2;
	trk2.SetPtEtaPhiM(trk2_ptr->pt(), trk2_ptr->eta(), trk2_ptr->phi(), PI_MASS);
	
	// hardcoded...
	if (trk2.DeltaR(p4mu1) < 0.01) continue;
	if (trk2.DeltaR(p4mu2) < 0.01) continue;
	if ((p4mu1 + p4mu2 + trk1 + trk2).M() - (p4mu1 + p4mu2).M() + JPSI_MASS > 4.1) continue; 


	// Now I've a pair of pions which matches the two selected muons
	pat::CompositeCandidate pipi_cand;
	auto trk1_p4=trk1_ptr->polarP4();
	auto trk2_p4=trk2_ptr->polarP4();
	trk1_p4.SetM(PI_MASS);
	trk2_p4.SetM(PI_MASS);

	// adding stuff for pre fit selection
	pipi_cand.setP4(trk1_p4 + trk2_p4);
	pipi_cand.addUserFloat("trk_deltaR", reco::deltaR(*trk1_ptr, *trk2_ptr));
	
	// Save candidates and indices
	pipi_cand.addUserCand("trk1", trk1_ptr );
	pipi_cand.addUserCand("trk2", trk2_ptr );
	pipi_cand.addUserInt("trk1_idx", trk1_idx );
	pipi_cand.addUserInt("trk2_idx", trk2_idx );

	// minimal selection 
	if( !pre_vtx_selection_(pipi_cand) ) continue;

	// Saver further quantities to be saved in the final ntuples
	pipi_cand.addUserInt("mumu_idx", ll_idx );  

	// save further quantities, to be saved in the final ntuples: pions before fit
	// Pions post fit are saved only after the very final B fit
	pipi_cand.addUserFloat("pi1_pt",  trk1_ptr->pt()); 
	pipi_cand.addUserFloat("pi1_eta", trk1_ptr->eta());
	pipi_cand.addUserFloat("pi1_phi", trk1_ptr->phi());
	pipi_cand.addUserFloat("pi2_pt",  trk2_ptr->pt());
	pipi_cand.addUserFloat("pi2_eta", trk2_ptr->eta());
	pipi_cand.addUserFloat("pi2_phi", trk2_ptr->phi());
	pipi_cand.addUserFloat("pi1_vx",  trk1_ptr->vx()); 
	pipi_cand.addUserFloat("pi1_vy",  trk1_ptr->vy()); 
	pipi_cand.addUserFloat("pi1_vz",  trk1_ptr->vz()); 
	pipi_cand.addUserFloat("pi2_vx",  trk2_ptr->vx()); 
	pipi_cand.addUserFloat("pi2_vy",  trk2_ptr->vy()); 
	pipi_cand.addUserFloat("pi2_vz",  trk2_ptr->vz()); 

	// Put in the event
	pipi_out->push_back(pipi_cand);

      } // Loop 1st pi
    } // Loop 2nd pi
  } // Loop muons  

  evt.put(std::move(pipi_out), "SelectedDiPions");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PiPiBuilder);

