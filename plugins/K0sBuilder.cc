///////////////////////// Code to produce K0s -> pi pi candidates ////////////////////////

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
#include "TLorentzVector.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
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
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"
#include "KinVtxFitterWithMassConstraint.h"

class K0sBuilder : public edm::global::EDProducer<> {

public:

  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  
  explicit K0sBuilder(const edm::ParameterSet &cfg):
    
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    svs_{consumes<reco::VertexCompositePtrCandidateCollection>( cfg.getParameter<edm::InputTag>("svSrc") )} {
      
      //output
      produces<pat::CompositeCandidateCollection>("SelectedK0s");
      produces<std::vector<KinVtxFitterWithMassConstraint> >("SelectedK0sKinVtxsWMC");
    }

  ~K0sBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; 
  const edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svs_;
};


void K0sBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &iSetup) const {

  // inputs  
  edm::Handle<reco::VertexCompositePtrCandidateCollection> v0Coll;
  evt.getByToken(svs_, v0Coll);
  
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
  
  // output
  std::unique_ptr<pat::CompositeCandidateCollection> k0s_out(new pat::CompositeCandidateCollection());
  std::unique_ptr<std::vector<KinVtxFitterWithMassConstraint> > kinVtx_out( new std::vector<KinVtxFitterWithMassConstraint> );

  // main loop
  for (reco::VertexCompositePtrCandidateCollection::const_iterator vee = v0Coll->begin(); vee != v0Coll->end(); vee++) {
    
    // Order pions in decreasing pT
    reco::TransientTrack pi1TT; 
    reco::TransientTrack pi2TT; 
    if (vee->daughter(0)->momentum().mag2() < vee->daughter(1)->momentum().mag2()) {
      pi1TT = theB->build( vee->daughter(0)->bestTrack() );  
      pi2TT = theB->build( vee->daughter(1)->bestTrack() );
    } else {
      pi1TT = theB->build( vee->daughter(1)->bestTrack() );
      pi2TT = theB->build( vee->daughter(0)->bestTrack() );
    }
    if(!pi1TT.isValid()) continue;
    if(!pi2TT.isValid()) continue;
    
    // K0s Candidate
    pat::CompositeCandidate k0s_cand;
    k0s_cand.addUserFloat("prefit_mass", vee->mass()); 

    // Vertex fit without mass constraint
    KinVtxFitter fitterNMC(
			   {pi1TT, pi2TT},
			   {PI_MASS, PI_MASS },
			   {PI_SIGMA, PI_SIGMA} 
			   );
    if ( !fitterNMC.success() ) continue;      
    k0s_cand.addUserFloat("fitted_mass_womc", fitterNMC.fitted_candidate().mass());   

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
    k0s_cand.addUserFloat("sv_prob", fitter.prob());    
    k0s_cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass() );
    if( !post_vtx_selection_(k0s_cand) ) continue;


    // For selected K0s add more infos
    k0s_cand.setP4(vee->polarP4());
    if (vee->daughter(0)->momentum().mag2() < vee->daughter(1)->momentum().mag2()) {
      k0s_cand.addUserFloat("trk1_pt",   vee->daughterPtr(0)->bestTrack()->pt());
      k0s_cand.addUserFloat("trk1_eta",  vee->daughterPtr(0)->bestTrack()->eta());
      k0s_cand.addUserFloat("trk1_phi",  vee->daughterPtr(0)->bestTrack()->phi());
      k0s_cand.addUserInt("trk1_charge", vee->daughterPtr(0)->bestTrack()->charge()); 
      k0s_cand.addUserFloat("trk2_pt",   vee->daughterPtr(1)->bestTrack()->pt());
      k0s_cand.addUserFloat("trk2_eta",  vee->daughterPtr(1)->bestTrack()->eta());
      k0s_cand.addUserFloat("trk2_phi",  vee->daughterPtr(1)->bestTrack()->phi());
      k0s_cand.addUserInt("trk2_charge", vee->daughterPtr(1)->bestTrack()->charge()); 
    } else {
      k0s_cand.addUserFloat("trk1_pt",   vee->daughterPtr(1)->bestTrack()->pt());
      k0s_cand.addUserFloat("trk1_eta",  vee->daughterPtr(1)->bestTrack()->eta());
      k0s_cand.addUserFloat("trk1_phi",  vee->daughterPtr(1)->bestTrack()->phi());
      k0s_cand.addUserInt("trk1_charge", vee->daughterPtr(1)->bestTrack()->charge()); 
      k0s_cand.addUserFloat("trk2_pt",   vee->daughterPtr(0)->bestTrack()->pt());
      k0s_cand.addUserFloat("trk2_eta",  vee->daughterPtr(0)->bestTrack()->eta());
      k0s_cand.addUserFloat("trk2_phi",  vee->daughterPtr(0)->bestTrack()->phi());
      k0s_cand.addUserInt("trk2_charge", vee->daughterPtr(0)->bestTrack()->charge()); 
    }
    
    // save quantities after fit
    k0s_cand.addUserFloat("fitted_pt", fitter.fitted_candidate().globalMomentum().perp() );
    k0s_cand.addUserFloat("fitted_eta", fitter.fitted_candidate().globalMomentum().eta() );
    k0s_cand.addUserFloat("fitted_phi", fitter.fitted_candidate().globalMomentum().phi() );
    
    // Put in the event
    k0s_out->push_back(k0s_cand);  
    kinVtx_out->push_back(fitter);

  } // Loop over V0s
  
  evt.put(std::move(k0s_out), "SelectedK0s");
  evt.put(std::move(kinVtx_out), "SelectedK0sKinVtxsWMC");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(K0sBuilder);

