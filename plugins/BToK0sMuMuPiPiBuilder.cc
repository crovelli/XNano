////////////////////// Code to produce B0 -> LL Pi Pi Pi Pi candidates /////////////////////////

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

class BToK0sMuMuPiPiBuilder : public edm::global::EDProducer<> {

public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  
  explicit BToK0sMuMuPiPiBuilder(const edm::ParameterSet &cfg):
    
    // selections
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},

    // inputs
    dimuons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dimuons") )},
    dimuons_kinVtxs_{consumes<std::vector<KinVtxFitter> >( cfg.getParameter<edm::InputTag>("dimuonKinVtxs") )},
    muons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("muonTransientTracks") )},

    k0short_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("k0short") )},
    k0short_kinVtxs_{consumes<std::vector<KinVtxFitter> >( cfg.getParameter<edm::InputTag>("k0shortKinVtxs") )},
    k0short_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("k0shortTransientTracks") )},

    dipions_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dipion") )},
    dipions_kinVtxs_{consumes<std::vector<KinVtxFitter> >( cfg.getParameter<edm::InputTag>("dipionKinVtxs") )},
    pions_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("pionTransientTracks") )},

    // wrt PV / BS
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )},
    vertex_src_{consumes<reco::VertexCollection>( cfg.getParameter<edm::InputTag>("offlinePrimaryVertexSrc") )}

    //isotracksToken_{consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))},
    //isolostTracksToken_{consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))},
    //isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")}
    {
      //output
      produces<pat::CompositeCandidateCollection>();
    }
  
  ~BToK0sMuMuPiPiBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:

  // selections
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_;  
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; 


  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuons_;
  const edm::EDGetTokenT<std::vector<KinVtxFitter> > dimuons_kinVtxs_;
  const edm::EDGetTokenT<TransientTrackCollection> muons_ttracks_;

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> k0short_;
  const edm::EDGetTokenT<std::vector<KinVtxFitter> > k0short_kinVtxs_;
  const edm::EDGetTokenT<TransientTrackCollection> k0short_ttracks_;

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dipions_;
  const edm::EDGetTokenT<std::vector<KinVtxFitter> > dipions_kinVtxs_;
  const edm::EDGetTokenT<TransientTrackCollection> pions_ttracks_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
  const edm::EDGetTokenT<reco::VertexCollection> vertex_src_;

  //const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  //const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;
  //const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 
};

void BToK0sMuMuPiPiBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  // inputs
  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  evt.getByToken(dimuons_, dimuons);  
  edm::Handle<std::vector<KinVtxFitter> > dimuons_kinVtxs;
  evt.getByToken(dimuons_kinVtxs_, dimuons_kinVtxs);
  edm::Handle<TransientTrackCollection> muons_ttracks;
  evt.getByToken(muons_ttracks_, muons_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> k0shorts;
  evt.getByToken(k0short_, k0shorts);  
  edm::Handle<std::vector<KinVtxFitter> > k0short_kinVtxs;
  evt.getByToken(k0short_kinVtxs_, k0short_kinVtxs);
  edm::Handle<TransientTrackCollection> k0short_ttracks;
  evt.getByToken(k0short_ttracks_, k0short_ttracks);   

  edm::Handle<pat::CompositeCandidateCollection> dipions;
  evt.getByToken(dipions_, dipions);  
  edm::Handle<std::vector<KinVtxFitter> > dipions_kinVtxs;
  evt.getByToken(dipions_kinVtxs_, dipions_kinVtxs);  
  edm::Handle<TransientTrackCollection> pions_ttracks;
  evt.getByToken(pions_ttracks_, pions_ttracks);

  // Others
  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  edm::Handle<reco::VertexCollection> pvtxs;
  evt.getByToken(vertex_src_, pvtxs);

  //for isolation
  //edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  //evt.getByToken(isotracksToken_, iso_tracks);
  //edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  //evt.getByToken(isolostTracksToken_, iso_lostTracks);
  //unsigned int nTracks     = iso_tracks->size();
  //unsigned int totalTracks = nTracks + iso_lostTracks->size();

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  
  // All k0s, muon pair and pipi pais already passed cuts; no need for more preselection

  // Loop over pi-pi candidates from X
  for (size_t pipi_idx = 0; pipi_idx < dipions->size(); ++pipi_idx) {

    edm::Ptr<pat::CompositeCandidate> pipi_ptr(dipions, pipi_idx);
    edm::Ptr<reco::Candidate> trk1_ptr= pipi_ptr->userCand("trk1");
    edm::Ptr<reco::Candidate> trk2_ptr= pipi_ptr->userCand("trk2");
    int trk1_idx = pipi_ptr->userInt("trk1_idx");
    int trk2_idx = pipi_ptr->userInt("trk2_idx");

    // Loop over mu-mu candidates from X    
    for(size_t ll_idx = 0; ll_idx < dimuons->size(); ++ll_idx) {
      
      edm::Ptr<pat::CompositeCandidate> ll_ptr(dimuons, ll_idx);
      edm::Ptr<reco::Candidate> l1_ptr = ll_ptr->userCand("l1");
      edm::Ptr<reco::Candidate> l2_ptr = ll_ptr->userCand("l2");
      int l1_idx = ll_ptr->userInt("l1_idx");
      int l2_idx = ll_ptr->userInt("l2_idx");

      // Loop over K0s candidates
      for (size_t k0short_idx = 0; k0short_idx < k0shorts->size(); ++k0short_idx) {

	edm::Ptr<pat::CompositeCandidate> k0short_ptr(k0shorts, k0short_idx);
	edm::Ptr<reco::Candidate> trk1k0s_ptr= k0short_ptr->userCand("trk1");
	edm::Ptr<reco::Candidate> trk2k0s_ptr= k0short_ptr->userCand("trk2");
	int trk1k0s_idx = k0short_ptr->userInt("trk1_idx");
	int trk2k0s_idx = k0short_ptr->userInt("trk2_idx");

	// We need 4 different tracks
	if ( (trk1k0s_idx==trk1_idx) || (trk1k0s_idx==trk2_idx) ) continue;
	if ( (trk2k0s_idx==trk1_idx) || (trk2k0s_idx==trk2_idx) ) continue;

	// B0 candidate
	pat::CompositeCandidate cand;
	cand.setP4( l1_ptr->p4() + l2_ptr->p4() + trk1k0s_ptr->p4() + trk2k0s_ptr->p4() + trk1_ptr->p4() + trk2_ptr->p4() );
	cand.setCharge( 0 );    // neutral

	// save daughters - unfitted
	cand.addUserCand("mu1", l1_ptr);
	cand.addUserCand("mu2", l2_ptr);
	cand.addUserCand("trk1Rho", trk1_ptr);
	cand.addUserCand("trk2Rho", trk2_ptr);
	cand.addUserCand("trk1k0s", trk1k0s_ptr);
	cand.addUserCand("trk2k0s", trk2k0s_ptr);
	cand.addUserCand("k0short", k0short_ptr);
	cand.addUserCand("dipion",  pipi_ptr);
	cand.addUserCand("dilepton", ll_ptr);

	// save indices
	cand.addUserInt("mu1_idx", l1_idx);
	cand.addUserInt("mu2_idx", l2_idx);
	cand.addUserInt("trk1Rho_idx", trk1_idx);
	cand.addUserInt("trk2Rho_idx", trk2_idx);
	cand.addUserInt("trk1k0s_idx", trk1k0s_idx);
	cand.addUserInt("trk2k0s_idx", trk2k0s_idx);
	cand.addUserInt("k0short_idx", k0short_idx);
	cand.addUserInt("dipion_idx",  pipi_idx);
	cand.addUserInt("dilepton_idx", ll_idx);

	auto dr_info = min_max_dr({l1_ptr, l2_ptr, trk1_ptr, trk2_ptr, trk1k0s_ptr, trk2k0s_ptr});  
	cand.addUserFloat("min_dr", dr_info.first);   
	cand.addUserFloat("max_dr", dr_info.second); 

	// check if pass pre vertex cut
	if( !pre_vtx_selection_(cand) ) continue;

	KinVtxFitter fitter( 
			    { k0short_ttracks->at(trk1k0s_idx), k0short_ttracks->at(trk2k0s_idx), pions_ttracks->at(trk1_idx), pions_ttracks->at(trk2_idx), muons_ttracks->at(l1_idx), muons_ttracks->at(l2_idx) },
			    {PI_MASS, PI_MASS, PI_MASS, PI_MASS, l1_ptr->mass(), l2_ptr->mass()},
			    {PI_SIGMA, PI_SIGMA, PI_SIGMA, PI_SIGMA, LEP_SIGMA, LEP_SIGMA}  
			     );
	if(!fitter.success()) continue; 


	// B0 position
	cand.setVertex( 
		       reco::Candidate::Point( 
					      fitter.fitted_vtx().x(),
					      fitter.fitted_vtx().y(),
					      fitter.fitted_vtx().z()
					       )  
			);

	// vertex vars
	cand.addUserFloat("sv_chi2", fitter.chi2());
	cand.addUserFloat("sv_ndof", fitter.dof());
	cand.addUserFloat("sv_prob", fitter.prob());

	// refitted kinematic vars
	cand.addUserFloat("fitted_k0short_mass",(fitter.daughter_p4(0) + fitter.daughter_p4(1)).mass() );
	cand.addUserFloat("fitted_Rho_mass",    (fitter.daughter_p4(2) + fitter.daughter_p4(3)).mass() );
	cand.addUserFloat("fitted_JPsi_mass"   ,(fitter.daughter_p4(4) + fitter.daughter_p4(5)).mass());
	cand.addUserFloat("fitted_X_mass"      ,(fitter.daughter_p4(2) + fitter.daughter_p4(3) + fitter.daughter_p4(4) + fitter.daughter_p4(5)).mass());

	auto fit_p4 = fitter.fitted_p4();
	GlobalPoint fit_pos = fitter.fitted_pos();

	cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
	cand.addUserFloat("fitted_eta" , fit_p4.eta());
	cand.addUserFloat("fitted_phi" , fit_p4.phi());
	cand.addUserFloat("fitted_mass", fit_p4.mass());      
	cand.addUserFloat("fitted_massErr", sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6))); 

	// refitted daughters (leptons/tracks)     
	std::vector<std::string> dnames{ "trk1k0s", "trk2k0s", "trk1Rho", "trk2Rho", "mu1", "mu2" };
      	
	for (size_t idaughter=0; idaughter<dnames.size(); idaughter++){
	  cand.addUserFloat("fitted_" + dnames[idaughter] + "_pt" ,fitter.daughter_p4(idaughter).pt() );
	  cand.addUserFloat("fitted_" + dnames[idaughter] + "_eta",fitter.daughter_p4(idaughter).eta() );
	  cand.addUserFloat("fitted_" + dnames[idaughter] + "_phi",fitter.daughter_p4(idaughter).phi() );
	}      

	// Select the PV with best cosAlphaXYb wrt the fitted B candidate
	int pv_idx = -1;
	Double_t lip = -1000.;
	for (size_t vtx_idx = 0; vtx_idx < pvtxs->size(); ++vtx_idx) {
	  
	  edm::Ptr<reco::Vertex> thisPV(pvtxs, vtx_idx);
	  Double_t dx = cand.vx() - thisPV->x();
	  Double_t dy = cand.vy() - thisPV->y();
	  Double_t dz = cand.vz() - thisPV->z();
	  Double_t cosAlphaXYb = ( fit_pos.x() * dx + fit_pos.y()*dy + fit_pos.z()*dz  )/( sqrt(dx*dx+dy*dy+dz*dz)* fit_pos.mag() );

	  if (cosAlphaXYb>lip) {
	    lip = cosAlphaXYb;
	    pv_idx = vtx_idx;
	  }
	} // Loop over PVs


	// chosen PV
	cand.addUserInt("pv_idx", pv_idx);	
	edm::Ptr<reco::Vertex> chosenPV(pvtxs, pv_idx);

	// other vars: B direction wrt distance PV-SV
	cand.addUserFloat(
			  "cosTheta2D_BS", 
			  cos_theta_2D(fitter, *beamspot, cand.p4())
			  );
	
	cand.addUserFloat(
			  "cosTheta2D_PV", 
			  cos_theta_2D(fitter, *chosenPV, cand.p4())
			  );
	
	cand.addUserFloat(
			  "fitted_cosTheta2D_BS", 
			  cos_theta_2D(fitter, *beamspot, fit_p4)
			  );

	cand.addUserFloat(
			  "fitted_cosTheta2D_PV", 
			  cos_theta_2D(fitter, *chosenPV, fit_p4)
			  );
	
	// B vertex displacement significance 
	auto lxy_BS = l_xy(fitter, *beamspot);
	auto lxy_PV = l_xy(fitter, *chosenPV);
	cand.addUserFloat("lxy_BS", lxy_BS.value());
	cand.addUserFloat("lxy_PV", lxy_PV.value());
	cand.addUserFloat("lxyUnc_BS", lxy_BS.error());
	cand.addUserFloat("lxyUnc_PV", lxy_PV.error());
	cand.addUserFloat("lxySign_PV", (lxy_PV.value()/lxy_PV.error()));
	
	// post fit selection
	if( !post_vtx_selection_(cand) ) continue;        
	


	// chiara: second fit - to be removed, pick-up just one
	KinVtxFitter fitter2( 
			     { dimuons_kinVtxs->at(ll_idx).fitted_candidate_ttrk(), k0short_kinVtxs->at(k0short_idx).fitted_candidate_ttrk(), dipions_kinVtxs->at(pipi_idx).fitted_candidate_ttrk() },
			     { JPSI_MASS, 0.497611, 0.77526},
			     { JPSI_SIGMA, KSHORT_SIGMA, RHO_SIGMA }
			      );
	if(!fitter2.success()) continue; 
	
	auto fit2_p4 = fitter2.fitted_p4();
	cand.addUserFloat("fitted2_mass", fit2_p4.mass());      
	// chiara: second fit - to be removed, pick-up just one 
	

	/*
	//compute isolation
	float l1_iso03  = 0;
	float l1_iso04  = 0;
	float l2_iso03  = 0;
	float l2_iso04  = 0;
	float tk1_iso03 = 0;
	float tk1_iso04 = 0;
	float tk2_iso03 = 0;
	float tk2_iso04 = 0;
	float b_iso03   = 0;
	float b_iso04   = 0;
	
	for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
	  
	  const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
	  // define selections for iso tracks (pT, eta, ...)
	  if( !isotrk_selection_(trk) ) continue;
	  // check if the track is the kaon or the pion
	  if (trk1k0s_ptr ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
	  if (trk2k0s_ptr ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
	  // check if the track is one of the two leptons
	  if (track_to_lepton_match(l1_ptr, iso_tracks.id(), iTrk) || 
	      track_to_lepton_match(l2_ptr, iso_tracks.id(), iTrk) ) continue;
	  
	  // add to final particle iso if dR < cone
	  float dr_to_l1  = deltaR(cand.userFloat("fitted_l1_eta"),  cand.userFloat("fitted_l1_phi"),  trk.eta(), trk.phi());
	  float dr_to_l2  = deltaR(cand.userFloat("fitted_l2_eta"),  cand.userFloat("fitted_l2_phi"),  trk.eta(), trk.phi());
	  float dr_to_tk1 = deltaR(cand.userFloat("fitted_trk1_eta"),cand.userFloat("fitted_trk1_phi"),trk.eta(), trk.phi());
	  float dr_to_tk2 = deltaR(cand.userFloat("fitted_trk2_eta"),cand.userFloat("fitted_trk2_phi"),trk.eta(), trk.phi());
	  float dr_to_b   = deltaR(cand.userFloat("fitted_eta"),     cand.userFloat("fitted_phi"),     trk.eta(), trk.phi());
	  
	  if (dr_to_l1 < 0.4){
	    l1_iso04 += trk.pt();
	    if ( dr_to_l1 < 0.3) l1_iso03 += trk.pt();
	  }
	  if (dr_to_l2 < 0.4){
	    l2_iso04 += trk.pt();
	    if (dr_to_l2 < 0.3)  l2_iso03 += trk.pt();
	  }
	  if (dr_to_tk1 < 0.4){
	    tk1_iso04 += trk.pt();
	    if (dr_to_tk1 < 0.3) tk1_iso03 += trk.pt();
	  }
	  if (dr_to_tk2 < 0.4){
	    tk2_iso04 += trk.pt();
	    if (dr_to_tk2 < 0.3) tk2_iso03 += trk.pt();
	  }
	  if (dr_to_b < 0.4){
	    b_iso04 += trk.pt();
	    if (dr_to_b < 0.3) b_iso03 += trk.pt();
	  }
	}
	cand.addUserFloat("l1_iso03" , l1_iso03);
	cand.addUserFloat("l1_iso04" , l1_iso04);
	cand.addUserFloat("l2_iso03" , l2_iso03);
	cand.addUserFloat("l2_iso04" , l2_iso04);
	cand.addUserFloat("tk1_iso03", tk1_iso03);
	cand.addUserFloat("tk1_iso04", tk1_iso04);
	cand.addUserFloat("tk2_iso03", tk2_iso03);
	cand.addUserFloat("tk2_iso04", tk2_iso04);
	cand.addUserFloat("b_iso03"  , b_iso03 );
	cand.addUserFloat("b_iso04"  , b_iso04 );
	*/

	ret_val->push_back(cand);

      } // pi-pi Loop
    }   // mu-mu Loop
  }     // K0s Loop

  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToK0sMuMuPiPiBuilder);
