// Merges the PFPackedCandidates and Lost tracks
// beam spot readout in case dcasig to be calculated wrt beam spot
// currently computed wrt triggeringMuon vertex

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/AssociationVector.h"

#include "helper.h"

class TrackMerger : public edm::global::EDProducer<> {


public:
  
  explicit TrackMerger(const edm::ParameterSet &cfg):
    beamSpotSrc_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
    tracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    lostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    muonToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
    vertexToken_(consumes<reco::VertexCollection> (cfg.getParameter<edm::InputTag>( "vertices" ))), 
    triggerBits_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("bits"))),  
    triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(cfg.getParameter<edm::InputTag>("objects"))), 
    HLTPaths_(cfg.getParameter<std::vector<std::string>>("HLTPaths")),      
    drForTriggerMatch_(cfg.getParameter<double>("drForTriggerMatch")),
    trkPtCut_(cfg.getParameter<double>("trkPtCut")),
    trkEtaCut_(cfg.getParameter<double>("trkEtaCut")),
    dcaSig_(cfg.getParameter<double>("dcaSig")),
    trkNormChiMin_(cfg.getParameter<int>("trkNormChiMin")),
    trkNormChiMax_(cfg.getParameter<int>("trkNormChiMax")) 
  {
    produces<pat::CompositeCandidateCollection>("SelectedTracks");  
    produces<TransientTrackCollection>("SelectedTransientTracks");  
  }
  
  ~TrackMerger() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;  
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;  

  // for trigger match 
  std::vector<std::string> HLTPaths_;  
  const double drForTriggerMatch_; 

  // selections                                                                 
  const double trkPtCut_;
  const double trkEtaCut_;
  const double dcaSig_;
  const int trkNormChiMin_;
  const int trkNormChiMax_;
};


void TrackMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &stp) const {

  // input
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  evt.getByToken(beamSpotSrc_, beamSpotHandle);
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("BToKstllProducer") << "No beam spot available from Event" ;
  }  
  const reco::BeamSpot& beamSpot = *beamSpotHandle;
  
  edm::ESHandle<MagneticField> bFieldHandle;
  stp.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  edm::Handle<pat::PackedCandidateCollection> tracks;
  evt.getByToken(tracksToken_, tracks);
  edm::Handle<pat::PackedCandidateCollection> lostTracks;
  evt.getByToken(lostTracksToken_, lostTracks);

  edm::Handle<pat::MuonCollection> muons;
  evt.getByToken(muonToken_, muons);

  edm::Handle<reco::VertexCollection> vertexHandle;
  evt.getByToken(vertexToken_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();

  edm::Handle<edm::TriggerResults> triggerBits;
  evt.getByToken(triggerBits_, triggerBits);
  
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  evt.getByToken(triggerObjects_, triggerObjects);

  // for lost tracks / pf discrimination
  unsigned int nTracks = tracks->size();
  unsigned int totalTracks = nTracks + lostTracks->size();
  
  // Outputs
  std::unique_ptr<pat::CompositeCandidateCollection> tracks_out      (new pat::CompositeCandidateCollection);
  std::unique_ptr<TransientTrackCollection>          trans_tracks_out(new TransientTrackCollection);

  std::vector< std::pair<pat::CompositeCandidate,reco::TransientTrack> > vectrk_ttrk; 


  // vectors for trigger
  std::vector<std::vector<int>> fires;
  std::vector<std::vector<float>> matcher; 
  std::vector<std::vector<float>> DR;

  
  // Loop over tracks and apply preselection
  std::vector<pat::PackedCandidate> preselTracks;

  for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
    const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*tracks)[iTrk] : (*lostTracks)[iTrk-nTracks];
    
    // arranging cuts for speed
    if (!trk.hasTrackDetails())  continue;    
    if (abs(trk.pdgId()) != 211) continue; 
    if (trk.pt() < trkPtCut_ )   continue;
    if (fabs(trk.eta()) > trkEtaCut_) continue;

    if ( (trk.bestTrack()->normalizedChi2() < trkNormChiMin_ &&
	  trkNormChiMin_>=0 ) ||
	 (trk.bestTrack()->normalizedChi2() > trkNormChiMax_ &&
	  trkNormChiMax_>0)  )    continue; 

    // high purity requirment applied only in packedCands
    if( iTrk < nTracks && !trk.trackHighPurity()) continue;

    // distance closest approach in x,y wrt beam spot
    const reco::TransientTrack trackTT( (*trk.bestTrack()) , &(*bFieldHandle));
    std::pair<double,double> DCA = computeDCA(trackTT, beamSpot);
    float DCABS = DCA.first;
    float DCABSErr = DCA.second;
    float DCASig = (DCABSErr != 0 && float(DCABSErr) == DCABSErr) ? fabs(DCABS/DCABSErr) : -1;
    if (DCASig >  dcaSig_  && dcaSig_ >0) continue;

    // Selected tracks
    preselTracks.push_back(trk);  
  }
  unsigned int numPresTracks = preselTracks.size();

  // std::cout << "totalTracks = " << totalTracks << ", presel = " << numPresTracks << std::endl;



  // First do trigger match, only for selected tracks
  for( unsigned int iTrk=0; iTrk<numPresTracks; ++iTrk ) {
    pat::PackedCandidate trk = preselTracks[iTrk];
    
    // These vectors have one entry per HLT path
    std::vector<int> frs(HLTPaths_.size(),0);              
    std::vector<float> temp_matched_to(HLTPaths_.size(),1000.);
    std::vector<float> temp_DR(HLTPaths_.size(),1000.);

    // Loop over trigger paths
    int ipath=-1;
    for (const std::string path: HLTPaths_){
      
      ipath++;
      
      // Here we loop over trigger objects
      float minDr = 1000.;
      float minPt = 1000.;
      
      const edm::TriggerNames &triggerNames = evt.triggerNames(*triggerBits);
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
	obj.unpackPathNames(triggerNames);

	// Is this HLT object matching the track? 
	bool hltMatchOffline = false;
	TVector3 p3triggerObj;
	p3triggerObj.SetPtEtaPhi(obj.pt(), obj.eta(), obj.phi());  
	TVector3 p3Track;                                           
	p3Track.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());    
	float dr = p3triggerObj.DeltaR(p3Track);
	if (dr<drForTriggerMatch_) hltMatchOffline = true;

	// Is this HLT object matching this path?
	bool hltMatchThisPath = false;
	if (hltMatchOffline==true) { 
	  char cstr[ (path+"*").size() + 1]; 
	  strcpy( cstr, (path+"*").c_str() ); 
	  bool isBoth = obj.hasPathName( cstr, true, true );
	  if (isBoth) hltMatchThisPath = true;
	}

	// Look for minDR between reco track and matched HLT object for this path
	if (hltMatchThisPath==true && hltMatchOffline==true) {
	  frs[ipath]=1;
	  if (dr<minDr) {
	    minDr = dr;
	    minPt = obj.pt();
	  }
	}
	
      } // Loop over trigger objects
    
      // Here we store the minimum between reco track and all its matched HLT objects for this HLT path
      temp_DR[ipath]=minDr;
      temp_matched_to[ipath]=minPt;
      
    } // Loop over paths
    
    // One vector per track. Each vector : one element per path (corresponding to the closest HLT object
    fires.push_back(frs);                 // This is used in order to see if a reco track fired a Trigger (1) or not (0).
    matcher.push_back(temp_matched_to);   // PT of the reco track matched to HLT object
    DR.push_back(temp_DR);
    
  } // Loop over reco track


  // Now check for different reco tracks that are matched to the same HLTObject.
  for(unsigned int path=0; path<HLTPaths_.size(); path++){

    for( unsigned int iTrk=0; iTrk<numPresTracks; iTrk++ ) {
      for(unsigned int itr=(iTrk+1); itr<numPresTracks; itr++){
	if(matcher[iTrk][path]!=1000. && matcher[iTrk][path]==matcher[itr][path]){
	  if(DR[iTrk][path]<DR[itr][path]){     // Keep the one that has the minimum DR with the HLT object
	    fires[itr][path]=0;
	    matcher[itr][path]=1000.;
	    DR[itr][path]=1000.;                       
	  }
	  else{
	    fires[iTrk][path]=0;
	    matcher[iTrk][path]=1000.;
	    DR[iTrk][path]=1000.;                       
	  }
	}              
      }
    }
  }

  
  // Loop over tracks and save all tracks passing the selection
  for( unsigned int iTrk=0; iTrk<numPresTracks; ++iTrk ) {

    pat::PackedCandidate trk = preselTracks[iTrk];
    const reco::TransientTrack trackTT( (*trk.bestTrack()) , &(*bFieldHandle));

    // clean tracks wrt to all muons
    int matchedToMuon       = 0;
    int matchedToLooseMuon  = 0;
    int matchedToSoftMuon   = 0;
    for (const pat::Muon &imutmp : *muons) {
      for (unsigned int i = 0; i < imutmp.numberOfSourceCandidatePtrs(); ++i) {
	if (! ((imutmp.sourceCandidatePtr(i)).isNonnull() && 
	       (imutmp.sourceCandidatePtr(i)).isAvailable())
	    )   continue;
	
	const edm::Ptr<reco::Candidate> & source = imutmp.sourceCandidatePtr(i);
	if (source.id() == tracks.id() && source.key() == iTrk){
	  matchedToMuon =1;
	  if (imutmp.isLooseMuon())    matchedToLooseMuon  = 1;
	  if (imutmp.isSoftMuon(PV))   matchedToSoftMuon   = 1;
	  break;
	}
      }
    }

    pat::CompositeCandidate pcand;
    pcand.setP4(trk.p4());
    pcand.setCharge(trk.charge());
    pcand.setVertex(trk.vertex());
    pcand.setPdgId(trk.pdgId());
    pcand.addUserFloat("dxy", trk.dxy());
    pcand.addUserFloat("dxyS", trk.dxy()/trk.dxyError());
    pcand.addUserFloat("dz", trk.dz()); 
    pcand.addUserFloat("dzS", trk.dz()/trk.dzError());
    pcand.addUserInt("isMatchedToMuon", matchedToMuon);
    pcand.addUserInt("isMatchedToLooseMuon", matchedToLooseMuon);
    pcand.addUserInt("isMatchedToSoftMuon",  matchedToSoftMuon);
    pcand.addUserInt("nValidHits", trk.bestTrack()->found());
    // trigger match
    for(unsigned int i=0; i<HLTPaths_.size(); i++){
      pcand.addUserInt(HLTPaths_[i],fires[iTrk][i]);
    }

    // in order to avoid revoking the sxpensive ttrack builder many times and still have everything sorted, we add them to vector of pairs
    vectrk_ttrk.emplace_back( std::make_pair(pcand,trackTT ) );   
  }
  
  // sort to be uniform with leptons
  std::sort( vectrk_ttrk.begin(), vectrk_ttrk.end(), 
             [] ( auto & trk1, auto & trk2) -> 
	     bool {return (trk1.first).pt() > (trk2.first).pt();} 
	     );
  
  // finnaly save ttrks and trks to the correct _out vectors
  for ( auto & trk: vectrk_ttrk){
    tracks_out -> emplace_back( trk.first);
    trans_tracks_out -> emplace_back(trk.second);
  }
  
  evt.put(std::move(tracks_out),       "SelectedTracks");
  evt.put(std::move(trans_tracks_out), "SelectedTransientTracks");
}


//define this as a plug-in
DEFINE_FWK_MODULE(TrackMerger);
