// Merges the PFPackedCandidates and Lost tracks

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "helper.h"
#include "TVectorD.h"    // for fixing tracks 
#include "TMatrixDSym.h" // for fixing tracks 

class TrackMerger : public edm::global::EDProducer<> {


public:
  
  explicit TrackMerger(const edm::ParameterSet &cfg):
    bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),  
    beamSpotSrc_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
    tracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    lostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    muonToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
    vertexToken_(consumes<reco::VertexCollection> (cfg.getParameter<edm::InputTag>( "vertices" ))), 
    trkPtCut_(cfg.getParameter<double>("trkPtCut")),
    trkEtaCut_(cfg.getParameter<double>("trkEtaCut")),
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

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

  // selections                                                                 
  const double trkPtCut_;
  const double trkEtaCut_;
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

  const auto& bField = stp.getData(bFieldToken_);

  edm::Handle<pat::PackedCandidateCollection> tracks;
  evt.getByToken(tracksToken_, tracks);
  edm::Handle<pat::PackedCandidateCollection> lostTracks;
  evt.getByToken(lostTracksToken_, lostTracks);

  edm::Handle<pat::MuonCollection> muons;
  evt.getByToken(muonToken_, muons);

  edm::Handle<reco::VertexCollection> vertexHandle;
  evt.getByToken(vertexToken_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();

  // for lost tracks / pf discrimination
  unsigned int nTracks = tracks->size();
  unsigned int totalTracks = nTracks + lostTracks->size();
  
  // Outputs
  std::unique_ptr<pat::CompositeCandidateCollection> tracks_out      (new pat::CompositeCandidateCollection);
  std::unique_ptr<TransientTrackCollection>          trans_tracks_out(new TransientTrackCollection);

  std::vector< std::pair<pat::CompositeCandidate,reco::TransientTrack> > vectrk_ttrk; 


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

    // high purity requirement applied only in packedCands
    if( iTrk < nTracks && !trk.trackHighPurity()) continue;

    // Selected tracks
    preselTracks.push_back(trk);  
  }
  unsigned int numPresTracks = preselTracks.size();

  
  // Loop over tracks and save all tracks passing the selection
  for( unsigned int iTrk=0; iTrk<numPresTracks; ++iTrk ) {

    pat::PackedCandidate trk = preselTracks[iTrk];
    const reco::TransientTrack trackTT((*(trk.bestTrack())), &bField);    

    // clean tracks wrt to all muons
    int matchedToMuon       = 0;
    int matchedToLooseMuon  = 0;
    int matchedToMediumMuon = 0;
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
	  if (imutmp.isMediumMuon())   matchedToMediumMuon = 1;
	  if (imutmp.isSoftMuon(PV))   matchedToSoftMuon   = 1;
	  break;
	}
      }
    }

    // Was used in Run2 trigger paths, keep in case it's useful for offline
    Basic3DVector<float> thepos( trk.bestTrack()->vertex());
    GlobalPoint thegpos( thepos);
    Basic3DVector<float> themom( trk.bestTrack()->momentum());
    GlobalVector thegmom( themom);
    GlobalTrajectoryParameters thepar( thegpos, thegmom, trk.bestTrack()->charge(), &bField);
    CurvilinearTrajectoryError theerr( trk.bestTrack()->covariance());
    FreeTrajectoryState InitialFTS( thepar, theerr); 
    TSCBLBuilderNoMaterial blsBuilder;
    TrajectoryStateClosestToBeamLine tscb( blsBuilder(InitialFTS, *beamSpotHandle));
    float d0sig=-1000.;
    if (tscb.isValid()) d0sig = tscb.transverseImpactParameter().significance();

    // To be stored for offline selection: min transverse impact parameter wrt all PVs in the event
    float maxD0PV=-10000.;
    float minD0PV= 10000.;
    for (size_t vtx_idx = 0; vtx_idx < vertexHandle->size(); ++vtx_idx) {
      edm::Ptr<reco::Vertex> thisPV(vertexHandle, vtx_idx);
      TrajectoryStateClosestToPoint tscPV;
      tscPV = trackTT.trajectoryStateClosestToPoint(GlobalPoint(thisPV->x(),thisPV->y(),thisPV->z()));
      float d0PV=-1000.;
      if (tscPV.isValid()) d0PV = fabs(tscPV.perigeeParameters().transverseImpactParameter());
      if (d0PV>maxD0PV) maxD0PV=d0PV;
      if (d0PV<minD0PV) minD0PV=d0PV;
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
    pcand.addUserInt("isMatchedToMediumMuon", matchedToMediumMuon);
    pcand.addUserInt("isMatchedToSoftMuon",  matchedToSoftMuon);
    pcand.addUserInt("nValidHits", trk.bestTrack()->found());
    pcand.addUserFloat("d0sig", d0sig);    
    pcand.addUserFloat("maxd0PV", maxD0PV);
    pcand.addUserFloat("mind0PV", minD0PV);

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
