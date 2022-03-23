#ifndef BPHNANO_CLASSES_H
#define BPHNANO_CLASSES_H

#include "DataFormats/Common/interface/Wrapper.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "PhysicsTools/XNano/plugins/KinVtxFitter.h"
#include "PhysicsTools/XNano/plugins/KinVtxFitterWithMassConstraint.h"

#include <vector>

namespace {
  struct dictionary {
      std::vector<reco::TransientTrack> ttv;
      edm::Wrapper<std::vector<reco::TransientTrack> > wttv; 
      edm::Wrapper<std::vector<KinVtxFitter> > wkv;
      edm::Wrapper<std::vector<KinVtxFitterWithMassConstraint> > wkvm;
  };
}

#endif  
