#ifndef PhysicsTools_XNano_KinVtxFitterWithMassConstraint
#define PhysicsTools_XNano_KinVtxFitterWithMassConstraint

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include <vector>

class KinVtxFitterWithMassConstraint {
public: 
  KinVtxFitterWithMassConstraint():
    fitted_vtx_{}, 
    fitted_state_{},
    fitted_particle_{} {};
      
  KinVtxFitterWithMassConstraint(float particlemass, float particlesigma, RefCountedKinematicTree vertexFitTree);  
  
  ~KinVtxFitterWithMassConstraint() {};
  
  bool success() const {return success_;}
  float chi2() const {return success_ ? fitted_vtx_->chiSquared() : 999;}
  float dof() const  {return success_ ? fitted_vtx_->degreesOfFreedom() : -1;}
  float prob() const {
    return success_ ? ChiSquaredProbability(chi2(), dof()) : 0.;
  }
  float kin_chi2() const {return kin_chi2_;} 
  float kin_ndof() const {return kin_ndof_;}
    
  const KinematicState fitted_candidate() const {
    return fitted_state_;
  }

  const RefCountedKinematicVertex fitted_refvtx() const {
    return fitted_vtx_;
  }

  const RefCountedKinematicParticle fitted_particle() const {     
    return fitted_particle_;
  }

  const math::PtEtaPhiMLorentzVector fitted_p4() const { 
    return math::PtEtaPhiMLorentzVector(
					fitted_state_.globalMomentum().perp(), 
					fitted_state_.globalMomentum().eta() ,
					fitted_state_.globalMomentum().phi() ,
					fitted_state_.mass()
					);
  }
  
  const GlobalPoint fitted_pos() const { 
    return GlobalPoint(
		       fitted_state_.globalMomentum().x(), 
		       fitted_state_.globalMomentum().y() ,
		       fitted_state_.globalMomentum().z() 
		       );
  }

  GlobalPoint fitted_vtx() const {
    return fitted_vtx_->position();
  }
  
  GlobalError fitted_vtx_uncertainty() const {
    return fitted_vtx_->error();
  }
  
private:
  float kin_chi2_ = 0.;
  float kin_ndof_ = 0.;
  bool success_ = false;

  RefCountedKinematicVertex fitted_vtx_; 
  KinematicState fitted_state_;
  RefCountedKinematicParticle fitted_particle_;
};
#endif
