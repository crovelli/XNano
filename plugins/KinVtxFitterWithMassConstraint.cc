#include "KinVtxFitterWithMassConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h" 
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"

KinVtxFitterWithMassConstraint::KinVtxFitterWithMassConstraint(
							       float particlemass,
							       float particlesigma,
							       RefCountedKinematicTree vertexFitTree) {
  
  KinematicConstraint * mass_c = new MassKinematicConstraint(particlemass, particlesigma);
  KinematicParticleFitter kcv_fitter;    

  RefCountedKinematicTree vtx_tree = kcv_fitter.fit(mass_c, vertexFitTree);  

  if (vtx_tree->isEmpty() || !vtx_tree->isValid() || !vtx_tree->isConsistent()) {
    success_ = false; 
    return;
  }

  vtx_tree->movePointerToTheTop(); 
  fitted_particle_ = vtx_tree->currentParticle();
  fitted_vtx_ = vtx_tree->currentDecayVertex();
  if (!fitted_particle_->currentState().isValid() || !fitted_vtx_->vertexIsValid()){ 
    success_ = false; 
    return;
  }
  fitted_state_ = fitted_particle_->currentState();

  success_ = true;
}
