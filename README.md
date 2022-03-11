# Code for X(3872)-related analyses

cmsrel CMSSW_9_4_7
cd CMSSW_9_4_7/src
cmsenv
git cms-init


# Add the modification needed to use the KinematicParticleVertexFitter
# Still to be checked, comment for me moment
# git cms-merge-topic -u CMSBParking:fixKinParticleVtxFitter


# Add the XNano package and build everything
git clone git@github.com:crovelli/XNano.git ./PhysicsTools/XNano
git cms-addpkg PhysicsTools/NanoAOD
scram b

# To run on a test file
cd PhysicsTools/XNano/test/
cmsenv 
cmsRun run_nano_cfg.py