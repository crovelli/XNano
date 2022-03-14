# Code for X(3872)-related analyses

## Getting started

```shell
cmsrel CMSSW_9_4_7
cd CMSSW_9_4_7/src
cmsenv
git cms-init
```


## Add the modification needed to use the KinematicParticleVertexFitter + tables for nanoAODs

```shell
git remote add crovelli-cmssw https://github.com/crovelli/cmssw.git -f -t forX-from-CMSSW_9_4_7
cp ~/public/X/sparse-checkout_X-9-4-7 .git/info/sparse-checkout
git checkout -b forX-from-CMSSW_9_4_7 crovelli-cmssw/forX-from-CMSSW_9_4_7
cp ~/public/X/sparse-checkout_clean .git/info/sparse-checkout
```

## Add the XNano package and build everything

```shell
git clone git@github.com:crovelli/XNano.git ./PhysicsTools/XNano
PhysicsTools/XNano
git checkout CMSSW_9_4_7
scram b
```


## To run on a test file
```shell
cd PhysicsTools/XNano/test/
cmsenv 
cmsRun run_nano_cfg.py
```