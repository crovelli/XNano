# Code for X(3872)-related analyses

## Getting started

```shell
cmsrel CMSSW_10_6_29
cd CMSSW_10_6_29/src
cmsenv
git cms-init
```


## Add the modification needed to use the KinematicParticleVertexFitter 

```shell
git remote add crovelli-cmssw https://github.com/crovelli/cmssw.git -f -t forX-from-CMSSW_10_6_29
cp ~/public/X/sparse-checkout_X-10-6-29 .git/info/sparse-checkout
git checkout -b forX-from-CMSSW_10_6_29 crovelli-cmssw/forX-from-CMSSW_10_6_29
cp ~/public/X/sparse-checkout_clean .git/info/sparse-checkout
```

## Add the XNano package and build everything

```shell
git clone git@github.com:crovelli/XNano.git ./PhysicsTools/XNano
cd PhysicsTools/XNano
git checkout CMSSW_10_6_29_slim_winter22_for2016
cd ../../
scram b
```


## To run on a test file
```shell
cd PhysicsTools/XNano/test/
cmsenv 
cmsRun run_nano_cfg.py
```