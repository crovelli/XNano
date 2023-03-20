# Code for X(3872)-related analyses

## Getting started

```shell
cmsrel CMSSW_12_4_11
cd CMSSW_12_4_11/src
cmsenv
git cms-init
```


## Add the modification needed to use the KinematicParticleVertexFitter 

```shell
git remote add crovelli-cmssw https://github.com/crovelli/cmssw.git -f -t forX-from-CMSSW_12_4_11
cp ~/public/X/sparse-checkout_X-12-4-11 .git/info/sparse-checkout
git checkout -b forX-from-CMSSW_12_4_11 crovelli-cmssw/forX-from-CMSSW_12_4_11
cp ~/public/X/sparse-checkout_clean .git/info/sparse-checkout
```

## Add the XNano package and build everything

```shell
git clone git@github.com:crovelli/XNano.git ./PhysicsTools/XNano
cd PhysicsTools/XNano
git checkout CMSSW_12_4_11_slim_winter22
cd ../../
scram b
```


## To run on a test file
```shell
cd PhysicsTools/XNano/test/
cmsenv 
cmsRun run_nano_cfg.py
```