#!/bin/sh

#source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc7_amd64_gcc700

cd ~/local_not_scp/CMSSW_10_6_4_patch1/src

<<<<<<< HEAD
cd ~/CMSSW_10_6_4_patch1/src

#cmsenv

=======
#cmsenv
>>>>>>> 8db22d7e64236e122b91a2227ad91ec6e4e26177
eval `scramv1 runtime -sh`
cd ~/for_jiahao/HeavyIonAnalysis/TrackAnalysis/batch
echo PWD: $PWD
../bin/PYTHIA_gen.exe ./pythia_lists/list_cor_$1 0 1
