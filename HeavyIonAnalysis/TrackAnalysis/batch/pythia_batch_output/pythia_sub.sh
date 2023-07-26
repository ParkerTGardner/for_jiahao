#!/bin/sh

#source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc7_amd64_gcc700

cd ~/second_checkout


#cmsenv

eval `scramv1 runtime -sh`
cd ~/for_jiahao/HeavyIonAnalysis/TrackAnalysis/batch
echo PWD: $PWD
../bin/PYTHIA_gen.exe  ./pythia_lists/list_cor_$1 0 1
