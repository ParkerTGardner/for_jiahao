#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters"
    echo "<Inputs File> <nJobs>"
    exit
fi

#pathvar=sept_out18


inDir=$1
nJobs=$2

#2d corr runs
#hostname > 'src/comp/'$pathvar'/host.log'
#date +"%m-%d-%y"     >> 'src/comp/'$pathvar'/host.log'
#date +"%T"     >> 'src/comp/'$pathvar'/host.log'
#echo "$inDir" >> 'src/comp/'$pathvar'/host.log'
#egrep "src/comp" src/MyClass.C | grep -v "//" >> 'src/comp/'$pathvar'/host.log'
#egrep "const int   trackbinbounds" src/MyClass.C | grep -v "//" >> 'src/comp/'$pathvar'/host.log'
#egrep "const float ptbinbounds" src/MyClass.C | grep -v "//" >> 'src/comp/'$pathvar'/host.log'


#hostname > 'src/'$pathvar'/host.log'
#date +"%m-%d-%y"     >> 'src/'$pathvar'/host.log'
#date +"%T"     >> 'src/'$pathvar'/host.log'
#echo "$inDir" >> 'src/'$pathvar'/host.log'
#egrep "src/major_dist" src/MyClass.C | grep -v "//" >> 'src/'$pathvar'/host.log'
#egrep "const int   trackbinbounds" src/MyClass.C | grep -v "//" >> 'src/'$pathvar'/host.log'
#egrep "const float ptbinbounds" src/MyClass.C | grep -v "//" >> 'src/'$pathvar'/host.log'


for i in $(seq 0 $(($nJobs-1)) );
do 
    ./bin/MyClass.exe $inDir $i $nJobs > 'wide_eta_'$i'.log' 2>&1 & 
    #disown -h  
done
