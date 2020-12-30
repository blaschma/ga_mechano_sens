#!/bin/bash


dispdir=$(pwd) # run this script in the directory, where lastdir is located (so either in disp_pos or disp_neg)
lastdir=$1 # provide path to starting point WITHOUT LEADING ZEROS! (e.g. 0003 -> 3)
stopdir=$2 # defines number of displacements to make (=stopdir-lastdir), WITHOUT LEADING ZEROS
displacement=$3 # Angstrom
config_file=$4 # path to config file
source $config_file


#set turbo path
. $turbopath

# set Intel environment
INTEL=$intel_path
. $INTEL/bin/compilervars.sh intel64




export PARA_ARCH=SMP
export SMPCPUS=$cpus_per_task




jobgen=$helper_files/jobgen.py

lastdirfp=$dispdir/$(printf "%04d" $lastdir)
if [ ! -d "$lastdirfp" ]; then
    echo "Path" $lastdirfp "not found!"
    exit 1
fi

echo "starting from" $lastdirfp
echo "displacement:" $displacement "Angstrom"

# read lead limits, prepare variables
cd $(printf "%04d" $lastdir)
jobex -c $relax_iterations -level $relax_level > jobex.log
ridft -> ridft.out
lower=$(sed '1q;d' limits )
upper=$(sed '2q;d' limits )
cd $dispdir
disphalf=$(echo $displacement/2|bc -l)
currentdir=$lastdir

echo $displacement

while [ $currentdir -lt $stopdir ]
do
    currentdir=$[$lastdir+1]

    # zero padding
    printf -v zpcurrentdir "%04d" $currentdir
    printf -v zplastdir "%04d" $lastdir

    echo "copying" $zplastdir "to" $zpcurrentdir
    rm -rf $zpcurrentdir
    cp -r $zplastdir $zpcurrentdir
    
    echo "displacing atoms with block limits" $lower "and" $upper
    
    python  $jobgen  -displace $zplastdir/coord $zpcurrentdir/coord $displacement $lower $upper
    lower=$(echo $lower - $disphalf|bc -l)
    upper=$(echo $upper + $disphalf|bc -l)

    # save new limits for further post processing
    cd $zpcurrentdir
    rm limits
    echo $lower >> limits
    echo $upper >> limits

    echo "starting turbomole calculation in" $zpcurrentdir
    jobex -c $relax_iterations -level $relax_level > jobex.log
    ridft -> ridft.out
    cd $dispdir
    lastdir=$currentdir
done

cd $dispdir

#notify the world when calc completed
parent=$(basename ${PWD%})
if [ "$parent" = "disp_pos" ]; then
   echo "This was disp_pos"
   touch ../DISP_POS_DONE
fi
if [ "$parent" = "disp_neg" ]; then
   echo "This was disp_neg"
   touch ../DISP_NEG_DONE
fi

#generation_individual
filename=$(basename ${PWD%/*/*})_$(basename ${PWD%/*})

#if calc for pos and neg displacement are ready -> evaluation
file=../DISP_POS_DONE
if test -f "$file"; then
    file=../DISP_NEG_DONE
    if test -f "$file"; then
        echo "now evaluate everything"
        #all the evaluation scripts
        . $helper_files/eval_stiffness.sh $(basename ${PWD%/*}) $filename $config_file


        #now everything is done
        echo "i reached the end...tell everyone"
        touch ../../${filename}_DONE
    fi
fi



