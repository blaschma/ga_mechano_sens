#!/bin/bash



dispdir=$(pwd) # run this script in the directory, where lastdir is located (so either in disp_pos or disp_neg)
lastdir=$1 # provide path to starting point WITHOUT LEADING ZEROS! (e.g. 0003 -> 3)
stopdir=$2 # defines number of displacements to make (=stopdir-lastdir), WITHOUT LEADING ZEROS
displacement=$3 # Angstrom
config_file=$4 # path to config file
source $config_file

#catch signal
trap  ". $helper_files/action_before_death.sh $dispdir $config_file"  SIGUSR1
#trap  "touch ./signal"  SIGUSR1

# set Intel environment
module load intel
module load openmpi
module load mkl

#set turbo path
. $turbopath







export PARA_ARCH=SMP

#export SMPCPUS=$cpus_per_task
cpus_per_task=$(head -n 1 ../complexity)
export SMPCPUS=$cpus_per_task

ulimit -s unlimited
ulimit -a > mylimits.out



jobgen=$helper_files/jobgen.py

lastdirfp=$dispdir/$(printf "%04d" $lastdir)
if [ ! -d "$lastdirfp" ]; then
    echo "Path" $lastdirfp "not found!"
    exit 1
fi

echo "starting from" $lastdirfp
echo "displacement:" $displacement "Angstrom"

# read lead limits, prepare variables, relax structures in 0000, cd to 0000
cd $(printf "%04d" $lastdir)

#check if prerelax (unfix sulfur, relax, align along z again and fix sulfur again) is wanted. 
if [ "$prerelax" == "T" ]; then
    echo "Prerelaxation...."
    define < $helper_files/build_calc > define.out
    #save fixed atoms, remove fix command from coord file    
    awk '{print $5}' coord > fixed
    grep -irn "f" fixed | cut -f1 -d: > fixed_lines
    awk '{print $1,$2,$3,$4}' coord > coord_unfixed
    rm -r coord 
    cp coord_unfixed coord
    rm -r coord_unfixed
    #relax
    jobex -c $relax_iterations -level $relax_level > jobex_prerelax.log
    #fix again
    paste coord fixed > coord_fixed_again
    rm -r coord
    cp coord_fixed_again coord
    rm -r coord_fixed_again
    #align anchors along z again and update limits
    python3 $helper_files/align_anchor_update_limits.py $dispdir/$(printf "%04d" $lastdir) $config_file        
   
fi
jobex -c $relax_iterations -level $relax_level > jobex.log
file=GEO_OPT_FAILED
if test -f "$file" ; then 
    echo "Geo opt failed"
    file=../../KILL_SIGNAL_SET
    if test -f "$file" ; then 
        echo "Kill signal set"
    else
    . $helper_files/action_before_death.sh $dispdir $config_file
    exit 0
    fi
fi
file=GEO_OPT_CONVERGED 
if test -f "$file" ; then  
    ls
    echo "geo opt converged... starting dft"
    t2x coord > coord.xyz
    ridft -> ridft.out
fi

lower=$(sed '1q;d' limits )
upper=$(sed '2q;d' limits )
cd $dispdir
disphalf=$(echo $displacement/2|bc -l)
currentdir=$lastdir

echo $displacement

while [ $currentdir -lt $stopdir ]
do
    file=../KILL_SIGNAL_SET
    if test -f "$file"; then        
        exit 0
    fi

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
    file=GEO_OPT_FAILED
    if test -f "$file" ; then 
        echo "Geo opt failed"
        file=../../KILL_SIGNAL_SET
        if test -f "$file" ; then 
            echo "Kill signal set"
        else
        . $helper_files/action_before_death.sh $dispdir $config_file
        exit 0
        fi
    fi
    file=GEO_OPT_CONVERGED 
    if test -f "$file" ; then  
    	ls
    	echo "geo opt converged... starting dft"
    	ridft -> ridft.out
    fi
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
        file=../KILL_SIGNAL_SET
        if test -f "$file"; then 
            echo "kill signal has been set"
        else
            echo "now evaluate everything"
            #all the evaluation scripts

            #stiffness evaluation
            GRANDDADDY="$(cd ../; pwd)"
            . $helper_files/eval_stiffness.sh $GRANDDADDY $filename $config_file

            #T estimation
            #take number of occupied orbital from 0000/ridft.out
            line=$(grep -o " number of occupied orbitals : .*" ./0000/ridft.out)
            homo=${line//$"number of occupied orbitals :"}
            python3 $helper_files/eval_propagator.py ../ $filename $homo $config_file


            #now everything is done
            echo "i reached the end...tell everyone"
            touch ../../${filename}_DONE
        fi
    fi
fi

#check if calculations of all individuals are ready but only if kill signal has not been set
num_finished=$(ls ../../ -1f | grep _DONE | wc -l) 
if [ "$num_finished" -eq "$population_size" ]; then
    file=../KILL_SIGNAL_SET
    if test -f "$file"; then 
        echo "kill signal has been set"
    else
        echo "Everybody seems to be ready"
        #eval fitness
        python3 $genetic_algorithm_path/scr/helper_files/eval_fitness.py $calculation_path"/generation_data/"$(basename ${PWD%/*/*}) $config_file

        #invoke next generation
        python3 $genetic_algorithm_path/scr/genetic/invoke_next_generation.py $config_file $calculation_path

        #plot fitness values
        python3 $genetic_algorithm_path/scr/helper_files/plot_fitness_values.py $calculation_path $config_file $calculation_path
    fi
fi


