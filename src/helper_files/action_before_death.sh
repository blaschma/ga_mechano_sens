#!/bin/bash
origin_path=$(pwd)
echo "action_before_death"
dispdir=$1
config_file=$2
echo $config_file
source $config_file




cd $dispdir

#set kill notification
touch ../KILL_SIGNAL_SET
cd $dispdir

#notify the world when calc completed
parent=$(basename ${PWD%})

touch ../DISP_POS_DONE
touch ../DISP_NEG_DONE


#generation_individual
filename=$(basename ${PWD%/*/*})_$(basename ${PWD%/*})
file=../DISP_POS_DONE
if test -f "$file"; then
    file=../DISP_NEG_DONE
    if test -f "$file"; then
        #copy error files
        parentdir="$(dirname "$dispdir")"
        cp -a $helper_files/error_files/. $parentdir
        touch ../../${filename}_DONE
    fi
fi

#check if calculations of all individuals are ready
num_finished=$(ls ../../ -1f | grep _DONE | wc -l) 
if [ "$num_finished" -eq "$population_size" ]; then
    echo "Everybody seems to be ready"
    #eval fitness
    python3 $genetic_algorithm_path/src/helper_files/eval_fitness.py $calculation_path"/generation_data/"$(basename ${PWD%/*/*}) $config_file

    #invoke next generation
    python3 $genetic_algorithm_path/src/genetic/invoke_next_generation.py $config_file $calculation_path
fi

cd $origin_path