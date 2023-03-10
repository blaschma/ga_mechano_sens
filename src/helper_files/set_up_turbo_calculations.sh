#this script sets up the turbomole calculations.
# $1: path to dir where turbo calc should be initialized. Must be prepared with coord and limits -> genome_to_molecule.py
# $2: path to config file
# $3: generation
# $4: individual

origin=$(pwd)
#load config file
config_file=$2
source $config_file

#set path
path=$1

#set cpus per task
cpus_per_task=$(head -n 1 $path/complexity)
echo $cpus_per_task

# set Intel environment
module load intel
module load openmpi
module load mkl

#set turbo path
source $turbopath


cd $path

#make dirs
mkdir disp_pos
cd disp_pos
#cp $helper_files/run_disp.sh ./
cp $helper_files/jobgen.py ./jobgen.py
cp $helper_files/run_disp.sh ./run_disp.sh
mkdir 0000
cd 0000
cp ../../coord ./
cp ../../limits ./
define < $helper_files/$define_input > define.out
cd ..

if [[ "$queuing" == "SLURM" ]]; then
    echo "Using Slurm"
    sbatch --job-name=gen$3id$4p --mem-per-cpu=$mem_per_cpu --partition=$partition --time=$max_time --ntasks=1 --cpus-per-task=$cpus_per_task --signal=B:SIGUSR1@$kill_time run_disp.sh 0 $num_stretching_steps_pos $displacement_per_step $config_file
elif [[ "$queuing" == "GE" ]]; then
    echo "Using Grid engine (GE)"
    qsub -N gen$3id$4p -cwd -q scc -pe openmp $cpus_per_task -l h_vmem=$mem_per_cpu run_disp.sh 0 $num_stretching_steps_pos $displacement_per_step $config_file
else
    echo "Unknwon queuing system"
    return -1
fi

cd ..

mkdir disp_neg
cd disp_neg
#cp $helper_files/run_disp.sh ./
cp $helper_files/jobgen.py ./jobgen.py
cp $helper_files/run_disp.sh ./run_disp.sh
mkdir 0000
cd 0000
cp ../../coord ./
cp ../../limits ./
define < $helper_files/build_calc > define.out
cd ..
if [[ "$queuing" == "SLURM" ]]; then
    echo "Using Slurm"
    sbatch --job-name=gen$3id$4n --mem-per-cpu=$mem_per_cpu --partition=$partition --time=$max_time --ntasks=1 --cpus-per-task=$cpus_per_task --signal=B:SIGUSR1@$kill_time run_disp.sh 0 $num_stretching_steps_neg -$displacement_per_step $config_file
elif [[ "$queuing" == "GE" ]]; then
    echo "Using Grid engine (GE)"
    qsub -N gen$3id$4n -cwd -q scc -pe openmp $cpus_per_task -l h_vmem=$mem_per_cpu run_disp.sh 0 $num_stretching_steps_neg -$displacement_per_step $config_file
else
    echo "Unknwon queuing system"
    return -1
fi
cd ..

cd $origin

