#collects specified .xyz files from stretching trace
path=$1 #path to dir where disp_pos and disp_neg is located
output_xyz=$2 #file where .xyz trajectory is stored
xyz_name=$3 #names of xyz_file to collect

start_dir=$(pwd)
cd $path
path=$(pwd)

touch $path/$output_xyz
#process disp_neg
cd disp_neg
old_dir=$(pwd)
#store subdirs in array
dirs=( $( ls -1p | grep / | sed 's/^\(.*\)/\1/') )
for ((i=${#dirs[@]}-1; i>=0; i--)); 
do
    cd ${dirs[$i]}
    echo $(pwd)
    t2x coord > tmp_coord.xyz
    cat $xyz_name >> $path/$output_xyz
    cd $old_dir
done    

cd ..
#process disp_pos
cd disp_pos
old_dir=$(pwd)
#store subdirs in array
dirs=( $( ls -1p | grep / | sed 's/^\(.*\)/\1/') )
arraylength=${#dirs[@]}
for (( i=1; i<${arraylength}; i++ ));    # list directories in the form "/tmp/dirname/"
do
    cd ${dirs[$i]}
    echo $(pwd)
    t2x coord > tmp_coord.xyz
    cat $xyz_name >> $path/$output_xyz
    cd $old_dir
done





cd $start_dir

