#this script collects the total energy and saves it to file

# $1 directory to process
# $2 filename
# $3 config file 

config_file=$3 # path to config file
source $config_file

# set Intel environment
INTEL=$intel_path
. $INTEL/bin/compilervars.sh intel64


exec 3<> ../$2_totalEnergy.dat

#collects total energy in given directory ($1)

start_dir=$(pwd)
cd $1

#process disp_pos


echo "processing disp_pos"
cd $start_dir/$1
cd "disp_pos"
old_dir=$(pwd)

for dir in ./*/     # list directories in the form "/tmp/dirname/"
do
    dir=${dir%*/}      # remove the trailing "/"
    echo ${dir##*/}    # print everything after the final "/"
    cd $dir

    converged=GEO_OPT_CONVERGED
    if test -f "$converged"; then
	#echo "Geo Opt converged"
	FILE=ridft.out
	if test -f "$FILE"; then
		#echo "$FILE exists."
			#find total energy in ridft.out
			line=$(grep -o "total energy      =.*" ridft.out)
			#echo $line
			totalEnergy=${line//$"total energy      ="}	
			#echo $totalEnergy
			totalEnergy=${totalEnergy//$"|"}
			#echo $totalEnergy
			#echo $dir
			displacement=${dir//$"./"}
			echo "$displacement	$totalEnergy" >&3

						
					
	fi
   else
	echo "Geo Opt not converged"
		
   fi

    cd $old_dir
done

echo "processing disp_neg"
cd $start_dir/$1
cd "disp_neg"
old_dir=$(pwd)


for dir in ./*/     # list directories in the form "/tmp/dirname/"
do
    dir=${dir%*/}      # remove the trailing "/"
    echo ${dir##*/}    # print everything after the final "/"
    cd $dir

    converged=GEO_OPT_CONVERGED
    if test -f "$converged"; then
	#echo "Geo Opt converged"
	FILE=ridft.out
	if test -f "$FILE"; then
		#echo "$FILE exists."
			#find total energy in ridft.out
			line=$(grep -o "total energy      =.*" ridft.out)
			#echo $line
			totalEnergy=${line//$"total energy      ="}	
			#echo $totalEnergy
			totalEnergy=${totalEnergy//$"|"}
			#echo $totalEnergy
			#echo $dir
			displacement=${dir//$"./"}
			echo "-$displacement	$totalEnergy" >&3		
	fi
   else
	echo "Geo Opt not converged"
		
   fi

    cd $old_dir
done


cd ..

cd $start_dir
num_atoms=$(wc -l < ../coord)

python3 $helper_files/plot_totalEnergy_single.py ../$2 $num_atoms

