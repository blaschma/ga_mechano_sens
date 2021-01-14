import os
import os.path
from os import path
import numpy as np
import sys

def load_stiffness_data(gen_dir):
	"""
	Loads stiffness of all individuals in given dir. 

	Args:
		param1 (String): dir to process

	Returns:
		list,list: (stiffness, std_stiffness)

	"""
	#load list of dirs and ensure it is a dir. Sort them 
	dirs = os.listdir(gen_dir)
	print(dirs)
	dirs = [i for i in dirs if os.path.isdir(gen_dir + "/" + i)]
	dirs = sorted(dirs)
	
	stiffness = list()
	std_stiffness = list()

	#read stiffness files and append stiffness and std to lists
	for i in range(len(dirs)):
		try:
			stiffness_file = open(gen_dir + "/" + dirs[i] + "/stiffness.dat")
		except OSError as e:
			print("stiffness file not found " + str(e))
			#todo :proper treatment
			stiffness.append(100)
			std_stiffness.append(100)
			continue
			return -1
		#skip first line 
		for line in stiffness_file:
			stiffness_line = line
		stiffness_file.close()

		stiffness_line = stiffness_line.strip().split("	")
		stiffness.append(float(stiffness_line[0]))
		std_stiffness.append(float(stiffness_line[3]))

	print(stiffness)
	print(std_stiffness)

	return stiffness, std_stiffness

def eval_fittness(stiffness, std_stiffness):
	"""
	Evaluates stiffness values. 

	Args:
		param1 (List): stiffness values (force constants)
		param2 (List): std stiffness values
	Returns:
		np.array: (fitness value)

	"""

	print(stiffness)
	for i in range(len(stiffness)):
		#negative stiffness is not possible -> set values to low fittness
		if(stiffness[i]<0.0):
			print("negative!")
			stiffness[i] = 100.0
		if(np.abs(std_stiffness[i]/stiffness[i])>0.05):
			stiffness[i] = 100.0

	print(stiffness)
	stiffness = np.asarray(stiffness)
	fittness = 1/(stiffness+0.2)

	return fittness

def write_fittness(fittness, path):
	"""
	Write fittness values to file

	Args:
		param1 (List): fittness

	Returns:
		

	"""
	file = open(path + "/fitness.dat", "w")
	for i in range(len(fittness)):
		file.write(str(fittness[i])+"\n")
	file.close()








if __name__ == '__main__':
	path = sys.argv[1]
	print(path)
	stiffness, std_stiffness = load_stiffness_data(path)
	fittness = eval_fittness(stiffness, std_stiffness)
	write_fittness(fittness,path)
