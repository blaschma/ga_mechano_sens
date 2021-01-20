import os
import fnmatch
import os.path
from os import path
import numpy as np
import sys
from turbomoleOutputProcessing import turbomoleOutputProcessing as top

def load_transmission_data(gen_dir):
	"""
	Loads transmission data of all individuals in given dir. 

	Args:
		param1 (String): dir to process

	Returns:
		list(np.ndarray) [individual][0=disp, 1=T_est, stretching step]

	"""	

	dirs = os.listdir(gen_dir)
	dirs = [i for i in dirs if os.path.isdir(gen_dir + "/" + i)]
	dirs = sorted(dirs)
	print(dirs)
	
	T_estimates = list()

	#read stiffness files and append stiffness and std to lists
	for i in range(len(dirs)):
		print("checking " + str(dirs[i]))
		for file in os.listdir(gen_dir + "/" + dirs[i]):
			if fnmatch.fnmatch(file, '*T_estimate.dat'):
				transmission_file = file
		print("transmission file " + str(transmission_file))
		transmission_file = gen_dir + "/" + dirs[i] + "/" + transmission_file
		dat_Content = top.read_plot_data(transmission_file)[0]
		print(dat_Content[1,:])
		T_estimates.append(dat_Content)

	return T_estimates



def load_stiffness_data(gen_dir):
	"""
	Loads stiffness of all individuals in given dir. 

	Args:
		param1 (String): dir to process

	Returns:
		list,list: (stiffness, std_stiffness)

	"""
	#load list of dirs and ensure it is a dir. Sort them 
	print("load fitness data")
	dirs = os.listdir(gen_dir)
	dirs = [i for i in dirs if os.path.isdir(gen_dir + "/" + i)]
	dirs = sorted(dirs)
	print(dirs)
	stiffness = list()
	std_stiffness = list()

	#read stiffness files and append stiffness and std to lists
	for i in range(len(dirs)):
		try:
			print("open " + str(dirs[i]))
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
		print(stiffness_line)
		stiffness.append(float(stiffness_line[0]))
		std_stiffness.append(float(stiffness_line[3]))


	return stiffness, std_stiffness

def eval_fittness(stiffness, std_stiffness, T_est):
	"""
	Evaluates stiffness values. 

	Args:
		param1 (List): stiffness values (force constants)
		param2 (List): std stiffness values
	Returns:
		np.array: (fitness value)

	"""
	print("eval fitness")
	#print(stiffness)
	mean = 0
	counter = 0
	for i in range(len(stiffness)):
		#negative stiffness is not possible -> set values to low fittness
		print(stiffness[i])

		if(stiffness[i]<0.0):
			print("negative!")
			stiffness[i] = 100.0
		elif(np.abs(std_stiffness[i]/stiffness[i])>0.25):
			print("std to big!")
			stiffness[i] = 100.0
		else:
			mean+=stiffness[i]
			counter+=1.
	mean = mean/counter

	print(stiffness)
	stiffness = np.asarray(stiffness)
	fittness = 1/(stiffness+mean)
	print(fittness)

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
	#stiffness, std_stiffness = load_stiffness_data(path)
	T_est = load_transmission_data(path)
	print(T_est)
	#fittness = eval_fittness(stiffness, std_stiffness)
	#write_fittness(fittness,path)
	
