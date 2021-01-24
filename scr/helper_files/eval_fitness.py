import os
import fnmatch
import os.path
from os import path
import numpy as np
import sys
from turbomoleOutputProcessing import turbomoleOutputProcessing as top
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

def load_transmission_data(gen_dir):
	"""
	Loads transmission data of all individuals in given dir. 

	Args:
		param1 (String): dir to process

	Returns:
		list(np.ndarray) [individual][0=disp, 1=T_est, stretching step]

	"""	

	dirs = os.listdir(gen_dir)
	dirs = [int(i) for i in dirs if os.path.isdir(gen_dir + "/" + i)]
	dirs = sorted(dirs)
	print(dirs)
	dirs = [str(i) for i in dirs]
	
	T_estimates = list()
	T_estimates_params = list()
	T_estimates_params_list = list()

	#read stiffness files and append stiffness and std to lists
	for i in range(len(dirs)):
		if(i==8):
			#continue
			pass
		print("checking " + str(dirs[i]))
		for file in os.listdir(gen_dir + "/" + dirs[i]):
			if fnmatch.fnmatch(file, '*T_estimate.dat'):
				transmission_file = file
			if fnmatch.fnmatch(file, '*_T_estimate_params.dat'):
				param_file = file
		print("transmission file " + str(transmission_file))
		transmission_file = gen_dir + "/" + dirs[i] + "/" + transmission_file
		dat_Content = top.read_plot_data(transmission_file)[0]
		param_file = open(gen_dir + "/" + dirs[i] + "/" + param_file)
		T_estimates_params = list()
		for line in param_file:
			T_estimates_params.append(float(line))
		#print(dat_Content[1,:])
		T_estimates.append(dat_Content)
		T_estimates_params_list.append(T_estimates_params)

	return T_estimates, T_estimates_params_list



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
	dirs = [int(i) for i in dirs if os.path.isdir(gen_dir + "/" + i)]
	dirs = sorted(dirs)
	print(dirs)
	dirs = [str(i) for i in dirs]
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
	T_est, T_estimates_params_list = load_transmission_data(path)
	"""
	file = open(path + "/T_est_params", "w")
	for i in range(len(T_estimates_params_list)):
		file.write(str(T_estimates_params_list[i][0]).replace(".", ",")+"	"+str(T_estimates_params_list[i][1]).replace(".", ",")+"	"+str(T_estimates_params_list[i][2]).replace(".", ",")+"	"+str(T_estimates_params_list[i][3]).replace(".", ",") + "\n")
	file.close()
	"""
	#plot all T estimates
	#"""
	#print(T_est)
	fig, ax = plt.subplots(1)
	NUM_COLORS = len(T_est)
	cm = plt.get_cmap('tab20')
	ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
	f = open(path + "/T_est_data.dat", "w")
	for i in range(len(T_est)):
		if(i != 6 and i!=7 and i!=8 or True):
			ax.plot(T_est[i][0,:], (T_est[i][1,:]), label=str(i))
			dx = 0.1
			dy = np.diff(T_est[i][1,:])/dx
			print(i)
			
			f.write(str(i)+"\n")
			f.write("median " + str(np.median(T_est[i][1,:]))+"\n")
			f.write("max derivate " + str(np.max(dy))+"\n")
			f.write("median derivate " + str(np.median(dy))+"\n")
			f.write("avg derivate " + str(np.average(dy))+"\n")
			f.write("min derivate " + str(np.min(dy))+"\n")
			f.write("max abs derivate " + str((np.max(np.abs(dy))))+"\n")
			f.write("median abs derivate " + str(np.median(abs(dy)))+"\n")
			f.write("min abs derivate " + str(np.min(abs(dy)))+"\n")

			f.write(".-.-.-.-."+"\n")
	ax.set_yscale('log')
	f.close()
	#ax.set_xlim(-0.1,0.1)
	#ax.set_ylim(-0.0000001,0.00004)
	ax.set_xlabel('Displacement ($\mathrm{\AA}$)',fontsize=20)
	ax.set_ylabel('$\mathrm{T}_{\mathrm{estimate}}$',fontsize=20)
	ax.legend()
	plt.savefig(path + "/T_estimates.pdf", bbox_inches='tight')

	#fittness = eval_fittness(stiffness, std_stiffness)
	#write_fittness(fittness,path)
	#"""
