import os
import fnmatch
import os.path
from os import path
import numpy as np
import sys
import tmoutproc as top
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
		for file in os.listdir(gen_dir + "/" + dirs[i]):
			if fnmatch.fnmatch(file, '*T_estimate.dat'):
				transmission_file = file
			if fnmatch.fnmatch(file, '*_T_estimate_params.dat'):
				param_file = file
		print("transmission file " + str(transmission_file))
		transmission_file = gen_dir + "/" + dirs[i] + "/" + transmission_file
		dat_Content = top.read_plot_data(transmission_file)
		param_file = open(gen_dir + "/" + dirs[i] + "/" + param_file)
		T_estimates_params = list()
		for line in param_file:
			#print("line " + line)
			T_estimates_params.append(float(line))
		#print(dat_Content[1,:])
		T_estimates.append(dat_Content)
		T_estimates_params_list.append(T_estimates_params)

	return T_estimates, T_estimates_params_list

def process_T_estimate_data(T_est, gen_path):
	"""
	Processes T estimate data. Plots of all T_estimates are saved 

	Args:
		param1 (List): T estimates
		param2 (String): path to generation data

	Returns:
		

	"""	
	fig, ax = plt.subplots(1)
	NUM_COLORS = len(T_est)
	cm = plt.get_cmap('nipy_spectral')
	ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
	f = open(gen_path + "/T_est_data.dat", "w")
	for i in range(len(T_est)):
		if(i != 6 and i!=7 and i!=8 or True):
			ax.plot(T_est[i][0,:], (T_est[i][1,:]), label=str(i))
			dx = 0.1
			dy = np.diff(T_est[i][1,:])/dx
			#print(i)
			
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
	ax.legend(loc='lower left', ncol = 6, bbox_to_anchor=(0.,1.02,1.,.102), mode="expand", borderaxespad=0., fontsize=10)
	plt.savefig(gen_path + "/T_estimates_summary.pdf", bbox_inches='tight')



def load_stiffness_data(gen_dir):
	"""
	Loads stiffness of all individuals in given dir. 

	Args:
		param1 (String): dir to process

	Returns:
		list,list: (stiffness, std_stiffness)

	"""
	#load list of dirs and ensure it is a dir. Sort them 
	print("load fitness data ")
	print(gen_dir)
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
			for file in os.listdir(gen_dir + "/" + dirs[i]):
				
				if fnmatch.fnmatch(file, '*_stiffness.dat'):
					
					stiffness_file = file
					stiffness_file = open(gen_dir + "/" + dirs[i] + "/" + stiffness_file)
			
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
		#print("stiffness_line")
		#print(float(stiffness_line[0]))
		#print(stiffness_line)
		stiffness.append(float(stiffness_line[0]))
		std_stiffness.append(float(stiffness_line[3]))


	return stiffness, std_stiffness

def eval_fittness(stiffness, std_stiffness, T_est, T_estimates_params):
	"""
	Evaluates stiffness values. 

	Args:
		param1 (List): stiffness values (force constants)
		param2 (List): std stiffness values
		param3 (List): T_estimates
		param4 (List(List)): T_estimates params [individual][0:fit_param, 1: std_fit_param, 2: T_est_min, 3: disp(T_est_min)]
	Returns:
		np.array: (fitness value)

	"""
	#print("eval fitness")
	#evaluate stiffness part
	#mean = 0
	counter = 0
	for i in range(len(stiffness)):

		#negative stiffness is not possible -> set values to low fittness
		print(str(i) + "   " + str(stiffness[i]))

		if(stiffness[i]<0.0):
			print("negative!")
			stiffness[i] = 100.0
		elif(np.abs(std_stiffness[i]/stiffness[i])>0.20):
			print("std to big!")
			stiffness[i] = 100.0
		else:
			#mean+=stiffness[i]
			counter+=1.
	#mean = mean/counter
	#print("stiffness")
	#print(stiffness)
	stiffness = np.asarray(stiffness)
	#print("mean " + str(mean))
	#fittness_stiffness = 1/(stiffness+mean)
	fittness_stiffness = 1/(stiffness+0.005)
	#print("fittness_stiffness")
	#print(fittness_stiffness)

	#evaluate T_estimate part
	min_min_T_est = 5.0
	fit_param = list()
	min_T_est = list()
	median_penalty = list()
	for i in range(len(T_estimates_params)):

		if(np.isfinite(float(T_estimates_params[i][0]))==True and T_estimates_params[i][2] !=0):
			if(np.abs(T_estimates_params[i][1]/(T_estimates_params[i][0]+1E-12)) > 1.0):
				fit_param.append(0)
				median_penalty.append(0)
				#print(np.abs(T_estimates_params[i][1]/(T_estimates_params[i][0]+1E-12)) )
				#print(T_estimates_params)
				#print("else fall oben ")
		
			else:
				if(np.median(T_est[i][1])<1):
					beta=1.2
					gamma=2.0
					nu=1.5
					#print(1/(1+np.exp(-beta*(-np.log(np.median(T_est[i][1])))-gamma)))
					median_penalty.append(1/(1+np.exp(beta*(-np.log(np.median(T_est[i][1])))-gamma)))
					#print("wurzel")
					#print(T_estimates_params[i][0])
					fit_param.append((np.sqrt(T_estimates_params[i][0]))/(np.min(T_est[i][1]))*(((np.median(T_est[i][1])))/(np.min(T_est[i][1])))**nu)
				else:
					fit_param.append(T_estimates_params[i][0])
					median_penalty.append(0)
					#print("else fall")

			min_T_est.append(T_estimates_params[i][2])
			if(T_estimates_params[i][2] < min_min_T_est):
				min_min_T_est = T_estimates_params[i][2]
		else:			
			fit_param.append(0)
			min_T_est.append(1)
			median_penalty.append(0)
		
	fittness_T_est = np.asarray(fit_param)
	median_penalty = np.asarray(median_penalty)
	print("fittness_T_est")
	print(len(fittness_T_est))
	print(fittness_T_est)
	print("median_penalty")
	print(len(median_penalty))
	print(median_penalty)
	if((len(fittness_stiffness) != len(fittness_T_est)) or (len(fittness_stiffness) != len(median_penalty))):
		raise ValueError('Lengths of fitness measures do not match')

	
	return fittness_stiffness,fittness_T_est,median_penalty

def write_fittness(fittness, path):
	"""
	Write fittness values to file

	Args:
		param1 (List): fittness

	Returns:
		

	"""
	file = open(path + "/fitness.dat", "w")
	for i in range(len(fittness[0])):
		print(fittness[2])
		file.write(str(fittness[0][i]*fittness[1][i]*fittness[2][i])+"\n")
	file.close()
	top.write_plot_data(path + "/fittness_contribution", fitness, "stiffness, T_est, median_penalty")







if __name__ == '__main__':
	# sys.argv[1] path to process
	# sys.argv[2] config path

	#"""
	#load all data
	path=sys.argv[1]
	stiffness, std_stiffness = load_stiffness_data(path)
	T_est, T_estimates_params_list = load_transmission_data(path)
	process_T_estimate_data(T_est, path)
	

	#eval fitness
	fitness = eval_fittness(stiffness, std_stiffness, T_est, T_estimates_params_list)
	#print("fitness unten")
	#print(fittness)
	write_fittness(fitness,path)

	"""
	"""
	file = open(path + "/T_est_params", "w")
	for i in range(len(T_estimates_params_list)):
		file.write(str(T_estimates_params_list[i][0]).replace(".", ",")+"	"+str(T_estimates_params_list[i][1]).replace(".", ",")+"	"+str(T_estimates_params_list[i][2]).replace(".", ",")+"	"+str(T_estimates_params_list[i][3]).replace(".", ",") + "\n")
	file.close()
	#"""
	#plot all T estimates
	#"""
	#print(T_est)
	

	#fittness = eval_fittness(stiffness, std_stiffness)
	#
	#"""
