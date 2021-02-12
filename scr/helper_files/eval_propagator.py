import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from turbomoleOutputProcessing import turbomoleOutputProcessing as top
import os.path
from multiprocessing import Pool
from functools import partial
import sys
import configparser
from scipy.optimize import curve_fit



def eval_T(disp_index, para):
	"""
	evaluates zeroth order Greens function to estimate the Transmission in wide band limit. Function is used in multiprocessing

	Args:
		param1 (int): disp_index. Specifies which index in dir list should be processed
		param2 (tuple): param2[0] list of directories, param2[1] minimal displacement, param2[2] n_occupied
	Returns:
		(float:estimate for t, int: processed displacement index) 
	"""		
	directories = para[0]
	min_displacement=para[1]
	n_occupied = para[2]

	directory = directories[disp_index+min_displacement]
	#print("index " + str())
	try:
		eigenvalues, eigenvectors = top.read_mos_file(directory + "/mos")
	except ValueError as e: 
		print("error " + str(e))
		return 0.0, float(disp_index)

	s_range = top.find_c_range_atom("s", 1, directory + "/coord")
	r_range = top.find_c_range_atom("s", 2, directory + "/coord")

	e_f = (eigenvalues[n_occupied-1]+eigenvalues[n_occupied])/2.

	eigenvectors = np.asmatrix(eigenvectors)

	denominator_list = list()
	numerator_list = list()
	G0_sr_k = list()
	G0_sr = 0
	delta = 1E-12
	#delta = 0.0
	
	for i in range(0,len(eigenvalues)):
		numerator = (eigenvectors[r_range[0]:r_range[1],i])*np.transpose(eigenvectors[s_range[0]:s_range[1],i])
		denominator = e_f - eigenvalues[i] + 1.0j* delta

		denominator_list.append(denominator)
		numerator_list.append(numerator)
		tmp = numerator/denominator
		#tmp = np.exp(-(i-(n_occupied-1))**2/(2*50^2))*tmp

		G0_sr += tmp

	T_est = np.real(np.trace(G0_sr*np.conj(np.transpose(G0_sr))))
	return T_est, float(disp_index)


def plot_T_vs_d(calc_path, molecule_name, n_occupied, config_path):
	"""
	evaluates zeroth order Greens function to estimate the Transmission in wide band limit. T estimate is plotted vs displacement.

	Args:
		param1 (String): Path to turbomole calculation
		param2 (String): molecule name generation_individual
		param3 (int): occupied mos -> define fermi energy
		param4 (String): Path to config file
	Returns:
		
	"""		
	
	print("config path " + str(config_path))
	#process disp_neg
	list_neg = os.listdir(calc_path + "/disp_neg/")
	list_neg = sorted(list_neg, reverse=True)		
	list_neg = [calc_path + "/disp_neg/" + s for s in list_neg if (os.path.isdir(calc_path + "/disp_neg/" + s)) and (s != '0000')]
	list_neg = [s for s in list_neg if os.path.exists(s+"/GEO_OPT_CONVERGED")]

	#process disp_pos
	list_pos = os.listdir(calc_path + "/disp_pos/")
	list_pos = sorted(list_pos)	
	list_pos = [calc_path + "/disp_pos/" + s for s in list_pos if os.path.isdir(calc_path + "/disp_pos/" + s)]
	list_pos = [s for s in list_pos if os.path.exists(s+"/GEO_OPT_CONVERGED")]

	#all folder which will be processed
	folders = list(list_neg)
	folders.extend(list_pos)

	displacement = list()
	g0_list = list()


	disp_indices = np.linspace(-len(list_neg), len(list_pos),len(list_pos)+len(list_neg), False, dtype=int)

	p = Pool(32)
	min_disp = len(list_neg)	
	eingabe = (folders, min_disp, n_occupied) 

	cfg = configparser.ConfigParser()
	cfg.read(config_path)
	displacement = float(cfg.get('Turbomole Calculations', 'displacement_per_step')) 

	#multiprocessing
	result = p.map(partial(eval_T, para=eingabe), disp_indices)

	T_est, disp = [], []
	for x, y in result:
	    T_est.append(x)
	    disp.append(float(y)*displacement)

	#fit data
	minimum = np.argmin(T_est)
	min_disp = disp[minimum]
	min_T_est = T_est[minimum]

	def func(x, a):
		return (1-np.exp(-a*np.abs(x-min_disp)+np.log(1-min_T_est)))
	try:
		popt, pcov = curve_fit(func, disp, T_est)
		print(pcov)

	except RuntimeError:
		print("Fit to wrong model")

	fitted_data = func(np.asarray(disp), popt[0])


	plt.plot(disp, T_est)
	plt.xlabel('Displacement [$\mathrm{\AA}$]',fontsize=20)
	plt.ylabel('$\mathrm{T}_{\mathrm{estimate}}$',fontsize=20)
	#plt.xlim(-2,4)
	#plt.ylim(15,50)
	#ACHTUNG !!!: x achse stimmt nicht ganz zuordnung!!
	plt.savefig(calc_path + "/" + molecule_name + "_T_estimate.pdf", bbox_inches='tight')
	plt.savefig(calc_path + "/" + molecule_name + "_T_estimate.svg", bbox_inches='tight')

	plt.yscale("log")
	#plt.ylim(10E-7,1)
	plt.savefig(calc_path + "/" + molecule_name + "_T_estimate_log.pdf", bbox_inches='tight')
	plt.savefig(calc_path + "/" + molecule_name + "_T_estimate_log.svg", bbox_inches='tight')
	top.write_plot_data(calc_path + "/" + molecule_name +  "_T_estimate.dat",(disp,T_est), "disp, T_est")
	plt.plot(disp, fitted_data)
	plt.savefig(calc_path + "/" + molecule_name + "_T_estimate_log_fit.pdf", bbox_inches='tight')
	plt.savefig(calc_path + "/" + molecule_name + "_T_estimate_log_fit.svg", bbox_inches='tight')
	file = open(calc_path + "/" + molecule_name + "_T_estimate_params.dat", "w")
	file.write(str(float(popt[0])))
	file.write("\n")
	file.write(str(float(pcov[0])))
	file.write("\n")
	file.write(str(min_T_est))
	file.write("\n")
	file.write(str(min_disp))
	file.close()

def plot_energy_levels(calc_path, molecule_name, n_occupied, config_path):
	"""
	Plots energy levels homo-4 ... Lumo+4 vs displacement

	Args:
		param1 (String): Path to turbomole calculation
		param2 (String): molecule name generation_individual
		param3 (int): occupied mos -> define fermi energy
		param4 (String): Path to config file
	Returns:
		
	"""	
	number_occupied = n_occupied

	#process disp_neg
	list_neg = os.listdir(calc_path + "/disp_neg/")
	list_neg = sorted(list_neg, reverse=True)	
	list_neg = [calc_path + "/disp_neg/" + s for s in list_neg if (os.path.isdir(calc_path + "/disp_neg/" + s)) and (s != '0000')]
	list_neg = [s for s in list_neg if os.path.exists(s+"/GEO_OPT_CONVERGED")]
	#process disp_pos
	list_pos = os.listdir(calc_path + "/disp_pos/")
	list_pos = sorted(list_pos)	
	list_pos = [calc_path + "/disp_pos/" + s for s in list_pos if os.path.isdir(calc_path + "/disp_pos/" + s)]
	list_pos = [s for s in list_pos if os.path.exists(s+"/GEO_OPT_CONVERGED")]
	#process disp_pos]
	reference_dir = list_pos[0]

	folders = list(list_neg)
	folders.extend(list_pos)
	displacement = list()

	disp_indices = np.linspace(-len(list_neg), len(list_pos),len(list_pos)+len(list_neg), False, dtype=int)
	homo_list = list()
	homo_list_m1 = list()
	homo_list_m2 = list()
	homo_list_m3 = list() 
	homo_list_m4 = list()
	lumo_list = list()
	lumo_list_p1 = list()
	lumo_list_p2 = list()
	lumo_list_p3 = list()
	lumo_list_p4 = list()

	e_fermi_list = list()

	cfg = configparser.ConfigParser()
	cfg.read(config_path)
	har2Ev = float(cfg.get('Building Procedure', 'har2Ev'))
	disp_per_step = float(cfg.get('Turbomole Calculations', 'displacement_per_step')) 
	

	for i in range(0, len(folders)):

		directory = folders[i]
		print(directory)
		eigenvalues, eigenvectors = top.read_mos_file(directory + "/mos")
		
		homo_list_m1.append(float(eigenvalues[number_occupied-1-1])*har2Ev)
		homo_list_m2.append(float(eigenvalues[number_occupied-1-2])*har2Ev)
		homo_list_m3.append(float(eigenvalues[number_occupied-1-3])*har2Ev)
		homo_list_m4.append(float(eigenvalues[number_occupied-1-4])*har2Ev)
		homo_list.append(float(eigenvalues[number_occupied-1])*har2Ev)
		lumo_list_p1.append(float(eigenvalues[number_occupied+1])*har2Ev)
		lumo_list_p2.append(float(eigenvalues[number_occupied+2])*har2Ev)
		lumo_list_p3.append(float(eigenvalues[number_occupied+3])*har2Ev)
		lumo_list_p4.append(float(eigenvalues[number_occupied+4])*har2Ev)
		#homo_list.append(float(eigenvalues[number_occupied-1])*27.211)
		lumo_list.append(float(eigenvalues[number_occupied])*har2Ev)

		e_fermi_list.append((float(eigenvalues[number_occupied-1])+float(eigenvalues[number_occupied]))*har2Ev/2.0)

		displacement.append(float(disp_indices[i])*disp_per_step)
	fig, ax = plt.subplots(1)

	ax.plot(displacement, lumo_list_p4, label="$E_{lumo+4}$")
	ax.plot(displacement, lumo_list_p3, label="$E_{lumo+3}$")
	ax.plot(displacement, lumo_list_p2, label="$E_{lumo+2}$")
	ax.plot(displacement, lumo_list_p1, label="$E_{lumo+1}$")
	ax.plot(displacement, lumo_list, label="$E_{lumo}$", color="blue")
	ax.plot(displacement, e_fermi_list, label="$E_{f}$", color="black", linestyle='dashed')
	ax.plot(displacement, homo_list, label="$E_{homo}$", color="green")
	ax.plot(displacement, homo_list_m1, label="$E_{homo-1}$")
	ax.plot(displacement, homo_list_m2, label="$E_{homo-2}$")
	ax.plot(displacement, homo_list_m3, label="$E_{homo-3}$")
	ax.plot(displacement, homo_list_m4, label="$E_{homo-4}$")	
	
	ax.set_xlabel('Displacement ($\mathrm{\AA}$)',fontsize=20)
	ax.set_ylabel('Energy (eV)',fontsize=20)
	ax.set_yscale('log')
	ax.legend()
	plt.savefig(calc_path + "/" + molecule_name + "_energy_levels.pdf", bbox_inches='tight')
	plt.savefig(calc_path + "/" + molecule_name + "_energy_levels.svg", bbox_inches='tight')
	top.write_plot_data(calc_path + "/" + molecule_name + "_energy_levels.dat", (np.round(displacement,2), homo_list_m4, homo_list_m3,homo_list_m2,homo_list_m1,homo_list, lumo_list,lumo_list_p1,lumo_list_p2,lumo_list_p3,lumo_list_p4), "displacement, homo-4 (eV) , ... , lumo +4")

if __name__ == '__main__':
	
	#argv[1] : calc path: where disp_pos and disp_neg is located
	#argv[2] : moleculename (prefix of plot files)
	#argv[3] : occupied orbitals (e_f = (E_homo+E_lumo)/2)
	#argv[4] : config_path
	#print(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
	plot_T_vs_d(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
	plot_energy_levels(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])