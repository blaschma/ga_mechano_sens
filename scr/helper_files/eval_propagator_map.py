import numpy as np
import matplotlib
from matplotlib.colors import LogNorm

matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from turbomoleOutputProcessing import turbomoleOutputProcessing as top
import os.path
from multiprocessing import Pool
from functools import partial
import sys
import configparser
from scipy.optimize import curve_fit

har2eV = 27.2114

def eval_T(disp_index, para):
	"""
	evaluates zeroth order Greens function to estimate the Transmission in wide band limit. Function is used in multiprocessing

	Args:
		param1 (int): disp_index. Specifies which index in dir list should be processed
		param2 (tuple): params
	Returns:
		(float:estimate for t, int: processed displacement index) 
	"""		
	directories = para[0]
	min_displacement=para[1]
	n_occupied = para[2]
	E_min = para[3]
	E_max = para[4]
	n_E_steps = para[5]

	E = np.linspace(E_min, E_max, n_E_steps)

	directory = directories[disp_index+min_displacement]
	try:
		eigenvalues, eigenvectors = top.read_mos_file(directory + "/mos")
	except ValueError as e: 
		print("error " + str(e))
		return 0.0, float(disp_index)

	s_range = top.find_c_range_atom("s", 1, directory + "/coord")
	r_range = top.find_c_range_atom("s", 2, directory + "/coord")

	eigenvectors = np.asmatrix(eigenvectors)

	delta = 1E-12
	T_est = list()
	numerator = np.asarray([(eigenvectors[r_range[0]:r_range[1],i])*np.transpose(eigenvectors[s_range[0]:s_range[1],i]) for i in range(0,len(eigenvalues))])
	denominator = np.asarray([(E - eigenvalues[i] + 1.0j* delta)  for i in range(0,len(eigenvalues))])
	for j, energy in enumerate(E):
		G0_sr = [np.asmatrix(numerator[i]/denominator[i][j]) for i in range(0,len(eigenvalues))]
		G0_sr = np.asmatrix(np.sum(G0_sr, axis=0))
		T_est_ = np.real(np.trace(G0_sr * np.conj(np.transpose(G0_sr))))
		T_est.append(T_est_)

	return T_est, float(disp_index)


def plot_T_vs_d_energy_resolved(calc_path, molecule_name, n_occupied, config_path, e_range=5, n_E_steps=1000):
	"""
	evaluates zeroth order Greens function to estimate the Transmission in wide band limit. T estimate is plotted vs displacement.

	Args:
		param1 (String): Path to turbomole calculation
		param2 (String): molecule name generation_individual
		param3 (int): occupied mos -> define fermi energy
		param4 (String): Path to config file
		param5 (float): Energy range from fermi energy (optional)
		param6 (float): n_energy_steps (optional)
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

	#determine Fermi energy
	try:
		eigenvalues, eigenvectors = top.read_mos_file(list_pos[0] + "/mos")
		e_f = (eigenvalues[n_occupied - 1] + eigenvalues[n_occupied]) / 2.
		print(f"Fermi energy={e_f}")
	except ValueError as e:
		print("error " + str(e))


	disp_indices = np.linspace(-len(list_neg), len(list_pos),len(list_pos)+len(list_neg), False, dtype=int)


	min_disp = len(list_neg)

	E_min = e_f-e_range/har2eV
	E_max = e_f+e_range/har2eV

	eingabe = (folders, min_disp, n_occupied, E_min, E_max, n_E_steps)

	cfg = configparser.ConfigParser()
	cfg.read(config_path)
	displacement = float(cfg.get('Turbomole Calculations', 'displacement_per_step')) 

	#multiprocessing
	p = Pool(16)
	result = p.map(partial(eval_T, para=eingabe), disp_indices)

	T_est, disp = [], []
	E = np.linspace(E_min, E_max, n_E_steps)

	for x, y in result:
		T_est.extend(x)
		disp.append(float(y)*displacement)
	T_est = np.reshape(T_est, (-1, len(E)))
	fig, ax1 = plt.subplots(1)
	cs = ax1.imshow(T_est.T, origin='lower', extent=[min(disp), max(disp), (E_min-e_f)*har2eV, (E_max-e_f)*har2eV], norm=LogNorm(),
					aspect='auto', interpolation='None')

	plt.xlabel(r'Displacement $\AA$', fontsize=20)
	plt.ylabel('Energy (eV)', fontsize=20)
	cb = fig.colorbar(cs, ax=ax1)
	cb.ax.tick_params(labelsize=15)
	cb.set_label(r"$|\mathrm{G_{l,r}^{r/a}}|^2$", fontsize=20)
	plt.tick_params(axis="x", labelsize=15)
	plt.tick_params(axis="y", labelsize=15)

	plt.savefig(calc_path + "/" + molecule_name + "_T_estimate_map.pdf", bbox_inches='tight')
	plt.savefig(calc_path + "/" + molecule_name + "_T_estimate_map.svg", bbox_inches='tight')


if __name__ == '__main__':
	
	#argv[1] : calc path: where disp_pos and disp_neg is located
	#argv[2] : moleculename (prefix of plot files)
	#argv[3] : occupied orbitals (e_f = (E_homo+E_lumo)/2)
	#argv[4] : config_path
	# argv[5] : e_range in eV (optional)
	# argv[6] : n_E_steps (optional)
	print(len(sys.argv))
	if(len(sys.argv)==5):
		plot_T_vs_d_energy_resolved(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
	elif(len(sys.argv)==6):
		plot_T_vs_d_energy_resolved(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], float(sys.argv[5]))
	elif(len(sys.argv)==7):
		plot_T_vs_d_energy_resolved(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], float(sys.argv[5]), int(sys.argv[6]))
	else:
		print("Check your call")



