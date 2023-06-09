import tmoutproc as top
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.optimize import curve_fit

def load_and_plot(molecule_name,n_atoms):
	"""
	loads total energy file and plots it vs displacement. fit is also done and written to STIFFNESS file

		param1 (String) : molecule_name
		param2 (int) : n_atoms for normalization

	Returns:
		
	"""

	#num atoms is 3 smaller beacuse keywords are counted, too
	n_atoms = n_atoms -3

	#load data
	try:
		data = top.read_plot_data(molecule_name + "_totalEnergy.dat")
	except ValueError as e:
		print(e)
		file = open(molecule_name + "_stiffness.dat", "w")
		error_value_fit = -100
		error_value_std = 1000
		file.write("a b c std(a) std(b) std(c)  -> a(x-b)**2+c (first line:normalized, second line: not normalized)" + "\n")
		file.write(str(error_value_fit) + "	" + str(error_value_fit) + "	" + str(error_value_fit) + "	" + str(error_value_std) + "	" + str(error_value_std) + "	" + str(error_value_std) + "\n")
		file.write(str(error_value_fit) + "	" + str(error_value_fit) + "	" + str(error_value_fit) + "	" + str(error_value_std) + "	" + str(error_value_std) + "	" + str(error_value_std))
		file.close()
	else:
		zipped_lists = zip(data[0], data[1])
		sorted_pairs = sorted(zipped_lists)
		tuples = zip(*sorted_pairs)

		disp, energy = [ list(tuple) for tuple in  tuples]
		
		disp = np.asarray(disp, dtype = float)
		energy = np.asarray(energy, dtype = float)


		#plotting and fitting
		interesting_energy_minimum = np.min(energy)
		energy = energy-interesting_energy_minimum
		disp = disp*0.1
		energy_normalized = energy * (1./float(n_atoms)+0.0)
		def func(x, a,b,c):
			return a * (x-b)**2 +c
		try:
			popt_normalized, pcov_normalized = curve_fit(func, disp, energy_normalized)
			popt, pcov = curve_fit(func, disp, energy)


		except RuntimeError:
			print("Fit to wrong model")
			fit = [-100,-100,-100]
			std = [1000,1000,1000]
			file = open(molecule_name + "_stiffness.dat", "w")
			file.write("a b c std(a) std(b) std(c)  -> a(x-b)**2+c (first line:normalized, second line: not normalized)" + "\n")
			file.write(str(fit[0]) + "	" + str(fit[1]) + "	" + str(fit[2]) + "	" + str(std[0]) + "	" + str(std[1]) + "	" + str(std[2]) + "\n")
			file.write(str(fit[0]) + "	" + str(fit[1]) + "	" + str(fit[2]) + "	" + str(std[0]) + "	" + str(std[1]) + "	" + str(std[2]))
			file.close()
		else:
			fit_normalized = popt_normalized
			fit = popt
			fitted_data_normalized = fit_normalized[0] * (disp-fit_normalized[1])**2+fit_normalized[2]
			std_normalized = np.sqrt(np.diag(pcov_normalized))
			std=np.sqrt(np.diag(pcov))
			print("a=" + str(fit_normalized[0]))

			#write fit data to file
			file = open(molecule_name + "_stiffness.dat", "w")
			file.write("a b c std(a) std(b) std(c)  -> a(x-b)**2+c (first line:normalized, second line: not normalized)" + "\n")
			file.write(str(fit_normalized[0]) + "	" + str(fit_normalized[1]) + "	" + str(fit_normalized[2]) + "	" + str(std_normalized[0]) + "	" + str(std_normalized[1]) + "	" + str(std_normalized[2]) + "\n")
			file.write(str(fit[0]) + "	" + str(fit[1]) + "	" + str(fit[2]) + "	" + str(std[0]) + "	" + str(std[1]) + "	" + str(std[2]))
			file.close()

			plt.plot(disp, energy_normalized)
			plt.plot(disp, fitted_data_normalized, '--')
			
			plt.title("a=" + str(fit_normalized[0]) + " std=" + str(std_normalized))
			plt.ylabel('($(E-E_0)/N$) [H]',fontsize=20)
			plt.xlabel('Displacement [$\mathrm{\AA}$]',fontsize=20)
			plt.savefig(molecule_name + "_totalEnergy.pdf", bbox_inches='tight')
			plt.savefig(molecule_name + "_totalEnergy.svg", bbox_inches='tight')

	


if __name__ == '__main__':
	#sys.argv[1]: path/moleculename
	#sys.argv[2]: n_atoms +3 -> because of counting of keywords
	load_and_plot(sys.argv[1],int(sys.argv[2]))

	
