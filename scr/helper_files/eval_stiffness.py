from turbomoleOutputProcessing import turbomoleOutputProcessing as top 
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
	data, headers = top.read_plot_data(molecule_name + "_totalEnergy.dat")
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
	energy = energy * (1./float(n_atoms)+0.0)
	def func(x, a,b,c):
		return a * (x-b)**2 +c
	try:
		popt, pcov = curve_fit(func, disp, energy)
	except RuntimeError:
		print("Fit to wrong model")
		fit = [-100,-100,-100]
		std = [1000,1000,1000]
	fit = popt
	fitted_data = fit[0] * (disp-fit[1])**2+fit[2]
	std = np.sqrt(np.diag(pcov))
	print("a=" + str(fit[0]))

	#write fit data to file
	file = open("../stiffness.dat", "w")
	file.write("a b c std(a) std(b) std(c)  -> a(x-b)**2+c" + "\n")
	file.write(str(fit[0]) + "	" + str(fit[1]) + "	" + str(fit[2]) + "	" + str(std[0]) + "	" + str(std[1]) + "	" + str(std[2]))
	file.close()

	plt.plot(disp, energy)
	plt.plot(disp, fitted_data, '--')
	
	plt.title("a=" + str(fit[0]) + " std=" + str(std))
	plt.ylabel('($(E-E_0)/N$) [H]',fontsize=20)
	plt.xlabel('Displacement [$\mathrm{\AA}$]',fontsize=20)
	plt.savefig(molecule_name + "_totalEnergy.pdf", bbox_inches='tight')
	plt.savefig(molecule_name + "_totalEnergy.svg", bbox_inches='tight')

if __name__ == '__main__':
	#sys.argv[1]: molecule name
	#sys.argv[2]: n_atoms +3 -> because of counting of keywords
	load_and_plot(sys.argv[1],int(sys.argv[2]))

	
