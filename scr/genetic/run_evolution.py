import configparser
import os
import os.path
from os import path
import run_generation as rg

def run_evolution(config_path, calculation_path):
	"""
	Runs evolution. Inits calculation dirs and invokes evaluation of first generation

	Args:
		param1 (String) : Path to config file
		param2 (String) : Path to calculation

	Returns:
		error codes (int)
	"""
	#check config file
	if(path.exists(config_path) == False):
		print("Faulty config file")
		return -1


	generation_data_path = calculation_path + "/" "generation_data"
	#create calculation files
	try:
		#create generation dir
		if(path.exists(generation_data_path) == False):
			os.mkdir(calculation_path)
			generation_zero_path = generation_data_path + "/0"
			os.mkdir(generation_zero_path)
	except OSError:
		print ("Creation of the directory %s failed" % generation_data_path)
		return -1

	rg.run_generation(0, config_path, calculation_path)

def test(a):
	"""

	Args:
		a (object):
	"""
	pass

if __name__ == '__main__':
	config_path = "/alcc/gpfs2/home/u/blaschma/genetic_run_25/config"
	calculation_path = "/alcc/gpfs2/home/u/blaschma/genetic_run_25/"
	run_evolution(config_path, calculation_path)