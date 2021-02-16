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
import configparser
from matplotlib.ticker import MaxNLocator


def read_population(generation, config_path, calculation_path):
	"""
    read current generation individuals and their fitness from calculation path

    Args:
    	param1 (int): number of generation which should be read
		param2 (String): Path to config file
		param3 (String): Path to calculation              
    Returns:
    	(population, fitness_values)
            
    """
	try:
		filename_population = calculation_path + "/curr_population.dat"	
		file_population = open(filename_population)
		filename_fitness = calculation_path + "/generation_data/" + str(generation) + "/fitness.dat"	
		file_fitness = open(filename_fitness)
	except OSError as e:
		print("Cannot open file " + str(e))
		return -1,-1

	#read population
	population = list()
	for line in file_population:
		#print(line)
		tmp = line.strip().split(", ")
		tmp = [int(tmp[i]) for i in range(0,len(tmp))]
		population.append(tmp)

	#read fitness values
	fitness_value = list()
	for line in file_fitness:
		#print(line)
		try:
			tmp = float(line)
		except ValueError as e:
			print("Faulty fitness file " + str(e))
			return -1,-1
		fitness_value.append(tmp)
		

	return population, fitness_value


def read_generation(config_path, calculation_path):
	"""
    Invokes the next generation. The number of generation is read from calculation_path/generation.dat

    Args:
    	param1 (String): Path to config file
		param2 (String): Path to calculation       

    Returns:
            
    """
	try:
		filename = calculation_path + "/generation.dat"	
		file = open(filename)
		for line in file:
			generation = int(line)

	except (OSError,ValueError) as e:
		print("generation file cannot be found or generation file is faulty " + str(e))

	return generation





if __name__ == '__main__':
	# sys.argv[1] path to process
	# sys.argv[2] config path
	# sys.argv[3] figure path
	calculation_path = sys.argv[1]
	config_path = sys.argv[2]
	#how many generations have been calculated
	generation = read_generation(config_path, calculation_path)
	if(generation == None):
		exit()
	
	generations_to_check = np.arange(0, generation)

	#load fitness values
	fitness_values = list()
	fitness_means=list()
	for i in generations_to_check:
		fitness_value = read_population(i, config_path, calculation_path)
		fitness_values.append(fitness_value[1])
		fitness_means.append(np.mean(fitness_value[1]))
	fig, ax = plt.subplots(1)
	num_individuals = len(fitness_value[1])
	#color & plotting stuff
	phi = np.linspace(0, 2*np.pi, num_individuals)
	rgb_cycle = np.vstack((            # Three sinusoids
    .5*(1.+np.cos(phi          )), # scaled to [0,1]
    .5*(1.+np.cos(phi+2*np.pi/3)), # 120° phase shifted.
    .5*(1.+np.cos(phi-2*np.pi/3)))).T # Shape = (60,3)

	for xe, ye in zip(generations_to_check, fitness_values):
		ax.scatter([xe] * len(ye), ye, c=rgb_cycle, s=num_individuals, marker="x")

	ax.plot(generations_to_check, fitness_means)
	ax.set_xlabel('Generation',fontsize=20)
	ax.set_ylabel('Fitness values',fontsize=20)
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	ax.set_yscale('symlog')

	plt.savefig( sys.argv[3] + "/fitness_values.pdf", bbox_inches='tight')