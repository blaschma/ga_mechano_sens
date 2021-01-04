import configparser
import os
import os.path
from os import path
import sys
import sys
import genetic_algorithm as ga
from functools import partial

def run_generation(generation : int, config_path, calculation_path):
	"""
    Runs evolution. Inits calculation dirs and invokes evaluation of first generation

    Args:
    		param1 (int) : number of generation
            param2 (String) : Path to config file
            param3 (String) : Path to calculation
            
            

    Returns:
            error codes (int)
    """

    #load specification from config file 
	cfg = configparser.ConfigParser()
	cfg.read(config_path)
	population_size = int(cfg.get('Genetic Algorithm', 'population_size'))
	generation_limit = int(cfg.get('Genetic Algorithm', 'generation_limit'))	
	genome_length = int(cfg.get('Genetic Algorithm', 'genome_length'))	 
	evaluation_methods_path = cfg.get('Genetic Algorithm', 'evaluation_methods_path')

	#check if default setting for evaluation methods
	if(evaluation_methods_path == "default"):
		#sys.path.insert(1, '../problem_specification/')
		sys.path.append(os.path.realpath('..'))
		from problem_specification import bbEV 
		ev = bbEV.bbEv(generation, 0) #-> todo individual processing
		pass
		
	#load specific evaluation methods
	else:
		#load methods for genetic algorithm
		sys.path.append(evaluation_methods_path)
		import evaluation_methods as evaluation_methods

	#first generation
	if(generation == 0):
		population, fitness_value = ga.run_generation(
		populate_func=partial(
			ev.generate_population, size=population_size, genome_length=genome_length
			),
		fitness_func=ev.fitness,
		selection_func=ev.selection_pair,			
		crossover_func=ev.crossover,			
		mutation_func=ev.mutation,			
		population=0,
		fitness_limit=5,
		generation_limit=generation_limit,
		population_size=population_size,
		generation=generation)
	#every other generation
	else:
		population = read_population(config_path, calculation_path)
		population, fitness_value = ga.run_generation(
		populate_func=partial(
			ev.generate_population, size=population_size, genome_length=genome_length
			),
		fitness_func=ev.fitness,
		selection_func=ev.selection_pair,			
		crossover_func=ev.crossover,			
		mutation_func=ev.mutation,			
		population=population,
		fitness_limit=5,
		generation_limit=generation_limit,
		population_size=population_size,
		generation=generation)



		#print(population)
		#print(fitness_value)
	write_generation(population, fitness_value, generation, config_path, calculation_path)

def write_generation(population, fitness_value,generation, config_path, calculation_path):
	"""
    write current generation to calculation_path/current_generation-dat and to generation_data/generation/generation_summary.dat.
    First file contains only information about the population. The second file contains the fittness value, too.

    Args:
    		param1 (Population) : population to write
    		param2 (List(float)) : Fittness values
    		param3 (int) : Generation
            param4 (String) : Path to config file
            param5 (String) : Path to calculation              
    Returns:
            
    """

	try:
		first_file_path = calculation_path + "/curr_population.dat"
		second_file_path = calculation_path + "/generation_data/" + str(generation)
		first_file = open(first_file_path, "w")
		if(path.exists(second_file_path) == False):
			os.mkdir(second_file_path)	
		second_file = open(second_file_path + "/summary.dat", "w")	
	except OSError as e:
		print("log file cannot be opened " + str(e))
		return 1	
		
	for i in range(0, len(population)):		
		first_file.write(str(population[i]).replace('[', "").replace("]", "")+ "\n")
		second_file.write(str(fitness_value[i]) + "		" + str(population[i]).replace('[', "").replace("]", "")+ "\n")
	first_file.close()
	second_file.close()

def read_population(config_path, calculation_path):
	"""
    read current generation from calculation path

    Args:
    		param1 (Population) : population to write
    		param2 (List(float)) : Fittness values
    		param3 (int) : Generation
            param4 (String) : Path to config file
            param5 (String) : Path to calculation              
    Returns:
            
    """
	filename = calculation_path + "/curr_population.dat"	
	file = open(filename)
	population = list()
	for line in file:
		#print(line)
		tmp = line.strip().split(", ")
		tmp = [int(tmp[i]) for i in range(0,len(tmp))]
		population.append(tmp)

	return population





if __name__ == '__main__':
	config_path= "C:/Users/bergb/OneDrive/Studium/Master/Masterarbeit/Code/genetic_algorithm/genetic_algorithm/test/config"
	calculation_path = "C:/Users/bergb/OneDrive/Studium/Master/Masterarbeit/Code/genetic_algorithm/genetic_algorithm/test/"
	#for i in range(0,10000):
	#run_generation(0, config_path, calculation_path)
	#print(population)
	#run_generation(2, config_path, calculation_path)
	#read_population(config_path, calculation_path)


