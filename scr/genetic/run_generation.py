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
		population = ga.run_generation(
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
		generation=generation,
		fitness_value=0)

	#every other generation
	else:
		#generation-1 because prevois generation should be read
		population, fitness_value = read_population(generation-1,config_path, calculation_path)
		print(fitness_value)
		population = ga.run_generation(
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
		generation=generation,
		fitness_value=fitness_value)



		#print(population)
		#print(fitness_value)
	write_generation(population,generation, config_path, calculation_path)

def write_generation(population, generation, config_path, calculation_path):
	"""
    write current generation to calculation_path/current_generation.dat and to generation_data/generation/generation_summary.dat.
    First file contains only information about the population. The second file contains the fittness value, too. In addition a file 
    calculation_path/generation.dat which contains the current number of generations

    Args:
    		param1 (Population) : population to write
    		param2 (int) : Generation
            param3 (String) : Path to config file
            param4 (String) : Path to calculation              
    Returns:
            
    """

	try:
		first_file_path = calculation_path + "/curr_population.dat"
		second_file_path = calculation_path + "/generation_data/" + str(generation)
		first_file = open(first_file_path, "w")
		generation_file = open(calculation_path + "/generation.dat", "w")
		if(path.exists(second_file_path) == False):
			os.mkdir(second_file_path)	
		second_file = open(second_file_path + "/summary.dat", "w")	
	except OSError as e:
		print("log file cannot be opened " + str(e))
		return 1	
		
	for i in range(0, len(population)):		
		first_file.write(str(population[i]).replace('[', "").replace("]", "")+ "\n")
		#second_file.write(str(fitness_value[i]) + "		" + str(population[i]).replace('[', "").replace("]", "")+ "\n")
		second_file.write(str(population[i]).replace('[', "").replace("]", "")+ "\n")

	generation_file.write(str(generation))
	generation_file.close()

	first_file.close()
	second_file.close()

def read_population(generation, config_path, calculation_path):
	"""
    read current generation individuals and their fitness from calculation path

    Args:
    	param1 (int) : number of generation which should be read
		param2 (String) : Path to config file
		param3 (String) : Path to calculation              
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

def next_generation(config_path, calculation_path):
	"""
    Invokes the next generation. The number of generation is read from calculation_path/generation.dat

    Args:
    	param1 (String) : Path to config file
		param2 (String) : Path to calculation              
    Returns:
            
    """
	try:
		filename = calculation_path + "/generation.dat"	
		file = open(filename)
		for line in file:
			generation = int(line)

	except (OSError,ValueError) as e:
		print("generation file cannot be found or generation file is faulty " + str(e))


	#check if generation was correctly processed
	file_to_check = calculation_path + "/generation_data/" + str(generation) + "/fitness.dat"
	if(os.path.exists(file_to_check) == False):
		print("Generation " + str(generation) + " was not processed correctly")
		return -1

	#increase number of genrations
	generation += 1
	run_generation(generation, config_path, calculation_path)


if __name__ == '__main__':
	#config_path= "C:/Users/bergb/OneDrive/Studium/Master/Masterarbeit/Code/genetic_algorithm/genetic_algorithm/test/config"
	#calculation_path = "C:/Users/bergb/OneDrive/Studium/Master/Masterarbeit/Code/genetic_algorithm/genetic_algorithm/test/"
	config_path = "/alcc/gpfs2/home/u/blaschma/test/config"
	calculation_path = "/alcc/gpfs2/home/u/blaschma/test/"
	#for i in range(0,10000):
	run_generation(0, config_path, calculation_path)
	#print(population)
	next_generation(config_path, calculation_path)
	#read_population(config_path, calculation_path)


