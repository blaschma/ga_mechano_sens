from random import choices, randint, randrange, random
from typing import List, Callable, Tuple
from collections import namedtuple
from functools import partial
import numpy as np
import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import copy


Genome = List[int]
Population = List[Genome]
FitnessFunc = Callable[[Genome, int, int], float]
PopulateFunc = Callable[[], Population]
SelectionFunc = Callable[[Population, FitnessFunc],Tuple[Genome, Genome]]
CrossoverFunc = Callable[[Genome,Genome], Tuple[Genome, Genome]]
MutationFunc = Callable[[Genome], Genome]


#"""
def run_evolution(
		populate_func: PopulateFunc,
		fitness_func: FitnessFunc,
		fitness_limit: int,
		selection_func : SelectionFunc,
		crossover_func : CrossoverFunc,
		mutation_func: MutationFunc,
		generation_limit: int = 100
	) -> Tuple[Population, int]:
	#list of fitness_values
	fitness_value = list()

	#initialize first population
	population = populate_func()

	#geration loop
	for i in range(generation_limit):
		print("generation " + str(i))

		#evaluate pupulation
		population = sorted(
			population,
			key = lambda genome: fitness_func(genome)			
			)
		fitness_value.append(fitness_func(population[0]))

		#check if fittness_limit is exceeded
		"""
		if fitness_func(population[0]) >= fitness_limit:
			break
		"""

		#take best individuals of generation and..
		next_generation = population[0:2]

		#...fill generation with mutated and cross over children
		for j in range(int(len(population)/2)-1):
			#select parent according to selection function
			parents = selection_func(population, fitness_func)
			#combine features of parents to generate offspring
			offspring_a, offspring_b = crossover_func(parents[0], parents[1])
			#mutate offspring
			offspring_a = mutation_func(offspring_a)
			offspring_b = mutation_func(offspring_b)
			#add offspring to generation
			next_generation += [offspring_a, offspring_b]
		population = next_generation

	#sort final population 
	population = sorted(
			population,
			key = lambda genome: fitness_func(genome)
			)
	return population, i,fitness_value



def run_generation(
		populate_func: PopulateFunc,
		fitness_func: FitnessFunc,		
		selection_func : SelectionFunc,
		crossover_func : CrossoverFunc,
		mutation_func: MutationFunc,
		population : Population,
		fitness_limit: int,
		generation_limit: int,
		population_size:int,
		generation: int,
		fitness_value
	) :
	#list of fitness_values

	

	#initialize first population
	if(generation == 0):
		population = populate_func()
		print("populated population " +str(population))
		#invoke fitness evaluation
		population_for_fitness_eval = copy.deepcopy(population)
		for i in range(0, population_size):
			fitness_func(population_for_fitness_eval[i],generation,i)
		print("population after fitness_func " +str(population))
		return population, ""


	#sort population
	zipped_lists = zip(fitness_value, population)
	sorted_pairs = sorted(zipped_lists, reverse=True)
	tuples = zip(*sorted_pairs)
	fitness_value, population = [ list(tuple) for tuple in  tuples]


	#check if fittness_limit is exceeded
	"""
	if fitness_func(population[0]) >= fitness_limit:
		break
	"""
	#family register
	family_register=""

	#take best individuals of generation and..
	next_generation = population[0:2]
	family_register += str(population[0]) + "\n"
	family_register += str(population[1]) + "\n"

	#...fill generation with mutated and cross over children
	for j in range(int(len(population)/2)-1):
		#select parent according to selection function
		parents = selection_func(population, fitness_value)#->todo too inefficent
		#combine features of parents to generate offspring
		print("parents[0] " + str(parents[0]) + " parents[1] " + str(parents[1]))
		offspring_a, offspring_b = crossover_func(parents[0], parents[1])
		offspring_a_save = str(offspring_a)
		offspring_b_save = str(offspring_b)

		#mutate offspring
		print("offspring_a " + str(offspring_a))
		print("offspring_b " + str(offspring_b))		
		offspring_a = mutation_func(offspring_a)
		offspring_b = mutation_func(offspring_b)

		#handle family register
		offspring_a_mutation=""
		if(str(offspring_a)!=offspring_a_save):
			offspring_a_mutation=str(offspring_a)
		offspring_b_mutation=""
		if(str(offspring_b)!=offspring_b_save):
			offspring_b_mutation=str(offspring_b)
		family_register+=offspring_a_save + " parents " + str(parents[0]) + "&" +  str(parents[1]) + " mutation " + str(offspring_a_mutation) + "\n"
		family_register+=offspring_b_save + " parents " + str(parents[0]) + "&" +  str(parents[1]) + " mutation " + str(offspring_b_mutation) + "\n"

		#add offspring to generation
		next_generation += [offspring_a, offspring_b]
	population = next_generation

	#reduce overlap by symmetry (no difference in 7-1 and 7-0 or small anchor)
	for i in range(len(population)):
		for j in range(len(population[i])-1):
			if(population[i][j] == 7):
				population[i][j+1] = 0
		population[i][0] = 0

	#find unique individuals
	individuals = list()
	for i in range(len(population)):
		if((population[i] in individuals)==False):
			individuals.append(population[i])
	unique_individuals = len(individuals)
	print("individuals " + str(individuals))

	#fill rest of generation randomly
	missing_individuals = population_size-unique_individuals
	individuals_to_add = populate_func()[0:missing_individuals]
	print("unique_individuals " + str(unique_individuals))
	family_register += "unique_individuals " + str(unique_individuals) + "\n"
	print(len(individuals_to_add))
	population = individuals
	population+=individuals_to_add

	#invoke evaluation of  new population
	population_for_fitness_eval = copy.deepcopy(population)
	for i in range(0, population_size):
		fitness_func(population_for_fitness_eval[i],generation,i)
	
	return population, family_register

