from random import choices, randint, randrange, random
from typing import List, Callable, Tuple
from collections import namedtuple
from functools import partial
import hamilton_helper as hahe
import numpy as np
import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt


Genome = List[int]
Population = List[Genome]
FitnessFunc = Callable[[Genome], int]
PopulateFunc = Callable[[], Population]
SelectionFunc = Callable[[Population, FitnessFunc],Tuple[Genome, Genome]]
CrossoverFunc = Callable[[Genome,Genome], Tuple[Genome, Genome]]
MutationFunc = Callable[[Genome], Genome]


Building_Block = namedtuple('Building_Block', ['abbrev', 'num_atoms', 'para_pos', 'meta_pos', 'ortho_pos', 'path'])
Coupling = namedtuple('Coupling', ['abbrev'])

benzene = Building_Block(abbrev="B", num_atoms=6, para_pos=3, meta_pos=4 ,ortho_pos=5, path="./hamiltionians/benzene.txt")
naphthalene = Building_Block(abbrev="N", num_atoms=10, para_pos=7, meta_pos=8 ,ortho_pos=9, path="./hamiltionians/naphthalene.txt")
anthracen = Building_Block(abbrev="A", num_atoms=14, para_pos=11, meta_pos=12 ,ortho_pos=13, path="./hamiltionians/anthracen.txt")
building_blocks = [benzene,naphthalene,anthracen]

para = Coupling(abbrev="p")
meta = Coupling(abbrev="m")
ortho = Coupling(abbrev="o")
couplings = [para, meta, ortho]

hamiltionians = 0


def generate_genome(max_length: int) -> Genome:
	num_building_blocks = randrange(2,max_length)
	indices_building_blocks = np.linspace(0,len(building_blocks),len(building_blocks),endpoint=False,dtype=int)   
	indices_couplings = np.linspace(0,len(couplings),len(couplings),endpoint=False,dtype=int)   
	selected_building_blocks = choices(indices_building_blocks, k=num_building_blocks)
	selected_couplings = choices(indices_couplings, k=num_building_blocks-1)

	genome=list()
	for i in range(0, num_building_blocks):
		genome.append(selected_building_blocks[i])
		if(i != num_building_blocks-1):
			genome.append(selected_couplings[i])
	return genome

	

def generate_population(size: int, genome_length: int) -> Population:
	return [generate_genome(genome_length) for _ in range(size)]

def fitness(genome: Genome) -> int:
	hamilton = genome_to_full_hamiltion(genome)
	eigenvalues, eigenvectors = np.linalg.eigh(hamilton)
	#print(eigenvalues)
	#eigenvalues = eigenvalues*-1
	length_pan = np.exp(-((float(len(genome))-3)**2/2.0))
	return eigenvalues[0]+length_pan

def genome_to_full_hamiltion(genome, t=1):
	hamiltonias = load_hamiltonians(building_blocks)
	length=0
	for i in range(0,len(genome)):
		if(i%2==0):
			length += building_blocks[genome[i]].num_atoms

	full_hamilton = np.zeros((length, length))
	handled_atoms = 0
	handled_coupling = 0
	t = -1
	for i in range(0, len(genome)):
		if(i%2==0):
			num_atoms = building_blocks[genome[i]].num_atoms
			#print(hamiltionians[genome[i]])
			full_hamilton[handled_atoms:handled_atoms+num_atoms,handled_atoms:handled_atoms+num_atoms] = hamiltonias[int(genome[i])]
			handled_atoms += num_atoms
		elif(i%2==1):
			#para
			if(genome[i]==0):
				coupling_index = building_blocks[genome[i-1]].para_pos
			#meta
			elif(genome[i]==1):
				coupling_index = building_blocks[genome[i-1]].meta_pos
			#ortho
			elif(genome[i]==2):
				coupling_index = building_blocks[genome[i-1]].ortho_pos
			else:
				raise ValueError("coupling seems to be funny")
			full_hamilton[handled_atoms-num_atoms+coupling_index, handled_atoms] = t
			full_hamilton[handled_atoms,handled_atoms-num_atoms+coupling_index] = t
	return full_hamilton


def load_hamiltonians(building_blocks):
	hamiltionians = list()
	for i in range(0, len(building_blocks)):
		data = hahe.load_hamilton(building_blocks[i].path)
		hamilton = hahe.replace_constants(data, e=5, t=-1,c=0)
		hamiltionians.append(hamilton)
	return hamiltionians


def selection_pair(population: Population, fitness_func: FitnessFunc) -> Population:
	return choices(
		population=population,
		weights=[fitness_func(genome) for genome in population],
		k=2
	)

def single_point_crossover(a:Genome, b:Genome) -> Tuple[Genome, Genome]:
	
	length_a = len(a)
	length_b = len(b)
	if(length_a<=1 or length_b<=1):
		return a, b
	#ensure that coupling and blocks alter
	cut_a = randrange(length_a)	
	if(cut_a%2 == 0):
		cut_b = randrange(int(length_b/2))
		cut_b = 2*cut_b
	else:
		cut_b = randrange(int(length_b/2))
		cut_b = 2*cut_b+1

	return a[0:cut_a] + b[cut_b:length_b], b[0:cut_b] + a[cut_a:length_a]

def mutation(genome: Genome, num: int=1, probability: float = 0.5) -> Genome:
	for _ in range(num):
		method = randrange(2)
		if(method == 0):
			return building_block_mutation(genome, probability)
		elif(method == 1):
			return coupling_mutation(genome, probability)
	return genome

def building_block_mutation(genome: Genome, probability: float = 0.5) -> Genome:	
	mutated_genome = list()
	if(random()<probability and len(genome)>=3):
		new_block = randrange(len(building_blocks))
		#print("new block " + str(new_block))
		block_to_mutate = randrange(int((len(genome)+1)/2))
		#print("block to mutate " + str(block_to_mutate))	
		block_to_mutate = block_to_mutate+block_to_mutate+1	
		
		
		mutated_genome.extend(genome[0:block_to_mutate-1])
		mutated_genome.append(new_block)
		mutated_genome.extend(genome[block_to_mutate:len(genome)])
		
		return mutated_genome
	return genome

def coupling_mutation(genome: Genome, probability: float = 0.5) -> Genome:	
	mutated_genome = list()
	if(random()<probability and len(genome)>=3):
		new_coupling = randrange(len(couplings))
		#print("new coupling " + str(new_coupling))		
		coupling_to_mutate = randrange(0,int((len(genome)-1)/2))
		#print("coupling to mutate " + str(coupling_to_mutate))			
		coupling_to_mutate = coupling_to_mutate+coupling_to_mutate+2
		
		#genome = genome[0:coupling_to_mutate-1] + couplings[new_coupling].abbrev + genome[coupling_to_mutate:len(genome)]
		
		mutated_genome.extend(genome[0:coupling_to_mutate-1])
		mutated_genome.append(new_coupling)
		mutated_genome.extend(genome[coupling_to_mutate:len(genome)])
		return mutated_genome
	return genome
#"""
def run_evolution(
		populate_func: PopulateFunc,
		fitness_func: FitnessFunc,
		fitness_limit: int,
		selection_func : SelectionFunc = selection_pair,
		crossover_func : CrossoverFunc = single_point_crossover,
		mutation_func: MutationFunc = mutation,
		generation_limit: int = 100
	) -> Tuple[Population, int]:
	fitness_value = list()
	population = populate_func()

	for i in range(generation_limit):
		print("generation " + str(i))
		population = sorted(
			population,
			key = lambda genome: fitness_func(genome),
			reverse = True
			)
		fitness_value.append(fitness_func(population[0]))
		if fitness_func(population[0]) >= fitness_limit:
			break

		next_generation = population[0:2]

		for j in range(int(len(population)/2)-1):
			parents = selection_func(population, fitness_func)
			offspring_a, offspring_b = crossover_func(parents[0], parents[1])
			offspring_a = mutation_func(offspring_a)
			offspring_b = mutation_func(offspring_b)

			next_generation += [offspring_a, offspring_b]
		population = next_generation

	population = sorted(
			population,
			key = lambda genome: fitness_func(genome),
			reverse = True
			)
	return population, i,fitness_value





population, generations, fitness_value = run_evolution(
	populate_func = partial(
			generate_population, size=20, genome_length=20
			),
	fitness_func = partial(
			fitness
			),
	fitness_limit = 1310,
	generation_limit = 25
	)

print(f"number of generations: {generations}")
print(f"list of things: {population[0]}")
plt.plot(fitness_value)
plt.show()

#"""
"""
for i in range(0,1000):
	genome_a = generate_genome(5)
	#genome_b = generate_genome(5)
	print(genome_a)
"""	

