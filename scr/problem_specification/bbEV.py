import configparser
from problem_specification.evaluation_methods import evaluation_methods
import numpy as np
from random import choices, randint, randrange, random
from collections import namedtuple
from typing import List, Callable, Tuple
from helper_files import genome_to_molecule as gtm
import copy

class bbEv(evaluation_methods.Evaluation):
	Genome = List[int]
	Population = List[Genome]
	FitnessFunc = Callable[[Genome, int, int], int]
	PopulateFunc = Callable[[], Population]
	SelectionFunc = Callable[[Population, FitnessFunc],Tuple[Genome, Genome]]
	CrossoverFunc = Callable[[Genome,Genome], Tuple[Genome, Genome]]
	MutationFunc = Callable[[Genome], Genome]

	def __init__(self, generation, individual, config_path):
		super().__init__(generation, individual)
		self.Building_Block = namedtuple('Building_Block', ['abbrev', 'num_atoms', 'para_pos', 'meta_pos', 'ortho_pos', 'path'])
		self.Coupling = namedtuple('Coupling', ['abbrev'])
		self.generation = generation
		self.individual = individual

		#self.benzene = Building_Block(abbrev="B", num_atoms=6, para_pos=3, meta_pos=4 ,ortho_pos=5, path="./hamiltionians/benzene.txt")
		#self.naphthalene = Building_Block(abbrev="N", num_atoms=10, para_pos=7, meta_pos=8 ,ortho_pos=9, path="./hamiltionians/naphthalene.txt")
		#self.anthracen = Building_Block(abbrev="A", num_atoms=14, para_pos=11, meta_pos=12 ,ortho_pos=13, path="./hamiltionians/anthracen.txt")
		#self.building_blocks = [benzene,naphthalene,anthracen]

		#todo : dependency on interchangeable module not good
		self.building_blocks=gtm.load_building_blocks("")

		self.para = self.Coupling(abbrev="p")
		self.meta = self.Coupling(abbrev="m")
		#self.ortho = self.Coupling(abbrev="o")
		self.couplings = [self.para, self.meta]

		self.config_path = config_path
		hamiltionians = 0


	def generate_genome(self, max_length: int) -> Genome:
		"""
		Generates genome mit maximum lengt of max_lengt (-> size random). Genome is generated from available couplings and building blocks

		Args:
			param1 (int) : maximum length of 
	                
	                

		Returns:
			Genome
	        """
	    #coupling must contain at least one block and two couplings
		#print("max_length " + str(max_length))
		if(max_length <= 1):
			print("sorry max_length too small")
			return -1
		num_building_blocks = randrange(1,max_length)
		
		indices_building_blocks = np.linspace(0,len(self.building_blocks),len(self.building_blocks),endpoint=False,dtype=int)   
		indices_couplings = np.linspace(0,len(self.couplings),len(self.couplings),endpoint=False,dtype=int)   
		selected_building_blocks = choices(indices_building_blocks, k=num_building_blocks)
		selected_couplings = choices(indices_couplings, k=num_building_blocks+2)

		genome=list()

		#add coupling to anchor
		genome.append(selected_couplings[0])
		for i in range(0, num_building_blocks):
			
			genome.append(selected_building_blocks[i])
			
			genome.append(selected_couplings[i+1])
			#print("genome " + str(genome))
		return genome

		

	def generate_population(self, size: int, genome_length: int) -> Population:
		return [self.generate_genome(genome_length) for _ in range(size)]

	def fitness(self, genome: Genome, generation: int, individual: int) -> float:
		genome_copy = copy.deepcopy(genome)
		gtm.process_genome(generation,individual,genome_copy,"/alcc/gpfs2/home/u/blaschma/genetic_run_2/")
		return random()




	def selection_pair(self, population: Population, fitness_value) -> Population:		
		return choices(
			population=population,
			weights=fitness_value,
			k=2
		)

	def crossover(self, a:Genome, b:Genome) -> Tuple[Genome, Genome]:
		#return (a,b)
		return self.single_point_crossover(a,b)

	def single_point_crossover_old(self, a:Genome, b:Genome) -> Tuple[Genome, Genome]:
		
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

	def single_point_crossover(self, a:Genome, b:Genome) -> Tuple[Genome, Genome]:
		
		length_a = len(a)
		length_b = len(b)
		minimim_length = np.min((length_a, length_b))
		if(length_a<=1 or length_b<=1):
			return a, b
			
		cut = randrange(minimim_length)			

		return a[0:cut] + b[cut:length_b], b[0:cut] + a[cut:length_a]

	def mutation(self, genome: Genome, num: int=2, probability: float = 0.5) -> Genome:		
		method = randrange(4)
		if(method == 0):
			return self.building_block_mutation(genome, probability)
		elif(method == 1):
			return self.coupling_mutation(genome, probability)
		elif(method == 2):
			return self.insert_mutation(genome, probability)
		elif(method == 3):
			return self.truncate_mutation(genome, probability)
		return genome

	def building_block_mutation(self, genome: Genome, probability: float = 0.5) -> Genome:	
		cfg = configparser.ConfigParser()
		cfg.read(self.config_path)
		probability = float(cfg.get('Genetic Algorithm', 'block_mutation_prob'))

		mutated_genome = list()
		if(random()<probability and len(genome)>=3):
			new_block = randrange(len(self.building_blocks))
			#print("new block " + str(new_block))
			block_to_mutate = randrange(int((len(genome)-1)/2))
			#print("block to mutate " + str(block_to_mutate))	
			block_to_mutate = block_to_mutate+block_to_mutate+2
			
			
			mutated_genome.extend(genome[0:block_to_mutate-1])
			mutated_genome.append(new_block)
			mutated_genome.extend(genome[block_to_mutate:len(genome)])
			print("block mutation!" + str(mutated_genome))
			return mutated_genome
		return genome

	def coupling_mutation(self, genome: Genome, probability: float = 0.5) -> Genome:	
		cfg = configparser.ConfigParser()
		cfg.read(self.config_path)
		probability = float(cfg.get('Genetic Algorithm', 'coupling_mutation_prob'))

		mutated_genome = list()
		if(random()<probability and len(genome)>=3):
			new_coupling = randrange(len(self.couplings))
			#print("new coupling " + str(new_coupling))		
			coupling_to_mutate = randrange(0,int((len(genome)+1)/2))
			#print("coupling to mutate " + str(coupling_to_mutate))			
			coupling_to_mutate = coupling_to_mutate+coupling_to_mutate+1
			
			#genome = genome[0:coupling_to_mutate-1] + couplings[new_coupling].abbrev + genome[coupling_to_mutate:len(genome)]
			
			mutated_genome.extend(genome[0:coupling_to_mutate-1])
			mutated_genome.append(new_coupling)
			mutated_genome.extend(genome[coupling_to_mutate:len(genome)])
			print("coupling mutation! " + str(mutated_genome))
			return mutated_genome
		return genome

	def insert_mutation(self, genome: Genome, probability: float = 0.5):
		cfg = configparser.ConfigParser()
		cfg.read(self.config_path)
		probability = float(cfg.get('Genetic Algorithm', 'insert_mutation_prob'))
		genome_length = float(cfg.get('Genetic Algorithm', 'genome_length'))

		mutated_genome = list()
		if(random()<probability and len(genome) < genome_length):

			new_coupling = randrange(len(self.couplings))
			new_block = new_block = randrange(len(self.building_blocks))

			insert_coupling=randrange(0,int((len(genome)+1)/2))
			insert_coupling=insert_coupling+insert_coupling+1
			to_add_at_end = genome[insert_coupling:len(genome)]
			mutated_genome.extend(genome[0:insert_coupling])
			mutated_genome.append(new_block)
			mutated_genome.append(new_coupling)
			mutated_genome.extend(to_add_at_end)
			print("insert mutation! " + str(mutated_genome))
			return mutated_genome

		return genome

	def truncate_mutation(self, genome: Genome, probability: float = 0.5):
		cfg = configparser.ConfigParser()
		cfg.read(self.config_path)
		probability = float(cfg.get('Genetic Algorithm', 'truncate_mutation_prob'))
		genome_length = float(cfg.get('Genetic Algorithm', 'genome_length'))

		mutated_genome = list()
		if(random()<probability and len(genome) > 3):

			block_to_truncate = randrange(int((len(genome)-1)/2))	
			block_to_truncate = block_to_truncate+block_to_truncate+2

			mutated_genome.extend(genome[0:block_to_truncate-1])
			mutated_genome.extend(genome[block_to_truncate+1:len(genome)])
			print("truncate mutation!" + str(mutated_genome))
			return mutated_genome

		return genome