
from helper_files import genome_to_molecules as gtm

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
	"""
	Generates genome mit maximum lengt of max_lengt (-> size random). Genome is generated from available couplings and building blocks

	Args:
		param1 (int) : maximum length of 
                
                

	Returns:
		Genome
        """
	num_building_blocks = randrange(1,max_length)
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

def fitness(genome: Genome, generation: int, individual: int) -> int:
	
	gtm.process_genome(generation,individual,genome,"/alcc/gpfs2/home/u/blaschma/test/")




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