import numpy as np
from random import choices, randint, randrange, random
from typing import List, Callable, Tuple
from collections import namedtuple
import tmoutproc as top
import configparser
import os
import os.path
from os import path
import subprocess
import shutil

Genome = List[int]
Point = namedtuple("Point", ['x','y', 'z'])
Angle = float
Building_Block = namedtuple('Building_Block', ['abbrev', 'num_atoms','origin', 'para_pos','para_angle', 'meta_pos','meta_angle', 'ortho_pos','ortho_angle','fixed_left','complexity', 'path'])
Coupling = namedtuple('Coupling', ['abbrev'])


def process_block_to_add(coupling_point: Point, coupling_angle :  Angle, conjugation_angle : Angle, cc_bond_length:float, block_to_add: Building_Block):
	"""
	Adds block_to_add to left_block. Method takes care of right alignment, shifting and rotation of block_to_add.

	Args:
		param1 (Point): coupling point
		param2 (Angle): coupling angle
		param3 (Angle): conjugation angle
		param4 (float): c-c bond length (Angstrom)
		param5 (Building_Block): block to be added

	Returns:
		np.ndarray: ([atom/x/y/z, line])

	"""

	#load data for block_to_add
	coord_xyz = top.read_xyz_file(block_to_add.path)
	coord_xyz_save = np.copy(coord_xyz)

	

	#rotate around z axis for conjugation, rotate around y for right orientation,shift to coupling point and shift by c-c bond length in right direction
	for i in range(0,coord_xyz.shape[1]):
		
		#conjugation x' = cos(phi)*x-sin(phi)*y
		coord_xyz[1,i] = coord_xyz_save[1,i]*np.cos(conjugation_angle)-np.sin(conjugation_angle)*coord_xyz_save[2,i]
		#conjugation y' = sin(phi)*x+cos(phi)*y
		coord_xyz[2,i] = coord_xyz_save[1,i]*np.sin(conjugation_angle)+np.cos(conjugation_angle)*coord_xyz_save[2,i]

		coord_xyz_save = np.copy(coord_xyz)

		#rotation around y: x' = cos(phi)*x+sin(phi)*z
		coord_xyz[1,i] = coord_xyz_save[1,i]*np.cos(coupling_angle)+np.sin(coupling_angle)*coord_xyz_save[3,i]
		#rotation around y: z' = cos(phi)*z-sin(phi)*x
		coord_xyz[3,i] = -coord_xyz_save[1,i]*np.sin(coupling_angle)+np.cos(coupling_angle)*coord_xyz_save[3,i]

		#shift to coupling point x -> x+x_c
		coord_xyz[1,i] = coord_xyz[1,i]+coupling_point.x
		#shift to coupling point y -> y+y_c
		coord_xyz[2,i] = coord_xyz[2,i]+coupling_point.y
		#shift to coupling point z -> z+z_c
		coord_xyz[3,i] = coord_xyz[3,i]+coupling_point.z

		#shift by C-C bond length in e_c direction sin=-0.866 cos=-0.499
		coord_xyz[1,i] = coord_xyz[1,i]+cc_bond_length*np.sin(coupling_angle)
		#shift by C-C bond length in e_c direction
		coord_xyz[3,i] = coord_xyz[3,i]+cc_bond_length*np.cos(coupling_angle)
	return coord_xyz




def construction_loop(genome : Genome, building_blocks, config_path, xyz_file_path):
	"""
	Construction loop Genome -> proper xyz file

	Args:
		param1 (Genome): Genome to build

	Returns:
		

	"""
	def determine_coupling_index(genome: Genome, index:int, building_blocks=building_blocks):
		"""
		determines coupling index (atom and corresponding line in xyz file of building block refered in genome[index]) and coupling angle

		Args:
			param1 (Genome): Genome to build
			param2 (int): index which block is processed and used as coupling point. Must be even -> odd indices are couplings

		Returns:
			(int,float): corresponding line in xyz file of building block refered in genome[index], coupling angle

		"""

		if(index > len(genome)-2 or index < 0):
			raise ValueError("index is out of proper range")

		# coupling after building_block of interest 
		i = index + 1
		
		#para
		if(genome[i]==0):
			coupling_index = building_blocks[genome[index]].para_pos
			coupling_angle = building_blocks[genome[index]].para_angle
		#meta
		elif(genome[i]==1):
			coupling_index = building_blocks[genome[index]].meta_pos
			coupling_angle = building_blocks[genome[index]].meta_angle
		#ortho
		elif(genome[i]==2):
			coupling_index = building_blocks[genome[index]].ortho_pos
			coupling_angle = building_blocks[genome[index]].ortho_angle
		else:
			raise ValueError("coupling seems to be funny")
		return coupling_index, coupling_angle

	def write_file_parts_to_file(xyz_file_parts, path, fixed_beginning, fixed_end,complexity, config_path):
		"""
		write xyz file parts to proper xyz file and turbomole coord file. Complexity is written to file 

		Args:
			param1 (List of np.ndarray): List of xyz files
			param2 (String): path
			param3 (int): fixed_beginning (index of atom in first block which should be fixed)
			param4 (int): fixed_end (index of atom in last block which should be fixed)
			param5 (int): complexity of whole molecule
			param6 (String): path to config file
		Returns:
			

		"""
		#load ang to bohr factor
		cfg = configparser.ConfigParser()
		cfg.read(config_path, encoding='utf-8')

		#write complexity to file
		with open(path+"/complexity", "w") as file_complexity:
			file_complexity.write(str(complexity))
		file_complexity.close()

		concat_xyz = np.concatenate(xyz_file_parts, axis=1)
		top.write_xyz_file(path+"/coord.xyz", concat_xyz)
		coord = top.x2t(concat_xyz)
		#fix right atoms
		coord[4,fixed_beginning] = "f"
		fixed_end = sum(np.array([xyz_file_parts[i].shape[1] for i in range(0,len(xyz_file_parts)-1)]))+fixed_end
		coord[4, fixed_end] = "f"
		top.write_coord_file(path+"/coord", coord)

		lower_limit = np.min(concat_xyz[3,:]) + 0.1
		upper_limit = np.max(concat_xyz[3, :]) - 0.1
		with open(path+"/limits", "w") as limits:
			limits.write(str(lower_limit) + "\n")
			limits.write(str(upper_limit))


	def determine_nearest_neighbor(datContent, coupling_index, atom_type):
		"""
		determines nearest neighbor of atom with index coupling index in dat content of atom type atom_type

		Args:
			param1 (List of np.ndarray): List of xyz files
			param2 (int): coupling_inxex
			param3 (string): atom_type of nearest neighbour
		Returns:
			int : index of nearest neighbour
		"""
		intersting_atoms = list()
		intersting_atoms_distance = list()
		for i in range(0, len(datContent[1,:])):
			if(datContent[0,i]==atom_type):
				intersting_atoms.append(i)
				distance = (float(datContent[1,i])-float(datContent[1,coupling_index]))**2+(float(datContent[2,i])-float(datContent[2,coupling_index]))**2+(float(datContent[3,i])-float(datContent[3,coupling_index]))**2			
				intersting_atoms_distance.append(distance)
		intersting_atoms = [x for _,x in sorted(zip(intersting_atoms_distance,intersting_atoms))]

		return intersting_atoms[0]

	def align_z_along_fixed_ends(xyz_file_parts, fixed_beginning, fixed_end):
		"""
		Align molecule z axis along fixed ends. This is done by rotation about the axis given by curl(vec(fixed_beginning->fixed_end), e_z) by the angle between vec(fixed_beginning-fixed_end) and e_z

		Args:
			param1 (List of np.ndarray): List of xyz files
			param2 (int): index in xyz_file_parts[0] of fixed beginning
			param3 (int): index in xyz_file_parts[-1] of fixed end
		Returns:
			int : (List of np.ndarray): List of xyz file
		"""

		molecule_axis = [xyz_file_parts[-1][1,fixed_end],xyz_file_parts[-1][2,fixed_end],xyz_file_parts[-1][3,fixed_end]]


		angle = np.arccos(molecule_axis[2]/np.linalg.norm(molecule_axis))
		theta = angle

		if(angle != 0):
			#calculate rotation axis
			rotation_axis = np.cross(molecule_axis, [0.0,0.0,1.0])
			rotation_axis = 1.0/np.linalg.norm(rotation_axis)*rotation_axis
			u = rotation_axis

			#calculate rotation_matrix
			rotation_matrix = [[np.cos(theta) + u[0]**2 * (1-np.cos(theta)), u[0] * u[1] * (1-np.cos(theta)) - u[2] * np.sin(theta), u[0] * u[2] * (1 - np.cos(theta)) + u[1] * np.sin(theta)],
	            [u[0] * u[1] * (1-np.cos(theta)) + u[2] * np.sin(theta), np.cos(theta) + u[1]**2 * (1-np.cos(theta)), u[1] * u[2] * (1 - np.cos(theta)) - u[0] * np.sin(theta)],
	            [u[0] * u[2] * (1-np.cos(theta)) - u[1] * np.sin(theta), u[1] * u[2] * (1-np.cos(theta)) + u[0] * np.sin(theta), np.cos(theta) + u[2]**2 * (1-np.cos(theta))]]

			for j in range(0, len(xyz_file_parts)):
				for i in range(0, len(xyz_file_parts[j][1,:])):
					 
					vector_to_rotate = [round(float(xyz_file_parts[j][1,i]),5),round(float(xyz_file_parts[j][2,i]),5),round(float(xyz_file_parts[j][3,i]),5)]
					rotated_vector = np.asmatrix(rotation_matrix)*np.asmatrix(vector_to_rotate).T
					xyz_file_parts[j][1,i] = round(rotated_vector[0,0],5)
					xyz_file_parts[j][2,i] = round(rotated_vector[1,0],5)
					xyz_file_parts[j][3,i] = round(rotated_vector[2,0],5)
			return xyz_file_parts
		else:
			return xyz_file_parts



	#load properties from config file 
	cfg = configparser.ConfigParser()
	cfg.read(config_path, encoding='utf-8')
	cc_bond_length = float(cfg.get('Building Procedure', 'CC_bond_lengt'))
	conjugation_angle_from_file = float(cfg.get('Building Procedure', 'conjugation_angle'))
	building_block_path = cfg.get('Building Procedure', 'building_block_path')

	#ensure that genome is not empty
	if(len(genome) < 1):
		print("Genome was emtpy")
		# TODO: proper treatment


	#add anchor to end -> couplings are missing 
	#add left anchor
	anchor_left, anchor_right = load_anchors_blocks(building_block_path)
	building_blocks.append(anchor_left)
	#para coupling
	genome.insert(0, len(building_blocks)-1)
	#add right anchor
	building_blocks.append(anchor_right)
	#para coupling
	genome.append(len(building_blocks)-1)

	#data content of every part of xyz file is stored in this list	
	xyz_file_parts = list()

	#first block as initialization directly added to list
	coupling_point = Point(x=0.0, y=0.0, z=0.0)
	coupling_angle = 0.0
	coupling_index = -1
	conjugation_angle = 0
	additional_angle = 0.0

	#indices for fixed atoms in beginning and end of chain
	fixed_beginning = 0
	fixed_end = 0

	#complexity measure of molecule
	complexity = 0
	for i in range(0, len(genome)):
		complexity += building_blocks[genome[i]].complexity
		#odd index -> coupling
		if(i%2==1):		
			#conclude coupling point
			x_c = float(xyz_file_parts[-1][1,coupling_index])
			y_c = float(xyz_file_parts[-1][2,coupling_index])
			z_c = float(xyz_file_parts[-1][3,coupling_index])
			coupling_point = Point(x=x_c, y=y_c, z=z_c)		


		#even index -> building block
		elif(i%2 == 0):			

			#handle rotation to process consecutive para or ortho couplings
			additional_angle += (-1)**(i/2+1)*np.pi
			additional_angle = 0
			
			#first block must not be shifted	
			if(i == 0):
				datContent = process_block_to_add(coupling_point, coupling_angle, conjugation_angle+additional_angle, 0.0, building_blocks[genome[i]])
				fixed_beginning = building_blocks[genome[i]].fixed_left
				if(building_blocks[genome[i]].fixed_left == -1):
					print("Error in first block: fixed atom not properly specified")
			else:
				datContent = process_block_to_add(coupling_point, coupling_angle, conjugation_angle+additional_angle, cc_bond_length, building_blocks[genome[i]])
				#find fix index of last block
				if(i == len(genome)-1):
					#para_pos is assumed to be right coupling point
					fixed_end = building_blocks[genome[i]].para_pos
					if(building_blocks[genome[i]].para_pos == -1):
						print("Error in last block: fixed atom not properly specified")



			#determine index of atom at origin
			origin = building_blocks[genome[i]].origin

			#if other block will be added -> hydrogen at c coupling atom must be removed
			if(i != len(genome)-1):				
				#determine coupling index and coupling angle
				coupling_index, coupling_angle_single = determine_coupling_index(genome,i,building_blocks)

				#handle sign to process consecutive para or ortho couplings
				#coupling_angle += (coupling_angle_single*(-1)**(i/2+1))
				coupling_angle += (coupling_angle_single)
				

				#remove hydrogen or other atom bonded to coupling atom
				nearest_neighbour = determine_nearest_neighbor(datContent, coupling_index, "H")
				datContent = np.delete(datContent,nearest_neighbour,1)
				
				#update coupling index and fixed beginning
				if(coupling_index>nearest_neighbour):				
					coupling_index -= 1
					if(i == 0 and fixed_beginning>nearest_neighbour):
						fixed_beginning -=1
				#update origin
				if(origin>nearest_neighbour):
					origin -=1


			#hydrogen bonded to C atom at origin must be removed, too (except for first atom)
			if(i != 0):						
				
				#remove hydrogen or other atom bonded to atom at origin
				nearest_neighbour = determine_nearest_neighbor(datContent, origin, "H")
				datContent = np.delete(datContent,nearest_neighbour,1)
				#update coupling index and fixed ending
				if(coupling_index>nearest_neighbour):					
					coupling_index = coupling_index -1
					if(i == len(genome)-1 and fixed_end>nearest_neighbour):
						fixed_end -=1
					pass

			xyz_file_parts.append(datContent)

			#alternating conjugation
			#conjugation_angle += (-1)**(i/2+1)*conjugation_angle_from_file
			conjugation_angle -= conjugation_angle_from_file

	#align molecule axis to z
	xyz_file_parts= align_z_along_fixed_ends(xyz_file_parts, fixed_beginning, fixed_end)

	#write xyz_file_parts to xyz file
	write_file_parts_to_file(xyz_file_parts, xyz_file_path, fixed_beginning, fixed_end, complexity, config_path)	

		
def load_building_blocks(path):
	"""
	load building blocks and set up Building_Block objects

	Args:
		param1 (path): path to dir where building_blocks are located
	Returns:
		list(Building_Block)
	"""		
	#TODO : automatization
	benzene = Building_Block(abbrev="B", num_atoms=6,origin=0, para_pos=3, para_angle=0, meta_pos=4 , meta_angle = -np.pi/3., ortho_pos=5, ortho_angle=-2.*np.pi/3, fixed_left = -1,complexity=1, path=path+"/benzene.xyz")
	napthtalene = Building_Block(abbrev="N", num_atoms=18,origin=0, para_pos=12, para_angle=0., meta_pos=11 , meta_angle = -np.pi/3., ortho_pos=10, ortho_angle=-2.*np.pi/3, fixed_left = -1,complexity=1, path=path+"/naphtalene.xyz")
	dbPc1 = Building_Block(abbrev="dbPc1", num_atoms=32,origin=13, para_pos=1, para_angle=0, meta_pos=0 , meta_angle = +np.pi/3., ortho_pos=0, ortho_angle=-2.*np.pi/3, fixed_left = -1,complexity=1, path=path+"/dbPc1_block.xyz")
	dbPc4 = Building_Block(abbrev="dbPc4", num_atoms=55,origin=22, para_pos=1, para_angle=0, meta_pos=0 , meta_angle = -np.pi/3., ortho_pos=0, ortho_angle=-2.*np.pi/3, fixed_left = -1,complexity=1, path=path+"/dbPc4.xyz")
	dbPc6 = Building_Block(abbrev="dbPc6", num_atoms=52,origin=17, para_pos=0, para_angle=0, meta_pos=1 , meta_angle = -np.pi/3., ortho_pos=0, ortho_angle=-2.*np.pi/3, fixed_left = -1,complexity=1, path=path+"/dbPc6.xyz")
	dbPc5 = Building_Block(abbrev="dbPc5", num_atoms=58,origin=12, para_pos=26, para_angle=0, meta_pos=20 , meta_angle = -np.pi/3., ortho_pos=0, ortho_angle=-2.*np.pi/3, fixed_left = -1,complexity=1, path=path+"/dbPc5.xyz")
	pseudo_para_naph_PCP = Building_Block(abbrev="pseudo-para_naph_PCP", num_atoms=44,origin=0, para_pos=18, para_angle=0, meta_pos=16 , meta_angle = -np.pi/3, ortho_pos=0, ortho_angle=-2.*np.pi/3, fixed_left = -1,complexity=1, path=path+"/pseudo-para_naph_PCP.xyz")
	line =Building_Block(abbrev="line", num_atoms=4,origin=0, para_pos=1, para_angle=0, meta_pos=1 , meta_angle = 0., ortho_pos=0, ortho_angle=-2.*np.pi/3, fixed_left = -1,complexity=1, path=path+"/line.xyz")
	#rot=Building_Block(abbrev="line", num_atoms=47,origin=6, para_pos=16, para_angle=0, meta_pos=20 , meta_angle = 0., ortho_pos=0, ortho_angle=-2.*np.pi/3, fixed_left = -1,complexity=2, path=path+"/rot.xyz")
	#stacked_anth=Building_Block(abbrev="stacked_anth", num_atoms=62,origin=3, para_pos=22, para_angle=0, meta_pos=30 , meta_angle = 0., ortho_pos=0, ortho_angle=-2.*np.pi/3, fixed_left = -1,complexity=2, path=path+"/stacked_anth.xyz")
	
	building_blocks = [benzene,napthtalene,dbPc1,dbPc4,dbPc6, dbPc5,pseudo_para_naph_PCP, line]

	return building_blocks

def load_anchors_blocks(path):
	"""
	load anchor blocks and set up Building_Block objects.

	Args:
		param1 (path): path to dir where anchors are located
	Returns:
		list(Building_Block)
	"""		
	#TODO : automatization
	left = Building_Block(abbrev="l", num_atoms=2,origin=0, para_pos=0, para_angle=0, meta_pos=0 , meta_angle = 0., ortho_pos=0, ortho_angle=-2.*np.pi/3, fixed_left = 0,complexity=1, path=path+"/anchor_small_left.xyz")
	right = Building_Block(abbrev="r", num_atoms=2,origin=0, para_pos=0, para_angle=0., meta_pos=0 , meta_angle = 0., ortho_pos=0, ortho_angle=-2.*np.pi/3, fixed_left = -1,complexity=1, path=path+"/anchor_small_right.xyz")
	
	anchors = [left,right]

	return anchors


		
def process_genome(generation : int, individual: int, genome:Genome, run_path):
	"""
	translates genome to xyz file. xyz file will be stored in $data/generation/individual and stretching and other calculations will be invoked. If geome has been processed in a previous generation, the data will be copied from the archive

	Args:
		param1 (int): generation
		param2 (int): individual in generation
		param3 (Genome): genome to process
		param4 (String): path of current run
	Returns:
		int : success (0), failure (-1)
	"""		
	#set up config path
	config_path = run_path + "/config"

	#check where building blocks are stored and generation data should be stored
	cfg = configparser.ConfigParser()
	cfg.read(config_path, encoding='utf-8')
	#TODO: set up correctly
	building_block_path = cfg.get('Building Procedure', 'building_block_path')
	generation_data_path = run_path + "/" + cfg.get('Building Procedure', 'generation_data_path')

	#create generation directory
	calc_path = generation_data_path + "/" + str(generation)
	try:
		#create generation dir
		if(path.exists(calc_path) == False):
			os.mkdir(calc_path)
	except OSError:
		print ("Creation of the directory %s failed" % calc_path)
		return -1

	#check if genome has been processed alreaddy 
	cfg = configparser.ConfigParser()
	cfg.read(config_path, encoding='utf-8')
	archive_path = cfg.get('Basics', 'archive_archive_path')
	print(archive_path)
	if(os.path.exists(archive_path)==False):
		print("No archive found")
	else:
		#read existing archinve
		archive_file = open(archive_path, "r")
		archive_population = list()
		archive_paths = list()
		for line in archive_file:
			if(len(line)<3):
				continue
			line = line.strip().split("	")
			tmp = line[0].replace("[", "").replace("]", "")
			tmp = tmp.split(",")
			tmp = [int(tmp[i]) for i in range(0,len(tmp))]
			archive_population.append(tmp)
			archive_paths.append(line[1])
		archive_file.close()

		#genome was procesed
		if (genome in archive_population)==True:
			print("copying existing genome")
			index = archive_population.index(genome)
			scr_dir = archive_paths[index] + "/."
			print("scr_dir " + str(scr_dir))
			dst_dir = generation_data_path + "/" + str(generation)+ "/" +str(individual) + "/"
			print("dst_dir " + str(dst_dir))
			if(path.exists(dst_dir) == True):
				print("Other job running ... Aborting")
				raise ValueError('Other job running ... Aborting')
			os.system("mkdir " + str(dst_dir))
			dst_dir += "."
			os.system("cp -R " + scr_dir + " " + dst_dir)
			#create DONE file
			DONE_file = generation_data_path + "/" + str(generation)+ "/" + str(generation)+ "_" +str(individual) + "_DONE"
			os.system("touch " + DONE_file)
			return 0




	#genome has not been calculated -> init calcs 
	#create directories for calculations  
	calc_path = generation_data_path + "/" + str(generation)
	try:
		#create generation dir
		if(path.exists(calc_path) == False):
			os.mkdir(calc_path)
		#create individual dir
		calc_path = generation_data_path + "/" + str(generation)+ "/" + str(individual)
		os.mkdir(calc_path)
	except OSError:
	    print ("Creation of the directory %s failed" % calc_path)
	    return -1
	
	#load building blocks
	building_blocks = load_building_blocks(building_block_path)	

	#construct molecule from genome
	construction_loop(genome, building_blocks, config_path, calc_path)

	#run next step -> invoke turbomole calculations
	set_up_turbo_calculations_path = cfg.get('Basics', 'helper_files') + "/set_up_turbo_calculations.sh"		
	os.system(set_up_turbo_calculations_path+" "+calc_path+" "+config_path + " " + str(generation) + " " + str(individual))




if __name__ == '__main__':
	process_genome(0,2,[0,1,0,2,1,4,0,5,0],"C:/Users/Matthias Blaschke/OneDrive - Universität Augsburg/Code/genetic_algorithm/test")



		