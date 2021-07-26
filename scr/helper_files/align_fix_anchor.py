import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from turbomoleOutputProcessing import turbomoleOutputProcessing as top
import os.path
from multiprocessing import Pool
from functools import partial
import sys
import configparser
from scipy.optimize import curve_fit





def find_fixed_atoms(coord):
	"""
	Finds index of fixed atoms in coord file stored in coord. index_left<index_right from building procedure

	Args:
		param1 (np.ndarray): coord file loaded with top.read_coord_file
	Returns:
		(index_left, index_right) 
	"""		
	index_1=-1
	index_2=-1
	for i in range(0,len(coord)):
		#if line has fixed symbol
		if(len(coord[i])==5):
			if(index_1==-1):
				index_1=i
			if(index_1!=-1):
				index_2=i
	if(index_1!=-1 and index_2!=-1):
		min=np.min((index_1,index_2))
		max=np.max((index_1,index_2))
		return(min,max)
	else:
		print("fixed atoms not foud")
		return (-1,-1)

def align_z_along_fixed_ends(coord, fixed_beginning, fixed_end):
	"""
	Align molecule z axis along fixed ends. This is done by rotation about the axis given by curl(vec(fixed_beginning->fixed_end), e_z) by the angle between vec(fixed_beginning-fixed_end) and e_z

	Args:
		param1 (List of np.ndarray): coord file loaded with top.load_coord_file
		param2 (int): index of fixed beginning (left)
		param3 (int): index fixed end (right)
	Returns:
		int : (List of np.ndarray): coord file
	"""

	#calculate vec(fixed_beginning->fixed_end)
	x_left = round(float(coord[fixed_beginning][0]),5)
	y_left = round(float(coord[fixed_beginning][1]),5)
	z_left = round(float(coord[fixed_beginning][2]),5)

	#shift to orgin
	for j in range(0, len(coord)):
		coord[j][0] = round(float(coord[j][0]),5)-x_left
		coord[j][1] = round(float(coord[j][1]),5)-y_left
		coord[j][2] = round(float(coord[j][2]),5)-z_left
	print("shift")
	print(coord)
	x_left = 0.0
	y_left = 0.0
	z_left = 0.0

	x_right = round(float(coord[fixed_end][0]),5)
	y_right = round(float(coord[fixed_end][1]),5)
	z_right = round(float(coord[fixed_end][2]),5)


	molecule_axis = [x_right-x_left, y_right-y_left,z_right-z_left]

	#calculate rotation angel

	angle = np.arccos(molecule_axis[2]/np.linalg.norm(molecule_axis))
	theta = angle
	
	if(angle != 0):
		#calculate rotation axis
		rotation_axis = np.cross(molecule_axis, [0.0,0.0,1.0])
		rotation_axis = 1.0/np.linalg.norm(rotation_axis)*rotation_axis
		u = rotation_axis
		#print("rotation axis " + str(rotation_axis))

		#calculate rotation_matrix
		rotation_matrix = [[np.cos(theta) + u[0]**2 * (1-np.cos(theta)), u[0] * u[1] * (1-np.cos(theta)) - u[2] * np.sin(theta), u[0] * u[2] * (1 - np.cos(theta)) + u[1] * np.sin(theta)],
            [u[0] * u[1] * (1-np.cos(theta)) + u[2] * np.sin(theta), np.cos(theta) + u[1]**2 * (1-np.cos(theta)), u[1] * u[2] * (1 - np.cos(theta)) - u[0] * np.sin(theta)],
            [u[0] * u[2] * (1-np.cos(theta)) - u[1] * np.sin(theta), u[1] * u[2] * (1-np.cos(theta)) + u[0] * np.sin(theta), np.cos(theta) + u[2]**2 * (1-np.cos(theta))]]



		for j in range(0, len(coord)):				 
			vector_to_rotate = [round(float(coord[j][0]),5),round(float(coord[j][1]),5),round(float(coord[j][2]),5)]
			rotated_vector = np.asmatrix(rotation_matrix)*np.asmatrix(vector_to_rotate).T

			coord[j][0] = round(rotated_vector[0,0],5)
			coord[j][1] = round(rotated_vector[1,0],5)
			coord[j][2] = round(rotated_vector[2,0],5)
		return coord
	else:
		return coord

def write_limits(path,coord,fixed_beginning, fixed_end, config_path):
	"""
	write limits for stretching. Done by extracting the z-coordinate of endings and adding 0.1Ang

	Args:
		param1 (string): path where limits file is sotred
		param2 (np.ndArray): coord file loaded with top.load_coord_file
		param3 (int): index of fixed beginning (left)
		param4 (int): index fixed end (right)
		param5 (string): path to config file 
	Returns:
		
	"""	

	#load ang to bohr factor
	cfg = configparser.ConfigParser()
	cfg.read(config_path)
	ang2Bohr = float(cfg.get('Building Procedure', 'ang2Bohr'))

	lower_limit = float(coord[fixed_beginning][2])/ang2Bohr
	upper_limit = float(coord[fixed_end][2])/ang2Bohr

	if(upper_limit==lower_limit):
		print("Faulty limits!")
		return
	if(upper_limit<lower_limit):
		upper_limit = a
		upper_limit = lower_limit
		lower_limit = a

	file = open(path+"/limits", "w")
	file.write(str(lower_limit) + "\n")
	file.write(str(upper_limit))




if __name__ == '__main__':
	
	#argv[1] : calc path: where coord and fixed file is located and where result is stored
	#argv[2] : config path
	#load coord file
	coord = top.read_coord_file(sys.argv[1]+ "/coord")
	#find achors
	fixed_index = find_fixed_atoms(coord)
	#rotate and shift
	coord = align_z_along_fixed_ends(coord, fixed_index[0], fixed_index[1])
	print(coord)
	#write coord filef
	top.write_coord_file(sys.argv[1]+ "/coord", coord)
	#write limits
	write_limits(sys.argv[1], coord, fixed_index[0], fixed_index[1], sys.argv[2])
