import sys
import tmoutproc as top

"""
This script can be used to strech coord file for stretching trajectories.
@author: Matthias Blaschke
"""

def strech_on_z_coord(coord, disp, lower, upper):
    """
    Displaces atoms with z_coord < lower -disp/2 and atoms with z_coord>lower +disp/s. All coordinates between lower
    and upper are proportionally stretched, like on a rubber band

    Args:
        coord (np.ndarray): Coord file from top
        disp (float): displacement in Bohr
        lower (float): Lower cut in Bohr
        upper (float): Upper cut in Bohr

    Returns:
        coord (np.ndarray): Stretched coord file
    """
    assert lower < upper, "Limits are not reasonable"
    for i in range(0, coord.shape[1]):
        if (coord[2,i] < lower):
            coord[2,i] -= disp / 2
        elif (coord[2,i] > upper):
            coord[2,i] += disp / 2
        else:  # stretch all atom coordinates in between limits
            norm = (2 * coord[2,i] - upper - lower) / (upper - lower)
            coord[2,i] += norm * disp / 2
    return coord

def stretch_on_indices(coord, disp, lower_index, upper_index):
    """
    Displaces atoms with z_coord <= z_coord[lower_index] -disp/2 and atoms with z_coord >= z_coord[lower_index] +disp/s.
    All coordinates between lower and upper are proportionally stretched, like on a rubber band
    Args:
        coord (np.ndarray): Coord file from top
        disp (float): displacement in Bohr
        lower (int): Index of lower atom
        upper (int): Index of uppper atom

    Returns:
        coord (np.ndarray): Stretched coord file
    """
    lower = coord[2,lower_index]
    upper = coord[2, upper_index]
    assert  lower < upper, "Indices are not reasonable"
    for i in range(0, coord.shape[1]):
        old = coord[2,i]
        if (coord[2,i] <= lower):
            coord[2,i] -= disp / 2
        elif (coord[2,i] >= upper):
            coord[2,i] += disp / 2
        else:  # stretch all atom coordinates in between limits
            norm = (2 * coord[2,i] - upper - lower) / (upper - lower)
            coord[2,i] += norm * disp / 2
        print(coord[2,i], old, coord[2,i])
    return coord

def find_fixed_atoms(coord):
    """
    Finds the first two fixed atoms in coord file
    Args:
        coord: coord file from top

    Returns:
        fixed_atoms (list): List with indices of fixed atoms
    """
    fixed_atoms = []
    for i in range(0, coord.shape[1]):
        if(coord[4,i] == "f"):
            fixed_atoms.append(i)
        if(len(fixed_atoms)==2):
            break
    return fixed_atoms

def find_atom_index_by_type(coord, type):
    """
    Finds the first two atoms of type type in coord file
    Args:
        coord: coord file from top
        type (String): Atom type to be found

    Returns:
        fixed_atoms (list): List with indices of fixed atoms
    """
    atoms = []
    for i in range(0, coord.shape[1]):
        if(coord[3,i] == type):
            atoms.append(i)
        if(len(atoms)==2):
            break
    return atoms

if __name__ == '__main__':
    coord_in_file_path = sys.argv[2]
    coord_out_file_path = sys.argv[3]

    disp = float(sys.argv[4]) * top.ANG2BOHR
    coord_in_file = top.read_coord_file(coord_in_file_path)

    #cuts are defined by zcoord
    if(sys.argv[1] == "-zcoord"):
        lower = float(sys.argv[5]) * top.ANG2BOHR
        upper = float(sys.argv[6]) * top.ANG2BOHR
        coord = strech_on_z_coord(coord_in_file, disp, lower, upper)
    #cuts are defined by two fixed atoms. There can be only two fixed atoms. Mode is made for isolated molecules
    elif(sys.argv[1] == "-fixed_atoms"):
        fixed_atoms = find_fixed_atoms(coord_in_file)
        assert len(fixed_atoms)==2, "Not two fixed atoms"
        coord = stretch_on_indices(coord_in_file, disp, fixed_atoms[0], fixed_atoms[1])
    #cuts are defined by two atoms of type sys.argv[5]
    elif(sys.argv[1] == "-type"):
        type = str(sys.argv[5])
        atoms = find_atom_index_by_type(coord_in_file, type)
        assert len(atoms)==2, "Not two atoms"
        coord = stretch_on_indices(coord_in_file, disp, atoms[0], atoms[1])
    else:
        print("Usage")
        print("jobgen.py [-zcoord] coords_in coords_out displacement[Ang] lower[Ang] upper[Ang]")
        print("or")
        print("jobgen.py [-fixed_sulfur] coords_in coords_out displacement[Ang]")
        print("or")
        print("jobgen.py [-type] coords_in coords_out displacement[Ang] type")

    top.write_coord_file(coord_out_file_path, coord)