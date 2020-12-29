import numpy
import os
import sys
import shutil
import pickle

Ang2Bohr = 1.889725989
Bohr2Ang = 1/Ang2Bohr

template_dir = './templates/'
tmtrans_in = 'tmtrans_control.in'

def replace(source_path, to_replace, replace_with, destination_path = ''):
    temp_path = ''
    if source_path == destination_path or destination_path == '':
        temp_path = source_path + '.temp'
        shutil.copyfile(source_path, temp_path)
        destination_path = source_path
        source_path = temp_path
    source = open(source_path)
    destination = open(destination_path, 'w')

    if isinstance(replace_with, str):
        for line in source:
            if to_replace in line:
                line = line.replace(to_replace, replace_with)
            destination.write(line)
    elif isinstance(replace_with, list):
        for line in source:
            if to_replace in line: break
            else: destination.write(line)
        for line in replace_with:
            destination.write(line)
        for line in source:
            destination.write(line)
    else:
        print('replace_with has to be either a string or a list of strings. Exiting ...')
        exit()

    source.close()
    destination.close()
    if temp_path != '' and os.path.exists(temp_path): os.remove(temp_path)

def read_coords(filename):
    fi = open(filename, 'r') # input file
    coords = []
    elements = []
    fixed = []
    fi.readline()
    for line in fi:
        if('$' in line): break
        sline = line.split()
        coords.append(list(map(float, sline[0:3])))
        elements.append(sline[3])
        if(len(sline) == 4):
            fixed.append(False)
        elif(len(sline) == 5 and sline[4] in ['f', 'F']):
            fixed.append(True)
        else:
            print('Check format of coord file! Exiting ...')
            exit()
    return numpy.array(coords), elements, fixed

def write_coords(filename, coords, elements, fixed):
    fo = open(filename, 'w') # output file
    fo.write('$coord\n')
    for i, element in enumerate(elements):
        for d in range(3):
            fo.write(str(coords[i][d]) + ' ')
        fo.write(element + ' ' + ('f' if fixed[i] else '') + '\n')
    fo.write('$user-defined bonds\n')
    fo.write('$end\n')

if len(sys.argv) > 1:
    request = sys.argv[1]
    if request == '-displace':
        job_type = 'displace'
        if len(sys.argv) > 2:
            coord_in_file = sys.argv[2]
            coord_out_file = sys.argv[3]
            disp = float(sys.argv[4]) * Ang2Bohr
            lower = float(sys.argv[5]) * Ang2Bohr
            upper = float(sys.argv[6]) * Ang2Bohr
    else:
        print('Unknown job type.')
        exit()
else:
    print('Usage:')
    print('jobge.py -displace coords_in coords_out displacement[Ang] lower[Ang] upper[Ang]')
    exit()

##################################################################################################

# moves all atoms within [-inf,lower] by -disp/2
# and all atoms defined by [upper,+inf] by +disp/2
# all other positions are not changed
elastic = True
# elastic stretching means that all coordinates between lower and upper are
# proportionally stretched, like on a rubber band
if job_type == 'displace':
    if(lower > upper):
        print('Check lowe/upper limits! ' + str(lower) + ' ' + str(upper))
        print('Exiting ...')
        exit()
    coords, elements, fixed = read_coords(coord_in_file)
    for i, element in enumerate(elements):
        if(coords[i][2] < lower):
            coords[i][2] -= disp/2
        elif(coords[i][2] > upper):
            coords[i][2] += disp/2
        elif(elastic): # stretch all atom coordinates in between limits
            norm = (2*coords[i][2] - upper - lower) / (upper - lower)
            coords[i][2] += norm * disp/2
    write_coords(coord_out_file, coords, elements, fixed)
    #replace(template_dir + tmtrans_in,
            #'#LOWER_LIMIT#', str(lower*Bohr2Ang), destination_path = tmtrans_in)
    #replace(tmtrans_in, '#UPPER_LIMIT#', str(upper*Bohr2Ang))
