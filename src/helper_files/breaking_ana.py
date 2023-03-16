import os
import numpy as np
import tmoutproc as top
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from more_itertools import sort_together
from scipy.optimize import curve_fit
import sys

__har2eV__ = 27.2114
__eV2J__ = 1.60218e-19
__bohr2Ang__ = 0.529177
F_breaking = 1.5E-9


def calculate_s_s_dist(path, filename, write_out=False):
    """
    Calculates distance of Sulfur atoms (fixed) in Angstroms for stretching procedure and stores it under {gen}_{ind}_s_s_dist.dat
    path: where disp_pos and disp_neg is located
    filename: prefix for output file
    """

    s_s_dist = list()
    dirs = list()
    rootdir = f'{path}/disp_neg'
    subdirs = [x[0] for x in os.walk(rootdir) if "transport" not in x[0]]
    subdirs = sorted(subdirs)[1:len(subdirs)]
    subdirs = subdirs[::-1]
    for j in range(0,len(subdirs)-1):
        dirs.append(subdirs[j])
        coord = top.read_coord_file(subdirs[j] + "/coord")
        lower = -1
        upper = -1
        for i in range(0,len(coord)):
            #if(len(coord[i]) == 5 and coord[i][4] == 'f'):
            if (len(coord[i]) == 4 and coord[i][3] == 's'):
                if(lower==-1):
                    lower = float(coord[i][2])
                else:
                    upper = float(coord[i][2])
                    s_s_dist.append(np.abs(upper-lower)/1.889725989)
                    break
    len_disp_neg = len(subdirs)
    assert len(s_s_dist)==len_disp_neg-1, "Not all anchors found"

    rootdir = f'{path}/disp_pos'
    subdirs = [x[0] for x in os.walk(rootdir) if "transport" not in x[0]]
    subdirs = sorted(subdirs)[1:len(subdirs)]
    found = list()
    for j in range(0,len(subdirs)):
        dirs.append(subdirs[j])
        coord = top.read_coord_file(subdirs[j] + "/coord")
        lower = -1
        upper = -1
        for i in range(0,len(coord)):
            #if(len(coord[i]) == 5 and coord[i][4] == 'f'):
            if (len(coord[i]) == 4 and coord[i][3] == 's'):
                if(lower==-1):
                    lower = float(coord[i][2])
                else:
                    upper = float(coord[i][2])
                    s_s_dist.append(np.abs(upper-lower)/1.889725989)
                    found.append(subdirs[j])
                    break

    assert len(s_s_dist) == len_disp_neg + len(subdirs) - 1, "Not all anchors found"
    if(write_out == True):
        top.write_plot_data(f"{path}/{filename}_s_s_dist.dat", (s_s_dist, s_s_dist), "dirs, s-s distance (ang)")
    return np.array(s_s_dist)


def calc_junction_extent(path, filename, write_out=False):
    """
    Calculates maximum extent of junction in z direction
    Args:
        path: path to calculation (where disp_pos and disp_neg is located)
        filename: prefix for output file
        write_out: Write_out of displacement data

    Returns:

    """
    s_s_dist = list()
    dirs = list()
    rootdir = f'{path}/disp_neg'
    subdirs = [x[0] for x in os.walk(rootdir) if "transport" not in x[0]]
    subdirs = sorted(subdirs)[1:len(subdirs)]
    subdirs = subdirs[::-1]

    for j in range(0,len(subdirs)-1):
        dirs.append(subdirs[j])
        coord = top.read_coord_file(subdirs[j] + "/coord")
        z_coord = list()
        for i in range(0,coord.shape[1]):
            if(coord[3,i] != 'h'):
                z_coord.append(coord[2,i])
        s_s_dist.append(np.abs(np.max(z_coord) - np.min(z_coord))*__bohr2Ang__)

    rootdir = f'{path}/disp_pos'
    subdirs = [x[0] for x in os.walk(rootdir) if "transport" not in x[0]]
    subdirs = sorted(subdirs)[1:len(subdirs)]

    for j in range(0,len(subdirs)):
        dirs.append(subdirs[j])
        coord = top.read_coord_file(subdirs[j] + "/coord")
        z_coord = list()
        for i in range(0,coord.shape[1]):
            if(coord[3,i] != 'h'):
                z_coord.append(coord[2,i])
        s_s_dist.append(np.abs(np.max(z_coord) - np.min(z_coord))*__bohr2Ang__)

    if(write_out == True):
        top.write_plot_data(f"{path}/{filename}_junction_extent.dat", (s_s_dist, s_s_dist), "dirs, s-s distance (ang)")
    return np.array(s_s_dist)



def breaking_ana(path, filename):
    """
    Calculates the force generated by the molecule. Therefore, the maximum extend of the junction is normalized to get the Displacement data. The total energy is taken from a prepared file (containing total energy from dft calculation)
    Args:
        path: path to calculation (where disp_pos and disp_neg is located)
        filename: Prefix for filename

    Returns:

    """


    #load energy and sort according to stretching step
    energy_data = top.read_plot_data(f"{path}/{filename}_totalEnergy.dat")
    stretching_step, energy = sort_together([energy_data[0, :], energy_data[1, :]])
    stretching_step = np.asarray(stretching_step, int)
    energy = np.asarray(energy, float)
    #remove doulbe entry for 0:
    #check if there are two entries for zero
    tmp = sorted(np.abs(stretching_step))
    if(tmp[0] == 0 and tmp[1] == 0):
        stretching_step = list(stretching_step)
        energy = list(energy)
        zero_index_to_remove = np.argmin(np.abs(stretching_step))
        energy.pop(zero_index_to_remove)
        stretching_step.pop(zero_index_to_remove)
        stretching_step = np.asarray(stretching_step)
        energy = np.asarray(energy)



    #load s_s dist
    #s_s_dist = calculate_s_s_dist(path, gen, ind)
    s_s_dist = calc_junction_extent(path, filename, write_out=True)

    #find energy_minimum, shift energy and convert to
    min_energy_index = np.argmin(energy)
    min_disp_index = np.argmin(np.abs(stretching_step))
    energy = energy - energy[min_energy_index]
    energy = energy*__har2eV__*__eV2J__
    #normalize s_s_dist and convert to si
    s_s_dist = s_s_dist * 1E-10
    shift = s_s_dist[min_energy_index]
    s_s_dist = s_s_dist - shift


    #fit to harmonic approx
    def harmonic_approximation(x, k):
        return 0.5 * k * x ** 2

    try:
        popt_symmetric, pcov_symmetric = curve_fit(harmonic_approximation, s_s_dist, energy)
        popt_asymmetric, pcov_asymmetric = curve_fit(harmonic_approximation, s_s_dist[min_energy_index:], energy[min_energy_index:])
    except RuntimeError:
        print("Fit to wrong model")
    energy_fitted_symmetric  = harmonic_approximation(s_s_dist, popt_symmetric)
    energy_fitted_asymmetric = harmonic_approximation(s_s_dist[min_energy_index:], popt_asymmetric)

    breaking_dist_symmetric = F_breaking / (popt_symmetric)
    breaking_dist_asymmetric = F_breaking / (popt_asymmetric)
    print(f"Breaking dist sym = {breaking_dist_symmetric-s_s_dist[min_disp_index]}")
    print(f"Force constant sym (N/m) = {popt_symmetric}")
    print(f"Breaking dist asym = {breaking_dist_asymmetric-s_s_dist[min_disp_index]}")
    print(f"Force constant sym (N/m) = {popt_asymmetric}")

    with open(f"{path}/{filename}_breaking_ana.dat", 'w') as f:
        f.write(f"Breaking dist sym = {breaking_dist_symmetric}\n")
        f.write(f"Force constant sym (N/m) = {popt_symmetric}\n")
        f.write(f"Breaking dist asym = {breaking_dist_asymmetric}\n")
        f.write(f"Force constant asym (N/m) = {popt_asymmetric}\n")

    numerical_force = np.gradient(energy, np.mean(np.diff(s_s_dist)))

    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot((s_s_dist-s_s_dist[min_disp_index])*1E10, energy, label="data", color="g", lw=0.1, marker="x")
    ax1.plot((s_s_dist-s_s_dist[min_disp_index])*1E10, energy_fitted_symmetric, label="fit (sym)", color="g", ls="--")

    ax1.plot((s_s_dist[min_energy_index:] - s_s_dist[min_disp_index]) * 1E10, energy_fitted_asymmetric, label="fit (asym)", color="orange", ls="--")

    ax1.set_xlabel(r"Displacement ($\AA$)", fontsize=15)
    ax1.set_ylabel("Energy (J)", fontsize=15)
    ax1.tick_params(axis="x", labelsize=12)
    ax1.tick_params(axis="y", labelsize=12)
    ax1.grid()


    ax2.plot((s_s_dist-s_s_dist[min_disp_index])*1E10, numerical_force*1E9, label="data", color="g")
    ax2.plot((s_s_dist-s_s_dist[min_disp_index])*1E10, popt_symmetric * s_s_dist * 1E9, label="fit (sym)", color="g", ls="--")
    ax2.plot((s_s_dist[min_energy_index:] - s_s_dist[min_disp_index]) * 1E10, popt_asymmetric * s_s_dist[min_energy_index:] * 1E9, label="fit (asym)", color="orange",
             ls="--")
    ax2.axhline(F_breaking*1E9, color="black", ls="--", label=r"$\mathrm{F_{max}^{Au-S}}$")

    ax2.legend(fontsize = 12)
    ax2.grid()
    ax2.set_xlabel(r"Displacement ($\AA$)", fontsize=15)
    ax2.set_ylabel("Force (nN)", fontsize=15)
    ax2.tick_params(axis="x", labelsize=12)
    ax2.tick_params(axis="y", labelsize=12)

    plt.tight_layout()
    plt.savefig(f"{path}/{filename}_breaking_ana.pdf", bbox_inches='tight')
    plt.savefig(f"{path}/{filename}_breaking_ana.svg", bbox_inches='tight')

    #plt.show()


if __name__ == '__main__':
    path = sys.argv[1]
    filename = sys.argv[2]

    breaking_ana(path, filename)