"""
geometry_analysis.py
This module contains the geometry analysis project from MolSSI workshop.
"""
import numpy 
import os
import sys

########################################################################
# HEPLER FUNCTIONS 
def open_xyz(input_file):
    """
    function to open and process an xyz coordinate file.
    Input : 
        input_file - name of the file to process
    Output :
        symbols - numpy array of chemical symbols
        coords - numpy array of 3D coordinates
    """
    # since we assume xyz, skip the first two lines
    data = numpy.genfromtxt(fname=input_file, dtype='unicode', skip_header=2)
    symbols = data[:,0]
    coords = data[:,1:].astype(numpy.float)
    return symbols, coords

# add a default value for min and max distance
def bond_check(distance, min_distance=0, max_distance=1.5):
    """
    Check if a given distance is between two cutoffs for a bond length.
    Input : 
        distance : bond length to check
        min_distance : minimum distance allowed - default 0. Angstroms 
        max_distance : maximum distance allowed - default 1.5 Angstroms
    Output :
        is_a_bond : Boolean (True or False) indicating if distance is between the cutoffs
    """
    # check the inputs
    if distance < 0.:
        raise ValueError("bond_check has detected a NEGATIVE distance! Check your inputs.")
    is_a_bond = False
    if (distance < max_distance) and (distance > min_distance):
        is_a_bond = True
    else:
        is_a_bond = False
    return is_a_bond

# calculate a distance for two atom positions
def calculate_distance(atom1_coord, atom2_coord):
    """
    function to calculate the distance between two atoms. 
    Inputs : 
        atom1_coord - 3D coordinates of atom 1
        atom2_coord - 3D coordinates of atom 2
    Outputs :
        distance between atom1 and atom2
    """
    # check the inputs
    if (len(atom1_coord) != 3) or (len(atom2_coord) != 3):
        raise ValueError("The shape of an atom's coordinates are incorrect. Double check your inputs")
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]    
    distance = numpy.sqrt(x_distance**2. + y_distance**2. + z_distance**2.)
    return distance
########################################################################

if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        raise NameError("Did NOT specify an input file. Add an argument for the file you want to analyze.")

    # get file name
    file_location = sys.argv[1]
    
    # extract symbols and coordinates
    symbols, coords = open_xyz(file_location)

    # process and print bonded atoms
    numcol = len(coords[0,:])
    for ii, coords_ii in enumerate(coords):
        for jj, coords_jj in enumerate(coords[ii:,:]):
            # calculate bond length
            bond_length12 = calculate_distance(coords_ii,coords_jj)
            if bond_check(bond_length12) is True:
                print(F'{symbols[ii]} to {symbols[jj]} : {bond_length12:.3f}')
