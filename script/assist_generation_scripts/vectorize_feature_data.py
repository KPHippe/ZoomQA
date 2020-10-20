import os 
import sys
import numpy as np


def vectorize_pdb_data(aa_data, hydro_data, mass_data, sol_data, iso_data, sequence):
    '''
    This method is responsible for vectorizing the pdb data extractedfrom the respective 'attribute'_change.py files
    and returning a dictionary containing all of the vector data

    Parameters:
    ----------
    aa_data: dictionary
        aa_data is a dictionary with keys mapping each index of the input sequence. The values is a nested dictionary with the keys being the radius and the values 
        being the relative density of each amino acid within the radius specified by the key. 

    hydro_data: dictionary
        hydro_data is a dictionary with the keys mapping to each index of the input sequence. The values are a dictionary with keys being the radius in range (5,25)
        and the values being the average hydrophobicity of the structure with that radius 

    mass_data: dictionary
        mass_data is a dictionary with the keys mapping to each index of the input sequence. The values are a dictionary with keys being the radius in range (5,25)
        and the values being the average mass of the structure with that radius 

    sol_data: dictionary
        sol_data is a dictionary with the keys mapping to each index of the input sequence. The values are a dictionary with keys being the radius in range (5,25)
        and the values being the average solvent accessibility of the structure with that radius

    iso_data: dictionary 
        iso_data is a dictionary with the keys mapping to each index of the input sequence. The values are a dictionary with keys being the radius in range (5,25)
        and the values being the isoelectric point of the structure with that radius

    sequence: list[char]
        sequence is a list of chars representing the amino acid sequence of the input PDB


    Returns: 
    Dictionary: 
        The return dictionary has the following keys: ['aa', 'hydro', 'mass', 'sol', 'localQA', 'sequenceAA']. 'aa' is a matrix
        representing the relative density of each amino acid at each radius level, 'hydro', 'mass', 'sol' are all length 21 vectors representing the 
        average 'feature' value of the structure at each radius is range (5,25)
    '''
    
    if len(aa_data) != len(hydro_data) or len(hydro_data) != len(mass_data) or len(mass_data) != len(sol_data):
        print("Incorrect length of input data, cannot vectorize")
        return None

    out_dictionary = {}
    for index in aa_data.keys():
        local_dict = {}

        target_acid = sequence[index]

        local_dict["aa_density_change"] = _vectorize_local_aa(aa_data, sequence, index)
        local_dict['hydro_change'] = _vectorize_local_hydro(hydro_data, sequence, index)
        local_dict['mass_change'] = _vectorize_local_mass(mass_data, sequence, index)
        local_dict['sol_change'] = _vectorize_local_sol(sol_data, sequence, index)
        local_dict['iso_change'] = _vectorize_local_iso(iso_data, sequence, index)

        out_dictionary[index] = local_dict

    return out_dictionary

def _vectorize_local_aa(acid_input, sequence, index):
    '''
    This method is repsonsible for getting the matrix that represents the relative density of each amino acid (columns) as the radius 
    increases (rows) 

    Parameters: 
    -----------
    acid_input: dictionary
        This is the relative density of each amino acid for a fragment generated by 'amino_acid_density_change.py' script 
        key-> radius
        value -> dictionary
            key -> amino acid letter code
            value -> relative density of that aminoa acid in this fragment 

    sequence: list[char]
        This is a list of chars representing the sequence of the sequence this PDB codes for

    index: int
        This is the index we are extracting data for in the sequence 

    Returns: 
    -------
    np.ndarray((21,21)): float64
        The return is a 21x21 matrix. The colums represent the relative density (between 0 and 1) of the amino acids (in alphabetical order) as the radius (rows)
        increases


    * Notes
        - we keep the total in the radius here, but that will most likely have to be changed, Dr. Cao and I have had some ideas, in my notes in the google docs
    '''
    
    acid_target = acid_input[index]
    target_acid = sequence[index]
    
    #vectorize the casp input
    aa_vector = []
    for radius, dist_dict in sorted(acid_target.items()):
        if type(radius) is str: 
            #total not relevant here
            continue
        local_radius_vector = []
        for acid, local_dist in sorted(dist_dict.items()):
            if 'total' in acid:
                # we describe this feature elsewhere
                continue
            local_radius_vector.append(local_dist)
        aa_vector.append(np.asarray(local_radius_vector))
    aa_vector = np.asarray(aa_vector)

    return aa_vector


def _vectorize_local_hydro(hydro_input,  sequence, index):
    '''
    This method is repsonsible for getting the vector that represents the average hydrophobicity of the structure specified 
    by the index as the radius increases from 5 to 25 angstroms

    Parameters: 
    -----------
    hydro_input: dictionary
        This is the average hydrophobicity change over radius increase data created for each server by the 'hydrophobicity_change.py' script 
        key-> radius
        value -> float, average hydrophobicity 

    sequence: list[char]
        This is a list of chars representing the sequence of the sequence this PDB codes for

    index: int
        This is the index we are extracting data for in the sequence 

    Returns: 
    -------
    np.ndarray((21,)): float64
        The return is a length 21 vector. Each element in the vectore represents the fragment with a radius (index + 5) angstroms average hydrophobicity


    '''
    
    hydro_target = hydro_input[index]
    target_acid = sequence[index]
    
    #vectorize the casp input
    hydro_vector = []
    for radius, ave_hydro in sorted(hydro_target.items()):
        if type(radius) is str: 
            #total, not relevant here
            continue
        hydro_vector.append(ave_hydro)
    hydro_vector = np.asarray(hydro_vector)
#     print(casp_vector)
#     print(casp_vector.shape)

    
    return hydro_vector


def _vectorize_local_mass(mass_input, sequence,  index):
    '''
    This method is repsonsible for getting the vector that represents the average mass of the structure specified 
    by the index as the radius increases from 5 to 25 angstroms 

    Parameters: 
    -----------
    mass_input: dictionary
        This is the average mass change over radius increase data created for each server by the 'mass_change.py' script 
        key-> radius
        value -> float, average mass for that radius
    sequence: list[char]
        This is a list of chars representing the sequence of the sequence this PDB codes for

    index: int
        This is the index we are extracting data for in the sequence 

    Returns: 
    -------
    np.ndarray((21,)): float64
        The return is a length 21 vector. Each element in the vectore represents the fragment with a radius (index + 5) angstroms average mass
         

    '''
    mass_target = mass_input[index]
    target_acid = sequence[index]
    
    #vectorize the casp input
    mass_vector = []
    for radius, ave_mass in sorted(mass_target.items()):
        if type(radius) is str: 
            #total, not relevant here
            continue
        mass_vector.append(ave_mass)
    mass_vector = np.asarray(mass_vector)
#     print(casp_vector)
#     print(casp_vector.shape)

    return mass_vector

def _vectorize_local_sol(sol_input, sequence, index):
    '''
    This method is repsonsible for getting the vector that represents the average solvent accessibility of the structure specified 
    by the index as the radius increases from 5 to 25 angstroms 

    Parameters: 
    -----------
    sol_input: dictionary
        This is the solvent change over radius increae data created for each server by the 'solvent_accessability_change.py' script 
        key-> radius
        value -> float, average solvent accesability for that radius

    sequence: list[char]
        This is a list of chars representing the sequence of the sequence this PDB codes for

    index: int
        This is the index we are extracting data for in the sequence 

    Returns: 
    -------
    np.ndarray((21,)): float64
        The return is a length 21 vector. Each element in the vectore represents the fragment with a radius (index + 5) angstroms average solvent accessibility
         

    '''
    
    sol_target = sol_input[index]
    traget_acid = sequence[index]
    
    #vectorize the casp input
    sol_vector = []
    for radius, sol in sorted(sol_target.items()):
        if type(radius) is str: 
            #total, not relevant here
            continue
        sol_vector.append(sol)
    sol_vector = np.asarray(sol_vector)
#     print(casp_vector)
#     print(casp_vector.shape)
    
    return sol_vector


def _vectorize_local_iso(iso_input, sequence, index):
    '''
    This method is repsonsible for getting the vector that represents the isoelectric point of the structure specified 
    by the index as the radius increases from 5 to 25 angstroms 

    Parameters: 
    -----------
    iso_input: dictionary
        This is the solvent change over radius increae data created for each server by the 'solvent_accessability_change.py' script 
        key-> radius
        value -> float, isoelectric point for that radius

    sequence: list[char]
        This is a list of chars representing the sequence of the sequence this PDB codes for

    index: int
        This is the index we are extracting data for in the sequence 

    Returns: 
    -------
    np.ndarray((21,)): float64
        The return is a length 21 vector. Each element in the vectore represents the fragment with a radius (index + 5) angstroms isoelectric point
         

    '''
    iso_target = iso_input[index]
    target_acid = sequence[index]
    
    #vectorize the casp input
    iso_vector = []
    for radius, local_iso in sorted(iso_target.items()):
        if type(radius) is str: 
            #not relevant here
            continue
        iso_vector.append(local_iso)
    iso_vector = np.asarray(iso_vector)

    return iso_vector