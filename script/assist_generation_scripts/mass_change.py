import os
import sys
import numpy as np

AA_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V','W', 'Y' ]

RADII = RADII = list(range(5, 56, 1))

monoisotopic_mass = {
    "A" : 71.03711,
    "C" : 103.00919,
    "D" : 115.02694,
    "E" : 129.04259,
    'F' : 147.06841,
    'G' : 57.02146,
    'H' : 137.05891,
    'I' : 113.08406,
    'K' : 128.09496,
    'L' : 113.08406,
    'M' : 131.04049,
    'N' : 114.04293,
    'P' : 97.05276,
    'Q' : 128.05858,
    'R' : 156.10111,
    'S' : 87.03203,
    'T' : 101.04768,
    'V' : 99.06841,
    'W' : 186.0793,
    'Y' : 163.06333
}

def get_category(distance):
    '''
    This method is responsible for finding the distance category a distance falls into

    Parameters:
    ----------
    distance: float
        This represents the distance in angstroms a target amino acid is away from the current center

    Returns:
    --------
    int/None
        This integer represents the radius category that the distance falls into (e.g. a radius of 6.7 would return 7 meaning
        that the target amino acid is within a radius of 7 angstroms). It will return None if the distance is greater
        than the threshhold set by RADII

    '''
    for radius in RADII:
        if distance <= radius:
            return radius
    return None

def zero_local_radius_data():
    '''
    This method initializes a blank dictionary for the change function to fill

    Returns:
    ----------
    dictionary:
        This dictionary is nested. The outer dictionary has keys equal to the radius and the innder dictionary has keys of each amino acid and
        is initialized to a value of 0

    '''
    blank_dict = {}

    for radius in RADII:
        blank_dict[radius] = []

    return blank_dict

def mass_change_from_json(target_data):
    '''
    This method takes one servers data and extracts the change over radius increase data for the average mass of the fragment structure
    as the radius increases.

    Parameters:
    -----------
    target_data: dictionary
        This dictionary is one server data from Dr. Cao's JSON database

    Returns:
    --------
    dictionary:
        This is a dictionary with the keys mapping to each index of the input sequence. The values are a dictionary with keys being the radius in range (5,25)
        and the values being the average mass of the structure with that radius


    '''
    cm = target_data['ContactMap']
    sequence = target_data['aa']
    pdb_radius_data = {}

    for row in range(len(cm)):
        cur_acid = sequence[row]
        local_radius_change = zero_local_radius_data()
        for col in range(len(cm[row])):
            #remove this comment hases below if you do not want to include the center amino acid mass
            # if row == col:
            #     continue
            comp_acid = sequence[col]

            radius_category = get_category(cm[row][col])
            if radius_category is not None:
                norm_mass = norm_mass = (monoisotopic_mass[comp_acid] - 57.02146)/(186.0793 - 57.02146)
                #TODO: the original stop was 26
                for rad_to_add in range(radius_category, 56, 1):
                    local_radius_change[rad_to_add].append(norm_mass)

        local_radius_dist = zero_local_radius_data()
        for radius, density_list in local_radius_change.items():
            local_radius_dist[radius] = np.mean(np.asarray(density_list))
        pdb_radius_data[row] = local_radius_dist
    return pdb_radius_data
