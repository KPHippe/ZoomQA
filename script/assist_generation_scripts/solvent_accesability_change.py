import os
import sys
import numpy as np

AA_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V','W', 'Y' ]

RADII = RADII = list(range(5, 56, 1))

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


def sol_change_from_json(target_data):
    '''
    This method takes one servers data and extracts the change over radius increase data for the average solvent accesability of the fragment structure
    as the radius increases.

    Parameters:
    -----------
    target_data: dictionary
        This dictionary is one server data from Dr. Cao's JSON database

    Returns:
    --------
    dictionary:
        This is a dictionary with the keys mapping to each index of the input sequence. The values are a dictionary with keys being the radius in range (5,25)
        and the values being the average solvent accesability of the structure with that radius


    '''
    cm = target_data['ContactMap']
    sequence = target_data['aa']
    pdb_radius_data = {}

    for row in range(len(cm)):
        cur_acid = sequence[row]
        local_radius_change = zero_local_radius_data()

        for col in range(len(cm[row])):
            #remove this comment hases below if you do not want to include the center amino acid solvent accesability
            # if row == col:
            #     continue
            comp_acid = sequence[col]

            radius_category = get_category(cm[row][col])

            if radius_category is not None:
                sol = target_data['sol'][col]
                norm_sol = np.clip(sol , 0, 300) / 300
                #TODO: the final radius was originally 26
                for rad_to_add in range(radius_category, 56, 1):
                    local_radius_change[rad_to_add].append(norm_sol)

        local_radius_dist = zero_local_radius_data()
        for radius, density_list in local_radius_change.items():
            local_radius_dist[radius] = np.mean(np.asarray(density_list))
        pdb_radius_data[row] = local_radius_dist
    return pdb_radius_data
