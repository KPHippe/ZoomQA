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
        is initialized to a count of 0

    '''
    blank_dict = {}
    for radius in RADII:
        radius_dict = {}
        for aa in AA_LIST:
            radius_dict[aa] = 0

        radius_dict["total"] = 0
        blank_dict[radius] = radius_dict

    return blank_dict

def apply_proportions(dist_dict):
    '''
    This method takes in a dictionary that contains the counts of each amino acid and the total number of amino acids found
    and returns a new dictionary with the relative proportions of each amino acid instead of its counts.

    Parameters:
    ------------
    dist_dict:
        This dictionary has keys of each amino acid as well as 'total' and the values represent the occurrences of each amino acid. Total
        represents the total number of amino acids represented in this dictinoary.

    Returns:
    --------
    dictionary:
        This method a new dictionary with the relative proportions of each amino acid instead of its counts. This is a value between 0 and 1

    '''
    total = dist_dict['total']
    out = {}
    for aa, count in dist_dict.items():
        out[aa] = count / total

    out['total'] = total

    return out

def aa_change_from_json(target_data):
    '''
    This method takes one servers data and extracts the change over radius increase data for the amino acid densities.

    Parameters:
    -----------
    target_data: dictionary
        This dictionary is one server data from Dr. Cao's JSON database

    Returns:
    --------
    dictionary:
        This is a dictionary with keys mapping each index of the input sequence. The values is a nested dictionary with the keys being the radius and the values
        being the relative density of each amino acid within the radius specified by the key.


    '''
    cm = target_data['ContactMap']
    sequence = target_data['aa']
    pdb_radius_data = {}

    for row in range(len(cm)):
        cur_acid = sequence[row]
        local_radius_change = zero_local_radius_data()
#         print(local_radius_change.keys())
        for col in range(len(cm[row])):
            if row == col:
                continue
            comp_acid = sequence[col]

            radius_category = get_category(cm[row][col])
#             print(radius_category, col)
            if radius_category is not None:
                #TODO: the original end was 26
                for rad_to_add in range(radius_category, 56, 1):
                    local_radius_change[rad_to_add][comp_acid] += 1
                    local_radius_change[rad_to_add]['total'] += 1
        local_radius_dist = zero_local_radius_data()
        for radius, counts in local_radius_change.items():
            local_radius_dist[radius] = apply_proportions(counts)
        pdb_radius_data[row] = local_radius_dist
    return pdb_radius_data
