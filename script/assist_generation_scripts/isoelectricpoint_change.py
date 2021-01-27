'''
Isoelectric point data taken from
https://www.anaspec.com/html/pK_n_pl_Values_of_AminoAcids.html
'''

import os
import sys
import json
import numpy as np
from os.path import join

from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA

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


def iso_change_from_json(target_data):
    '''
    This method takes one servers data and extracts the change over radius increase data for the isoelectric point of a fragment as the
    radius increases.

    Parameters:
    -----------
    target_data: dictionary
        This dictionary is one server data from Dr. Cao's JSON database

    Returns:
    --------
    dictionary:
        This is a dictionary with the keys mapping to each index of the input sequence. The values are a dictionary with keys being the radius in range (5,25)
        and the values being the isoelectric point of the structure with that radius


    '''
    cm = target_data['ContactMap']
    sequence = target_data['aa']
    pdb_radius_data = {}

    for row in range(len(cm)):
        cur_acid = sequence[row]
        local_radius_change = zero_local_radius_data()
#         print(local_radius_change.keys())
        for col in range(len(cm[row])):
            #remove this comment hases below if you do not want to include the center amino acid hydrophobicity
            # if row == col:
            #     continue
            comp_acid = sequence[col]

            radius_category = get_category(cm[row][col])
#             print(radius_category, col)
            if radius_category is not None:
                #TODO: the original was 26
                for rad_to_add in range(radius_category, 56, 1):
                    local_radius_change[rad_to_add].append(comp_acid)

        local_radius_iep = zero_local_radius_data()
        for radius, aa_list in local_radius_change.items():
            radius_seq = "".join([acid for acid in aa_list])
            radius_aa_content = PA(radius_seq).count_amino_acids()
            temp_protein = IP(radius_seq, radius_aa_content)
#             print(temp_protein.pi())
            temp_protein_pi = temp_protein.pi()
            norm_protein_pi = np.clip([(temp_protein_pi - 2.98) / (10.76 - 2.98)], 0.0, 1.0)[0]
            local_radius_iep[radius] = norm_protein_pi

        pdb_radius_data[row] = local_radius_iep
    return pdb_radius_data




if __name__ == "__main__":
    pathToData = '/media/kyle/IronWolf/CASP_ALL/'

    target = json.load(open(join(pathToData, sorted(os.listdir(pathToData))[10])))


    for server, data in target.items():
        server_data = iso_change_from_json(data)
        # print(server_data.keys())
        # print(server_data[100])

        vectorized = _vectorize_local_iso(server_data, ''.join([i for i in data['aa']]), 0)

        print(len(vectorized))
        print(vectorized)

        vectorized = _vectorize_local_iso(server_data, ''.join([i for i in data['aa']]), 110)

        print(len(vectorized))
        print(vectorized)
        break
