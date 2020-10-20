'''
This file is reponsible for compiling the data about the center amino acid, this information is not any change analysis, just physical/chemical properties
'''

import os 
import sys
import numpy as np

import json
from os.path import join, getsize

aa_one_hot_encode = {
    "A" : [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    "C" : [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    "D" : [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    "E" : [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    "F" : [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    "G" : [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    "H" : [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
    "I" : [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
    "K" : [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
    "L" : [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
    "M" : [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
    "N" : [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
    "P" : [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
    "Q" : [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
    "R" : [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
    "S" : [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
    "T" : [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
    "V" : [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
    "W" : [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
    "Y" : [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
}
ss_one_hot_encode = {
    "H" : [1,0,0],
    "E" : [0,1,0],
    "C" : [0,0,1]
}
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
hydrophobicity = {
    "A" : 47,
    "C" : 52,
    "D" : -18,
    "E" : 8,
    'F' : 92,
    'G' : 0,
    'H' : 8,
    'I' : 100,
    'K' : -37,
    'L' : 100,
    'M' : 74,
    'N' : -41,
    'P' : -46,
    'Q' : -18,
    'R' :  47,
    'S' : -7,
    'T' : 13,
    'V' : 79,
    'W' : 84,
    'Y' : 49
}

isoelectric_point = {
    #format: pK alpha-CO2H, pK NH3, pK R-group (not always present, if not -> 0), pI (overall isoelectric point at 25C)
    #taken from: https://www.anaspec.com/html/pK_n_pl_Values_of_AminoAcids.html
    "A" : [2.35, 9.87, 0.0, 6.11],
    "C" : [1.71, 10.78, 8.33, 5.02],
    "D" : [1.88, 9.60, 3.65, 2.98],
    "E" : [2.19, 9.67, 4.25, 3.08],
    'F' : [2.58, 9.24, 0.0, 5.91],
    'G' : [2.34, 9.60, 0.0, 6.06],
    'H' : [1.78, 8.97, 5.97, 7.64],
    'I' : [2.32, 9.76, 0.0, 6.04],
    'K' : [2.20, 8.90, 10.28, 9.47],
    'L' : [2.36, 9.60, 0.0, 6.04],
    'M' : [2.28, 9.21, 0.0, 5.74],
    'N' : [2.18, 9.09, 13.2, 10.76],
    'P' : [1.99, 10.60, 0.0, 6.30],
    'Q' : [2.17, 9.13, 0.0, 5.65],
    'R' : [2.18, 9.09, 13.2, 10.76],
    'S' : [2.21, 9.15, 0.0, 5.68],
    'T' : [2.15, 9.12, 0.0, 5.60],
    'V' : [2.29, 9.74, 0.0, 6.02],
    'W' : [2.38, 9.39, 0.0, 5.88],
    'Y' : [2.20, 9.11, 10.07, 5.63]
}

def get_non_change_features(casp_server_input, index):
    '''
    This method is responsible for aquiring the center amino acid physical and chemical information

    Parameters: 
    -----------
    casp_server_input: dictionary
        This dictionary is a single server from a single target (one of Dr. Cao's JSON files in the database)
    index: int
        This is the index of the center amino acid we are looking at. 

    Return: 
    ----------
    dictionary: 
        This returns a dictionary with keys being feature names, and values being thier corresponding feature 
        key -> feature name
        value -> feature value

    '''
    local_amine_return = {}

    local_qa = _normalize_lqa(casp_server_input['localQA'][index])
    sequence_aa = casp_server_input['aa'][index]
    local_ss = casp_server_input['ss'][index]
    local_psi, local_phi = int(float(casp_server_input['Angles']['psi_im1'][index])), int(float(casp_server_input['Angles']['phi'][index]))
    local_sol = casp_server_input['sol'][index]

    #normalize the mass
    norm_mass = (monoisotopic_mass[sequence_aa] - 57.02146)/(186.0793 - 57.02146)
    norm_hydro = (hydrophobicity[sequence_aa] + 46) / 146 
    norm_sol = np.clip(local_sol , 0, 300) / 300
    norm_iso = _normalize_iso(isoelectric_point[sequence_aa])
    norm_psi, norm_phi = ((local_psi + 180)/360), ((local_phi + 180)/360)

    local_amine_return['aa_mass'] = norm_mass
    local_amine_return['aa_hydro'] = norm_hydro
    local_amine_return['aa_sol'] = norm_sol
    local_amine_return['aa_iso'] = norm_iso
    local_amine_return['psiphi'] = (norm_psi, norm_phi)
    local_amine_return['aa_encoded'] = aa_one_hot_encode[sequence_aa]
    local_amine_return['ss_encoded'] = ss_one_hot_encode[local_ss.upper()]
    local_amine_return['local_qa'] = local_qa

    return local_amine_return
 

def _normalize_lqa(qa_score):
    '''
    Function to represent localqa as a score between 0 and 1

    Parameters: 
    ----------
    qa_score: float
        This is the float value representing the qa score of a specific amino acid

    Returns: 
    float: (0 <= x <= 1)
        The return is a normalized value between 0 and 1 corresponding to the input qa score, but normalized

    '''
    return 1/(1 + (qa_score * qa_score/12) )


def _normalize_iso(iso_vector):
    '''
    Function to normalize the pKa values and pI values for an amino acid

    Parameters: 
    ------------
    iso_vector: list([oxalic pKa, amine pKa, r-group pKa, pI])
        This is a list containing the relevant oxalix pKa, amine pKa, r-group pKa, and isoelectric point of the target amino acid

    Return: 
    -------
    list: [norm_oxalic pKa, norm amine pKa, norm r-group pKa, norm pI]
        Returns the normalized values in the same order as the input

    '''
    oxalic_norm = (iso_vector[0] - 1.71) / (2.58 - 1.71)
    amine_norm = (iso_vector[1] - 8.90) / (10.78 - 8.90)
    r_norm = (iso_vector[2] - 0.0)  / (13.20 - 0.0)
    pI_norm = (iso_vector[3] - 2.98) / (10.76 - 2.98)

    return (oxalic_norm, amine_norm, r_norm, pI_norm)


if __name__ == "__main__":
    pathToCASP = '/media/kyle/IronWolf/CASP_ALL/'

    target_paths = []

    for potential_path in os.listdir(pathToCASP):
        if 'tmp' in potential_path:
            continue
        target_paths.append(join(pathToCASP, potential_path))


    data = json.load(open(sorted(target_paths)[0]))


    for server_name, server_data in data.items():
        #test feature generation here 
        ret = get_non_change_features(server_data, 0)
        print(ret)
        break