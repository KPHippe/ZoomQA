import os
import sys
import json
import numpy as np

from os.path import join

# THRESHOLD = 10.0 # 10 angstrum threshold, we can change this later
AA_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V','W', 'Y' ]

RADII = RADII = list(range(5, 56, 1))

def get_protein_contact_frequeny(casp_input, index):
    '''
    This method gets a 21x20 matrix representing the weighted frequency of contacts from the center amino acid to other amino acids based on the letter code/
    Row 0 represents the weighted contact frequency of the center in relation to all the amino acids (['A', 'C'...'Y']the columns) at radius 5. Row 1 is at radius 6 and so on.

    Parameters:
    --------------
    casp_input: dictionary
        This dictionary comes from one server prediction for one target (one JSON file) from Dr. Cao's CASP JSON database

    index: int
        This is the index of the target acid in the sequence


    Returns:
    ----------
    np.ndarray((21,20))
        This method gets a 21x20 matrix representing the weighted frequency of contacts from the center amino acid to other amino acids based on the letter code/
        Row 0 represents the weighted contact frequency of the center in relation to all the amino acids (['A', 'C'...'Y']the columns) at radius 5. Row 1 is at radius 6 and so on.

    '''
    contacts = _get_contacts(casp_input, index)

    contact_occurence_matrix = np.zeros((51,20)) #radius x amino acids
    contact_frequency_matrix = np.zeros((51,20)) #radius x amino acids
    for radius, contact_list in sorted(contacts.items()):
        fragment_indices = _get_fragment_indices(contact_list, index)
        # print(fragment_indices)
        for contact_pair in contact_list:
            contact_aa = contact_pair[1]
            contact_sequence_index = contact_pair[2]
            contact_occurence_index = AA_LIST.index(contact_aa)
            for fragment in fragment_indices:
                if contact_sequence_index in fragment:
                    if index in fragment:
                        contact_occurence_matrix[radius-5][contact_occurence_index] += 1
                    else:
                        contact_occurence_matrix[radius-5][contact_occurence_index] += 2

    for radius_row in range(len(contact_occurence_matrix)):
        total_row_contacts = np.sum(contact_occurence_matrix[radius_row])
        for col in range(len(contact_occurence_matrix[radius_row])):
            if total_row_contacts > 0:
                contact_frequency_matrix[radius_row][col] = contact_occurence_matrix[radius_row][col] / total_row_contacts

    return contact_frequency_matrix

def _get_fragment_indices(contact_list, center_index):
    '''
    This method is reponsible for determing the fragments present within a radius.

    We consider a fragment being any string of amino acids that are within the radius and are not seperated by any indices not present in the fragment.

    Parameters:
    ---------
    contact_list: list[(center_aa, contact_aa, contact_index)]
        This is a list of tuples, the length of however many residues are in contact with the targer aminoa acid at the center index. The tuples are as follows:
        the center aa letter code, the contact amino acid letter code, and the index in the sequence of the contact amino acid

    center_index: int
        This is the index of the target amino acid in the sequence

    Returns:
    --------
    list[list]:
        This method returns a list of lists. Each sub-list represents a fragment within the radius, all values represent residues in the sequence that are in contact
        with the center


    '''
    contact_indices = []
    for pair in contact_list:
        contact_indices.append(pair[2])

    contact_indices = sorted(contact_indices)
    fragment_indices = [[center_index]]
    for i in range(len(contact_indices)):
        comp_index = contact_indices[i]

        inside_frag = False
        for frag_ind in range(len(fragment_indices)):
            if comp_index -1 in fragment_indices[frag_ind]:
                fragment_indices[frag_ind].append(comp_index)
                inside_frag = True
                break
        #not in the current fragments
        if not inside_frag:
            fragment_indices.append([comp_index])
    return fragment_indices

def _get_contacts(casp_input, row):
    '''
    This method compiles a dictionary with every relevant data for the contacts for every radius for a target residue in the sequence

    Parameters:
    ------------
    aa_adjacency_list: list[length of sequence]
        This is a row of the contact list, index 0 represents the contact distance for the acid we are considering to the residue at index 0. Index 1 is the contact distance between
        the residue we are considering and residue at index 2

    row: int
        This represents the target amino acid index, it also represents the row of the contact map we are looking at

    Returns:
    ---------
    dictionary:
        Returns a dictionary with the keys being the radius and the values being a list of tuples. Each tuple is the following, the first index is the center amino acid letter code,
        the second index is the letter code for the amino acid in contact with the center, and the third value is the index in the sequence the contact resides.


    '''

    contacts= {}
    cm = casp_input['ContactMap']
    sequence = casp_input['aa']
    center_aa = sequence[row]
    for radius in RADII:
        contacts[radius] = []
        for col in range(len(cm[row])):
            contact_distance = cm[row][col]
            contact_aa = sequence[col]
            if contact_distance <= radius and row != col:
                contacts[radius].append((center_aa, contact_aa, col))

    return contacts

if __name__ == "__main__":
    # pathToCASP = '/media/kyle/IronWolf/CASP_ALL/'
    pathToCASP = '/Users/kylehippe/Documents/Summer2020/CASP14/Data/Neighborhoods/CASP12_json'

    target_paths = []

    for potential_path in os.listdir(pathToCASP):
        if 'tmp' in potential_path:
            continue
        target_paths.append(join(pathToCASP, potential_path))


    data = json.load(open(sorted(target_paths)[1]))


    for server_name, server_data in data.items():
        #test feature generation here
        for i in range(len(server_data['aa'])):

            ret = get_protein_contact_frequeny(server_data, 100)
            print(ret)
            break
        break
