'''
Author: Kyle Hippe
        July 22, 2020
        NSSURP Summer Research

This program creates and saves the change-over-radius-increase data of CASP predictions in Dr. Cao's CASP JSON database.

Currently there are only 4 features, amino acid density of the structure per radius, average hydophobicity of the structure as the radius increases,
average mass of the structure as the radius increases, and average solvent accessibility of the structure as the radius increases.

The output of the program follows the following structure.
CASP_Fragment_Structures
    CASP6
        T0196
            ServerName_01.pkl
            ServerName_02.pkl
            …
            ServerName_n.pkl
        …(all other targets in CASP 6) T0224
            ServerName_01.pkl
            ServerName_02.pkl
            …
            ServerName_n.pkl
    CASP7
        T0322
            ServerName_01.pkl
            ServerName_02.pkl
            …
            ServerName_n.pkl
        …(all other targets in CASP 7)

        T0372
            ServerName_01.pkl
            ServerName_02.pkl
            …
            ServerName_n.pkl

    …


'''

import os
import sys
import copy
import json
import pickle
import traceback
import numpy as np
from os.path import join, getsize
from concurrent.futures import ProcessPoolExecutor, as_completed

from paths import PATHS

sys.path.insert(1, join(PATHS.sw_install, './script/assist_generation_scripts'))

from amino_acid_density_change import *
from hydrophobicity_change import *
from mass_change import *
from solvent_accesability_change import *
from isoelectricpoint_change import *
from vectorize_feature_data import *
from make_random_forest_predictions import *
from non_change_features import *
from contact_statistics import *
from structure_contact import *


def process_target(target_path, pathToSave):
    '''
    This method compiles all of the data from the scripts in assist_generation_scripts
    and compiles them into a dictionary with the following structure:

    Parameters:
    ----------
    target_path: string
        path to a json target file from Dr. Cao's JSON Database

    File: (dictionary)
        Keys -> index in sequence
        Values -> dictionary with relevant data
            Keys -> ['aa', 'hydro', 'mass', 'sol', 'localQA', 'sequenceAA'] update with new keys
            Values ->
                - 'aa_density_change' is a 21x 20 matrix. Each row in the matrix represents the radius we are considering
                  (eg. row 0 is a radius of 5 angstroms and row 20 is 25 angstrums radius) Each column represents the relative
                  density of the amino acids in alphabetical order within the radius structure.
                - 'hydro_change' is a vector with length 21. Each value in the matrix represents the average hydrophobicity
                  of the structure with a radius that corresponds to the index in the matrix (e.g. index 0 is the average
                  hydrophobicity of the structure with radius 5) This includes the aminoa acid hydrophobicity of the molecule
                  it is based around
                - 'mass_change' is a vector with length 21. Each value in the matrix represents the average mass of the structure
                   with a radius that corresponds to the index in the matrix (e.g. index 0 is the average mass of the structure with radius 5)
                   This includes the amino acid mass of the molecule it is based around
                - 'sol_change' is a vector with length 21. Each value in the matrix represents the average solvent accessibility
                   of the structure with a radius that corresponds to the index in the matrix (e.g. index 0 is the average solvent accessibility
                   of the structure with radius 5) This includes the amino acid sol of the molecule it is based around
                - 'iso_change' is a vector with length 21. Each value in the matrix represents the isoelectric point of the structure with a
                   radius that corresponds to the index in the matrix (e.g. index 0 is the average solvent accessibility of the structure with radius 5)
                - 'rf_predictions' is the random forest predictions for the stability of the amino acid based off of its torsion angles.
                   The first value corresponds to the stability regardless of secondary structure. The second value corresponds to the
                   stability if the secondary structure was a helix, the third if the secondary structure was a sheet, and the final if the secondary structure is a coil.
                - 'aa_mass' The normalized mass for the center amino acid
                - 'aa_hdyro' The normalized hydrophobicity for the center amino acid
                - 'aa_sol' The normalized solvent accessibility for the center amino acid
                - 'aa_iso' The normalized isoelectric point and pKa values for the center amino acid.
                   Index 0 is the pKa of the oxalic group. The second value is the pKa of the amino group.
                   The 3rd value is the pKa of the R-group if it exists, if not it is 0.0. The final value is the pI (isoelectric point) all values are normalized
                - 'psiphi' The normalized psi and phi angles for the center amino acid
                - 'aa_encoded' One-hot-encoded representation of the center amino acid
                - 'ss_encoded' One-hot-encoded representation of the secondary structure of the center amino acid
                - 'localQA' is a value between 0 and 1 representing the localQA of the residue at the center of the structure. 1 is perfect.
                - 'average_distance' a length 21 vector representing the average distance between the center and all acids in contact with the center.
                   This is normalized by the radius of gyration (the longest distance still considered by the contact)
                - 'std_dev_distance' a length 21 vector representing the std_deviation of the distance between the center and all acids in contact
                   with the center. The distances used for std dev calculation are normalized by the radius of gyration (the longest distance still 
                   considered by the contact), so the std dev will be between 0 and 1
                - 'percent_contact' a length 21 vector representing the total percentage of the protein that is in contact with the center as the
                   radius increases (e.g. index 0 is the percentage of the protein in contact with the center amino acid with radius 5)
                - 'structure_contact_matrix' a 21x20 matrix. Row 0 represents the weighted contact frequency of the center in relation to all the
                   amino acids (the columns) at radius 5. Row 1 is at radius 6 and so on.

    This is then saved to the pathToSave location
    '''
    json_data = load_json_file(target_path)
    casp_name = target_path.split("/")[-1].split("_")[0]
    target_name = target_path.split("_")[-1]
    create_file(join(pathToSave, casp_name))
    create_file(join(pathToSave, casp_name, target_name))
    for server, server_data in json_data.items():
        try:

            server_name = server.split(":")[-1]
            server_save = join(pathToSave, casp_name, target_name, f"{server_name}.pkl")

            sequence = server_data['aa']
            aa_data = aa_change_from_json(server_data)
            hydro_data = hydro_change_from_json(server_data)
            mass_data = mass_change_from_json(server_data)
            sol_data = sol_change_from_json(server_data)
            iso_data = iso_change_from_json(server_data)

            server_vectors = vectorize_pdb_data(aa_data, hydro_data, mass_data, sol_data, iso_data, sequence)

            for index in server_vectors.keys():
                # add a few comments here to describe what it adds

                local_qa = server_data['localQA'][index]
                sequence_aa = server_data['aa'][index]
                local_ss = server_data['ss'][index]

                local_psi, local_phi = int(float(server_data['Angles']['psi_im1'][index])), int(
                    float(server_data['Angles']['phi'][index]))
                rf_predictions = get_prediction((local_psi + 180), (local_phi + 180),
                                                sequence_aa)  # have to add 180 because its in a ramachandran plot
                server_vectors[index]['rf_predictions'] = rf_predictions

                non_change_data = get_non_change_features(server_data, index)
                for key, data_values in non_change_data.items():
                    server_vectors[index][key] = data_values

                radius_change_statistics = get_contact_stats(server_data, index)
                for key, data_values in radius_change_statistics.items():
                    server_vectors[index][key] = data_values

                structure_contact_matrix = get_protein_contact_frequeny(server_data, index)
                server_vectors[index]['structure_contact_matrix'] = structure_contact_matrix

            pickle.dump(server_vectors, open(server_save, 'wb'))
            print(f"Saved {server_name} to {server_save}")
        except Exception as e:
            print(f"Error creating {target_name}")


def load_json_file(target_path):
    return json.load(open(target_path))


def main(pathToData, pathToRandomForestPredictions, pathToSave):
    # load the random forest models so we don't have to distribute a list of them
    load_RF_predictions(pathToRandomForestPredictions)

    targets_list = os.listdir(pathToData)
    path_list = []
    for target in targets_list:
        if "tmp" in target:
            # this is a tmprun file, does not actually contain data
            continue
        path_list.append(join(pathToData, target))

    create_file(pathToSave)
    print('Saving data...')
    with ProcessPoolExecutor(max_workers=int(os.cpu_count() * 0.70)) as executor:
        executor.map(process_target, path_list, [pathToSave] * len(path_list))

    # for target_path in path_list:
    #     process_target(target_path, pathToSave)
    #     break

    print("\n\nData generation complete")


def create_file(pathToFile):
    try:
        os.mkdir(pathToFile)
    except FileExistsError:
        print(f"{pathToFile} already exists...")
    except:
        print(f"Fatal error making {pathToFile}...")
        sys.exit()


if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("Not enough arguemnts, example command: ")
        print(
            f"python {sys.argv[0]} /data/shared/databases/CASP_ALL_JSON /data/summer2020/Kyle/CASP14/Data/Angles/AminoAcid_RF/RF_Predictions /data/summer2020/Kyle/CASP14/Data/Graphs/CASP_Fragment_Databse/ ")

        sys.exit()

    pathToData = sys.argv[1]
    pathToRandomForestPredictions = sys.argv[2]
    pathToSave = sys.argv[3]

    main(pathToData, pathToRandomForestPredictions, pathToSave)
