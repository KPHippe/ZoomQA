import os
import sys
import copy 
import json
import pickle
import traceback
import numpy as np
from os.path import join, isfile, isdir, getsize


from sklearn.ensemble import RandomForestRegressor

allstruct_predictions = {}
helix_predictions = {}
sheet_predictions = {}
coil_predictions = {}


def get_prediction(target_psi, target_phi, aa):
    # TODO: update this documentation 
    '''
    This model makes a prediction on the input psi, phi data 

    Parameters: 
    -----------
    input_data: tuple(int, int)
        The input data is a tuple of integers. The first index represents the psi angle, the second index is the phi angle 

    aa: string -> length 1 
        This is the letter code representation of the amino acid we are testing 

    Return: 
    ---------
    tuple: (allstructure_stability, helix_stability, sheet_stability, coil_stability)
        This tuple is the stability predictions for the allstructure models, helix models, sheet models, and coil models that correspond to the input amino acid

    
    '''
    target_all_struct = allstruct_predictions[aa][target_psi][target_phi]
    target_helix = helix_predictions[aa][target_psi][target_phi]
    target_sheet = sheet_predictions[aa][target_psi][target_phi]
    target_coil = coil_predictions[aa][target_psi][target_phi]

    return (target_all_struct, target_helix, target_sheet, target_coil)


def load_RF_predictions(pathToPredictions):
    # TODO: redo the documentation on this
    '''
    This method established the models in the global variable 'models' It is done this way 
    so that when we multiprocess, we do not have to have 1000+ instances of the models, the exist
    in the persistent version of this file that gets distributed through the multiprocessing. 

    Parameters: 
    ----------
    pathToModels: string
        This is the string representation to the models folder 

    Returns: 
    ---------
    None
        This function does not return models, but establishes the global models variable with a structure as follows: 

        models: dictionary
        keys -> structure ('allstruct', 'helix', 'sheet', 'coil')
        values -> dictionary
            key -> amino acid ('A', 'C', ... 'Y')
            value -> random forest model corresponding to the structure and amino acid of the keys    
    
    '''
    global allstruct_predictions
    global helix_predictions
    global sheet_predictions
    global coil_predictions
    for pred_file in os.listdir(pathToPredictions):
        if 'allstruct' in pred_file:
            allstruct_predictions = pickle.load(open(join(pathToPredictions, pred_file), 'rb'))
        if 'helix' in pred_file:
            helix_predictions = pickle.load(open(join(pathToPredictions, pred_file), 'rb'))
        if 'sheet' in pred_file:
            sheet_predictions = pickle.load(open(join(pathToPredictions, pred_file), 'rb'))
        if 'coil' in pred_file:
            coil_predictions = pickle.load(open(join(pathToPredictions, pred_file), 'rb'))



if __name__ == "__main__":
    #change the paths
    load_RF_predictions('/media/kyle/Samsung860Evo/Summer2020/CASP14/Data/Angles/AminoAcid_RF/RF_Predictions')
    pathToData = '/media/kyle/IronWolf/CASP_ALL/'

    target = json.load(open(join(pathToData, sorted(os.listdir(pathToData))[10])))

    print(allstruct_predictions.keys())
    print(helix_predictions.keys())
    print(sheet_predictions.keys())
    print(coil_predictions.keys())


    for server, data in target.items():
        if '.pdb' not in server: 
            continue
        print(server)
        index = 0
        psi, phi = int(float(data['Angles']['psi_im1'][index])), int(float(data['Angles']['phi'][index]))
        print(psi, phi)
        print(data['ss'][index])
        print(make_prediction((psi+180), (phi + 180), data['aa'][index]))



    
    