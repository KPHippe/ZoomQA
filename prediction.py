import os
import re
import sys
import json
import math
import pickle
import shutil
import subprocess
from os.path import join, isdir, isfile
from timeit import default_timer as timer

import numpy as np
from sklearn.svm import SVR

from script.generate_formatted_SVR_input import parse_server_data

PYTHON_INSTALL = 'python'
SW_INSTALL = './'

TOP_N = 100
ZOOMQA = '''\


  _____                      ___      _    
 |__  /___   ___  _ __ ___  / _ \    / \   
   / // _ \ / _ \| '_ ` _ \| | | |  / _ \  
  / /| (_) | (_) | | | | | | |_| | / ___ \ 
 /____\___/ \___/|_| |_| |_|\__\_\/_/   \_\
                                           


'''
def preprocess_input(pathToInput, pathToSave):
    '''
    This method is responsible for taking the pdb input files and extract
    all of the necesary features into the pickle files that can be easily
    transformed into the correct format for the model input, saves it to a temp
    directory in the output folder

    Parameters:
    ----------------
    pathToInput: string
        This is a string representation to the path to the input data

    pathToSave: string
        This is a string representation to the path to the save folder, so we can
        make a temp folder to store the intermediary steps

    Return:
    ---------------
    type: string
        The return is the path to the final step (step2_generate_casp_fragment_structures)
        output so that the next step can use it as the input

    '''
    pattern = re.compile(r"T\d{4}[a-zA-Z]*[0-9]*")
    target_name = re.search(pattern, pathToInput)

    if target_name is not None: 
        target_name = str(target_name[0])
    else: 
        target_name = 'Target'

    pathToTempDirectory = join(pathToSave, 'tmp')
    create_folder(pathToTempDirectory)


    pathToStep0 = join(pathToTempDirectory, 'step_0')
    pathToJSON = join(pathToTempDirectory, 'JSON_Data')
    pathToZoomQAInputData = join(pathToTempDirectory, 'ZoomQA_Input')

    print("Processing input data...")
    # chain_add_command = f'python ./script/step0_prepare_add_chain_to_folder.py ./script/assist_add_chainID_to_one_pdb.pl {pathToInput} {pathToStep0} > {join(pathToTempDirectory, "step0_log.txt")} 2>&1'
    create_folder(pathToStep0)
    pathToStep0OUT = join(pathToStep0, target_name)
    #this was originally step0_prepare_add_chain_to_folder.py 
    step0_location = join(SW_INSTALL, 'script/step0_prepare_add_chain_to_folder_v2.py')
    chain_add_location = join(SW_INSTALL, 'script/assist_add_chainID_to_one_pdb.pl')
    chain_add_command = f'{PYTHON_INSTALL} {step0_location} {chain_add_location} {pathToInput} {pathToStep0OUT} >/dev/null 2>&1'
    os.system(chain_add_command)

    #change it to _linux for linux run, mac for mac run
    # json_command = f'python ./script/step1_create_json_from_PDB.py ./script/stride_mac {pathToStep0} {pathToJSON} > {join(pathToTempDirectory, "step1_log.txt")} 2>&1'
    step1_location = join(SW_INSTALL, 'script/step1_create_json_from_PDB.py')
    stride_location = join(SW_INSTALL, 'script/stride_linux')
    json_command = f'{PYTHON_INSTALL} {step1_location} {stride_location} {pathToStep0} {pathToJSON} >/dev/null 2>&1'
    os.system(json_command)

    step2_location = join(SW_INSTALL, 'script/step2_generate_casp_fragment_structures.py')
    rfpredictions_locations = join(SW_INSTALL, 'script/assist_generation_scripts/RF_Predictions/')
    frag_structure_command = f'{PYTHON_INSTALL} {step2_location} {pathToJSON} {rfpredictions_locations} {pathToZoomQAInputData} >/dev/null 2>&1'
    os.system(frag_structure_command)

    return f'{pathToZoomQAInputData}'


def load_input_data(pathToData):
    '''
    This method is repsonsible for loading the output from preprocess_input and getting
    it staged for going into the model

    Parameters:
    ---------------
    pathToData: string
        This is a path to the tmp folder returned from preprocess_input(), it holds
        all of the input data

    Returns:
    -----------
    dictionary {string, dictionary{feature_name: list}}
        The return is a dictionary where the keys are the server names, and the values
        are dictionaries of input features, will get formatted just before predictions are made

    '''
    #check if this path is to the data, or just to the folder containing it
    '''I need a better way to do this, this is messy and unreliable, but works for now'''
    pathToServers = pathToData
    folder_contents = os.listdir(pathToServers)
    while isdir(join(pathToServers, folder_contents[0])):
        pathToServers = join(pathToServers, folder_contents[0])
        folder_contents=os.listdir(pathToServers)

    server_names = os.listdir(pathToServers)
    # print(pathToServers)

    input_data = {}
    for server_name in server_names:
        pathToServer = join(pathToServers, server_name)
        data = pickle.load(open(pathToServer, 'rb'))
        input_data[server_name] = data

    return input_data

def load_model(pathToModel):
    '''
    This method loads a pre-trained model for QA predictions

    Parameters:
    ---------------
    pathToModel: string
        This is a path to the pre-trained model

    Return:
    ---------
    SVR model:
        This is the SVR trained on the top 100 features

    '''

    with open(pathToModel, 'rb') as f:
        model = pickle.load(f)

    print("Model parameters: ")
    print(model)

    return model

def make_predictions(model, input_data):
    '''
    This method is responsible for making predictions on the input data

    Parameters:
    --------------
    model: SVR model
        The pretrained model for QA prediction

    input_data: dictionary {string:dictionary:{feature_name:list}}
        this is a dictionary with the server name as the key, and the value is
        the input feature dictionary, with keys of feature names and values are lists with feature Values

    Return:
    dictionary: {string: list[float]}
        The return is a dictionary with keys being the name of the input server, and the
        values are a list that correspond to the distance in angstroms of the predicted
        amino acid compared to an unknown ground truth

    '''

    def _un_norm_qa(score):
        np.clip(score, 0, 1)
        # sqrt(((1/norm)-1) * 12)
        try:
            return math.sqrt(((1/np.clip(score, 1e-5, 1))-1)*12)
        except:
            print(score)
    predictions = {}
    for server_name, whole_target_data in input_data.items():
        server_prediction = []
        #turn data into correct input form
        server_X, server_y = parse_server_data(whole_target_data, TOP_N)

        #get predictions
        server_prediction_normalized = model.predict(server_X)
        #convert scores to distance
        server_prediction_distance = []
        for prediction in server_prediction_normalized:
            server_prediction_distance.append(_un_norm_qa(prediction))

        predictions[server_name] = server_prediction_distance

    return predictions


def write_predictions(prediction_data, pathToSave, target_name):
    '''
    This method writes out the predictions in CASP format

    Parameters:
    -------------
    prediction_data: dictionary {string: list[float]}
        The prediction data a dictionary with keys being the name of the input server, and the
        values are a list that correspond to the distance in angstroms of the predicted
        amino acid compared to an unknown ground truth

    pathToSave: string
        This is a string representation to the path to the save folder, this holds the output.txt file

    target_name: string
        The name of the input target, taken from the input data folder names


    Return:
    ------------
    None:
        This method just writes to an output text file


    '''

    '''
    Header format:
    PFRMAT QA
    TARGET T0999
    AUTHOR 1234-5678-9000
    REMARK Error estimate is CA-CA distance in Angstroms
    METHOD Description of methods used
    MODEL 1
    '''

    with open(join(pathToSave, f'{target_name}.txt'), 'w+') as f:
        #set up the header
        f.write('PFRMAT\tQA\n')
        f.write(f'TARGET\t{target_name}\n')
        f.write(f"AUTHOR\tYOUR_NAME_HERE\n")
        f.write("REMARK\tNone\n")
        f.write(f"METHOD\tSVR\n")
        f.write("MODEL\t1\n")
        f.write("QMODE\t2\n")

        #write results
        for server_name_unformatted, server_predictions in prediction_data.items():
            server_name_formatted = server_name_unformatted.split('.')[0]
            # prediction_string = ' '.join([str(round(pred, 3)) for pred in server_predictions])
            prediction_string = ''
            i = 0
            for prediction in server_predictions:
                prediction_string += str(round(prediction, 3)) + " "
                if i % 25 == 0 and i != 0:
                    prediction_string += "\n"
                i += 1
            f.write(f"{server_name_formatted}\t{prediction_string}")
            f.write("\n")

        #write ending
        f.write("END\n")



def main(pathToInput, pathToSave):
    start = timer()
    pathToModel = './model/LocalQA_rbf_1_1.mdl'
    
    pattern = re.compile(r"T\d{4}[a-zA-Z]*[0-9]*")
    target_name = re.search(pattern, pathToInput)

    if target_name is not None: 
        target_name = str(target_name[0])
    else: 
        target_name = 'Target'

    #create save folder
    create_folder(pathToSave)

    #make the data the proper input format for the model
    model_input_path = preprocess_input(pathToInput, pathToSave)
    print("Input data created...")

    #load input data
    #remove when done testing
    # model_input_path = './TEST_OUT/tmp/ZoomQA_Input/step/T1096/'
    input_data = load_input_data(model_input_path)
    print("Input data loaded...")

    #load models
    model = load_model(pathToModel)

    # make predictions
    target_predictions = make_predictions(model, input_data)

    #write predictions
    write_predictions(target_predictions, pathToSave, target_name)
    print(f"Prediction saved to {pathToSave}")
    #remove tmp folder
    print("Cleaning up...")
    folder_to_remove = join(pathToSave, 'tmp')
    os.system(f'rm -rf {folder_to_remove}')
    end = timer()
    total_t = end-start
    print(f"Prediction complete, elapsed time: {total_t}")


def create_folder(pathToFolder):
    '''
    Method to create folder if does not exist, pass if it does exist,
    cancel the program if there is another error (e.g writing to protected directory)
    '''
    try:
        os.mkdir(pathToFolder)
    except FileExistsError:
        print(f"{pathToFolder} already exists...")
    except:
        print(f"Fatal error making {pathToFolder}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print('Not enough arguments... example command: ')
        print(f'python {sys.argv[0]} /path/To/Input/folder/ /path/to/output/save')
        sys.exit()

    print(ZOOMQA)
    #sys.exit()
    pathToInput = sys.argv[1]
    pathToSave = sys.argv[2]

    main(pathToInput, pathToSave)
