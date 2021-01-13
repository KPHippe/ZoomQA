import os
import sys 
import json
import argparse
from os.path import isfile, isdir, join

import numpy as np 

def get_gdt(lqa_scores): 

    newQA = np.array(lqa_scores)

    return (1/(1+newQA*newQA/12)).mean()


def get_gdt_scores(scores): 
    gdt_lqa_scores = {}

    for server, lqa_scores in scores.items(): 
        gdt_score = get_gdt(lqa_scores)
        gdt_lqa_scores[server] = lqa_scores
        gdt_lqa_scores[server].insert(0, gdt_score)

    return gdt_lqa_scores

def process_file(file_path): 

    def _parse_casp_file(data):
        #remove header 
        new_lines = [i for i, e in enumerate(data) if e == '\n']
        data = data[new_lines[6]+1:]
        #remove footer 
        data = data.replace("\nEND\n", '')

        data = data.split('\n')
        #print(data)
        res_dict = {}
        target_name = None
        scores = []
        for val in data: 
            if '\t' in val:
                if len(scores) != 0: 
                    res_dict[target_name] = scores
                    scores = []
                target_name = val.split('\t')[0]
                scores.extend([float(e) for e in val.strip().split('\t')[-1].split(" ")])
            else: 
                if val == '':
                    continue
                scores.extend([float(e) for e in val.strip().split(" ")]) 
        
        return res_dict

    raw_data = open(file_path, 'r').read() 
    scores = _parse_casp_file(raw_data)
    gdt_lqa_scores = get_gdt_scores(scores)
    _write_to_file(file_path, gdt_lqa_scores)
    print(f"Processed {file_path}")

def _write_to_file(file_path, scores): 
    target_name = file_path.split('/')[-1].split('.')[0]
    with open(file_path, 'w+') as f:
        #set up the header
        f.write('PFRMAT\tQA\n')
        f.write(f'TARGET\t{target_name}\n')
        f.write(f"AUTHOR\tYOUR_NAME_HERE\n")
        f.write("REMARK\tNone\n")
        f.write(f"METHOD\tSVR\n")
        f.write("MODEL\t1\n")
        f.write("QMODE\t2\n")

        #write results
        for server_name, server_predictions in scores.items():
            # prediction_string = ' '.join([str(round(pred, 3)) for pred in server_predictions])
            prediction_string = ''
            i = 0
            for prediction in server_predictions:
                prediction_string += str(round(prediction, 3)) + " "
                if i % 25 == 0 and i != 0:
                    prediction_string += "\n"
                i += 1
            f.write(f"{server_name}\t{prediction_string}")
            f.write("\n")

        #write ending
        f.write("END\n")


def main(input=None): 
    if input is None: 
        print('No input provided...')
        sys.exit()


    all_file_paths = _find_out_files(input)
    
    for file in all_file_paths: 
        process_file(file)

def _find_out_files(input): 
    if isfile(input): 
        return [input]
    all_paths = []
    for path, dirs, files in os.walk(input):
        for file in files:
            if '.txt' in file[-4:]:
                all_paths.append(join(path, file))

    return all_paths


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to add global scores to the lqa output of ZoomQA')

    parser.add_argument('-i', '--input',
            help='The path to the input directory/file/ to add global qa to. This will determine the filestructure and overwite existign files with the added GDT')


    args=parser.parse_args()

    main(**vars(args))
