import os
from os import listdir
from os.path import isfile, join
import sys
import numpy
from os import path



if __name__ == "__main__":
    if(len(sys.argv)<4):
       print("modified by Dr. Cao on May 22, 2020, it would handle one target each time")
       print("Modified by Renzhi on 7/30/2019, it can handle parallel because there are too many NMA_PISCES_normal")
       print("This script is going to add chain A to all pdbs because the LGA software may need it")
       print("For example:\nYou should directly use the wrapper 3 program for casp targets unless you just want to use this one, uncomment to see the examples")
       print("python "+ sys.argv[0]+" "+"assist_add_chainID_to_one_pdb.pl ../../data/NMA_PISCES_normal_20/0 ../../result/chainAddedAll/NMADecoys_ALL/NMA_PISCES_0")
       print("python "+ sys.argv[0]+" "+"assist_add_chainID_to_one_pdb.pl ../../data/NMA_PISCES_normal_20/1 ../../result/chainAddedAll/NMADecoys_ALL/NMA_PISCES_1")
       #print("python "+sys.argv[0]+" "+"assist_add_chainID_to_one_pdb.pl "+"../../data/CASP_filtered_same_seq/CASP12"+" ../../result/chainAddedAll/CASP_ALL/CASP12")
       #print("python "+sys.argv[0]+" "+"assist_add_chainID_to_one_pdb.pl ../../data/NMADecoys_CASP/CASP12 ../../result/chainAddedAll/NMADecoys_ALL/CASP12")
       sys.exit(0)
    tool_add = sys.argv[1] 
    dir_input = sys.argv[2]
    dir_out = sys.argv[3]
    
    try:
        os.stat(dir_out)
    except:
        os.mkdir(dir_out)

    # now list all folders and try to add chains to new output folder
    inputTarget = dir_input
    outTarget = dir_out
    for PDBfile in listdir(inputTarget):
       inPDB = inputTarget+"/"+PDBfile
       outPDB = outTarget+"/"+PDBfile
       os.system("perl "+tool_add+" "+inPDB + " " +outPDB)
