# Stride has 7 classes:
#   H	    Alpha helix
#   G	    3-10 helix
#   I	    PI-helix
#   E	    Extended conformation
#   B or	b   Isolated bridge
#   T	    Turn
#   C	    Coil (none of the above)
#
# However, eight types of secondary structures are too many for the existing
# methods of secondary structure prediction. Instead usually only three states
# are predicted: helix(H), extended(beta-sheet)(E) and coil(C). There are many
# different methods to translate the eight-letter DSSP alphabet into the 
# three-letter code. The translation used in the CASP experiment is as follows:
#   H,G,I -> H
#   E,B   -> E
#   T,S,C -> C
#
# Stride Output format
# ASG    Detailed secondary structure assignment
#    Format:  6-8  Residue name
#             10-10 Protein chain identifier
#	      12-15 PDB	residue	number
#	      17-20 Ordinal residue number
#	      25-25 One	letter secondary structure code
#	      27-39 Full secondary structure name
#	      43-49 Phi	angle
#	      53-59 Psi	angle
#	      65-69 Residue solvent accessible area


import sys
import os
import numpy
from os import path
import Bio.PDB
from Bio.PDB import *
import math
from PeptideBuilder import Geometry
import PeptideBuilder
from os import listdir
try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
import json

resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
	    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
	    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
	    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

def extract_backbone_model(pdb_path):
    parser=PDBParser()
    structure=parser.get_structure('sample', pdb_path)
    model=structure[0]
    chain=model['A']
    prev="0"
    N_prev="0"
    CA_prev="0"
    CO_prev="0"
    ##O_prev="0"
    prev_res=""
    rad=180.0/math.pi

    result = dict() 
    result['CA_C_N_angle'] = []
    result['C_N_CA_angle'] = []
    result['CA_N_length'] = []
    result['CA_C_length'] = [] 
    result['peptide_bond'] = []
    result['psi_im1'] = []
    result['omega'] = []
    result['phi'] = []
    result['CA_N_length'] = []
    result['CA_C_length'] = []
    result['N_CA_C_angle'] = []
    ### now first print the headers, please filter them when you really load the file ###
    #headers = "# residue CA_C_N_angle C_N_CA_angle CA_N_length CA_C_length peptide_bond psi_im1 omega phi CA_N_length CA_C_length N_CA_C_angle\n"
    
    for res in chain:
        if(res.get_resname() in resdict.keys()):
            geo=Geometry.geometry(resdict[res.get_resname()])
            if(prev=="0"):
                 N_prev=res['N']
                 CA_prev=res['CA']
                 C_prev=res['C']
                 ##O_prev=res['O']
                 prev="1"
            else:
                 n1=N_prev.get_vector()
                 ca1=CA_prev.get_vector()
                 c1=C_prev.get_vector()
                 ##o1=O_prev.get_vector()

                 ##O_curr=res['O']
                 C_curr=res['C']
                 N_curr=res['N']
                 CA_curr=res['CA']

                 ##o=O_curr.get_vector()
                 c=C_curr.get_vector()
                 n=N_curr.get_vector()
                 ca=CA_curr.get_vector()

                 geo.CA_C_N_angle=calc_angle(ca1, c1, n)*rad
                 geo.C_N_CA_angle=calc_angle(c1, n, ca)*rad
                 geo.CA_N_length= CA_curr-N_curr
                 geo.CA_C_length= CA_curr-C_curr
                 geo.peptide_bond= N_curr-C_prev

                 psi= calc_dihedral(n1, ca1, c1, n) ##goes to current res
                 omega= calc_dihedral(ca1, c1, n, ca) ##goes to current res
                 phi= calc_dihedral(c1, n, ca, c) ##goes to current res

                 geo.psi_im1=psi*rad
                 geo.omega=omega*rad
                 geo.phi=phi*rad

                 geo.CA_N_length= CA_curr - N_curr
                 geo.CA_C_length= CA_curr - C_curr
                 ##geo.C_O_length= C_curr - O_curr

                 geo.N_CA_C_angle= calc_angle(n, ca, c)*rad
                 ##geo.CA_C_O_angle= calc_angle(ca, c, o)*rad

                 ##geo.N_CA_C_O= calc_dihedral(n, ca, c, o)*rad

                 N_prev=res['N']
                 CA_prev=res['CA']
                 C_prev=res['C']
                 ##O_prev=res['O']
            # now add angles to result
            result['CA_C_N_angle'].append(str(geo.CA_C_N_angle))
            result['C_N_CA_angle'].append(str(geo.C_N_CA_angle))
            result['CA_N_length'].append(str(geo.CA_N_length))
            result['CA_C_length'].append(str(geo.CA_C_length))
            result['peptide_bond'].append(str(geo.peptide_bond))
            result['psi_im1'].append(str(geo.psi_im1))
            result['omega'].append(str(geo.omega))	
            result['phi'].append(str(geo.phi))
            result['CA_N_length'].append(str(geo.CA_N_length))
            result['CA_C_length'].append(str(geo.CA_C_length))
            result['N_CA_C_angle'].append(str(geo.N_CA_C_angle))
    return result

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            if row <= col:
                answer[row, col] = calc_residue_dist(residue_one, residue_two)
            else:
                answer[row, col] = answer[col, row]
    return answer

def extract_contacts_model(pdb_path):
    parser=PDBParser()
    structure=parser.get_structure('sample', pdb_path)
    model=structure[0]
    chain=model['A']
    dist_matrix = calc_dist_matrix(model["A"], model["A"])
    return dist_matrix

#   H,G,I -> H
#   E,B   -> E
#   T,S,C -> C
def convert_8_to_3(inputss):
    if inputss == 'H' or inputss == 'G' or inputss == 'I':
       return 'H'
    if inputss == 'T' or inputss == 'S' or inputss == 'C':
       return 'C'
    else:
       return 'E'
def extract_ss(pdb_path, tool_path):
    out = os.popen(tool_path + " "+pdb_path).read().split("\n")
    #print(out)
    ss = []
    aa = []
    sol = []
    for line in out:
        tem = line.strip().split()
        if len(tem) < 10:
            continue
        if tem[0] != "ASG":
            continue 
        tem[1] = tem[1].upper() 
        if tem[1] not in resdict:
            print("Warning, the residue "+tem[1]+" is not recognized, skip it")
            continue
        aa.append(resdict[tem[1]])
        ss.append(convert_8_to_3(tem[5].upper()))
        sol.append(float(tem[9]))
    return (ss, aa, sol)


if __name__ == "__main__":
    if(len(sys.argv)<3):
       print("Version 3, we save json result for each target, instead of the whole CASP dataset, so each of them would be much smaller!")
       print("This script would import all information for one folder like CASP5, it will generate all information and save it to a json file")
       print("This script need three inputs, the first is the Stride exe file, the second is directory for all targets like CASP5, the second is the output directory for json file. \n")
       print("For example:\n")
       print("python "+sys.argv[0]+" ./stride ../test/CASP5 ../test/json_CASP5")
       sys.exit(0)
    strideTool = sys.argv[1] 
    inputDir = sys.argv[2]
    #LGAScoreDir = sys.argv[3]
    outputDir = sys.argv[3]
    #outFilePath = outputDir+"/"+inputDir.split('/')[-1]+".json"
    #DB = dict()
    try:
       os.stat(outputDir)
    except:
       os.mkdir(outputDir)
    for targetName in listdir(inputDir):
        outFilePath = outputDir+"/"+inputDir.split('/')[-1]+"_"+targetName+".json"
        DB = dict()
        checkRun = outFilePath+".tmpRun"
        if os.path.exists(checkRun):
            continue     # someone already run this script 
        else:
            fh1 = open(checkRun,"w")
            fh1.write("I am running ...")
            fh1.close()
        #GDTdict = loadGDT(GDTallPath)  # we don't know the GDT score for prediction
        for modelName in listdir(inputDir+"/"+targetName):
            pdbPath = inputDir+"/"+targetName+"/"+modelName
            print("processing "+pdbPath)
            uniqKey = targetName + ":" + modelName   # don't use tuple as key for the dictionary for keep all of our information, use X:X because we use NMA tool to expand models, there would be duplicated targetName, but those two would be unique, only the native casp pdb would overlap, but I guess we could keep one of them, it's fine
            F_GDT = -1
            # extract all secondary structure, amino acid, and solvent accessibility
            (F_ss, F_aa, F_sol) = extract_ss(pdbPath, strideTool)
            F_localQA = []
            for j in range(len(F_ss)):
               F_localQA.append(-1)     # we don't know the local QA score, put -1 
            try:
               F_dis_matrix = extract_contacts_model(pdbPath).tolist()
            except:
               print("Error to extract contact map, use 0 "+pdbPath)
               F_dis_matrix = numpy.zeros((len(F_ss), len(F_ss)), numpy.float)
            # now we need to get all angles information, and we are done for this model!
            try:
                F_backboneAngles = extract_backbone_model(pdbPath)
            except:
                print("This model "+pdbPath+" may only contains CA, we skip those kind of models for now")
                continue
            print("Adding ...") 
            DB[uniqKey] = dict()
            DB[uniqKey]['GDT'] = F_GDT
            DB[uniqKey]['localQA'] = F_localQA
            DB[uniqKey]['ss'] = F_ss
            DB[uniqKey]['aa'] = F_aa
            DB[uniqKey]['sol'] = F_sol
            DB[uniqKey]['ContactMap'] = F_dis_matrix
            DB[uniqKey]['Angles'] = F_backboneAngles
            #print(pdbPath)
            #print(F_GDT)
            #print(F_localQA)
            #print(F_ss)
            #print(F_aa)
            #print(F_sol)
            #print(F_dis_matrix)
            #print(F_backboneAngles)
            #print(len(F_localQA))
            #print(len(F_ss))
            #print(len(F_aa))
            #print(len(F_sol))
            #print(len(F_dis_matrix))
            #print(len(F_backboneAngles['CA_C_N_angle']))
            #sys.exit(0)
        # now store everything in a json file
        #print("Now we have the following:")
        #for each in DB:
        #   print(each)
        #   print(DB[each]) 
        with open(outFilePath, 'w') as fp:
           json.dump(DB, fp)
        #pickle.dump(str(DB), fp, protocol=pickle.HIGHEST_PROTOCOL)
    # now you could load it back using : json.load
