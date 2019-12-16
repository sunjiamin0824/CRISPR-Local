#Calculates the Rule set 2 score for the given 30-mer
#Input: 1. 30mer sgRNA+context sequence, NNNN[sgRNA sequence]NGGNNN
#       2. Amino acid cut position, for full model prediction only
#       3. Percent peptide, for full model prediction only
#Output: Rule set 2 score

import pandas as pd
import csv, argparse, sys
import pickle
import model_comparison
import os
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, help='input 30-mer sgRNA fasta file')
    parser.add_argument('--output',	type=str, help='output 20-mer sgRNA-scored fasta file')
    parser.add_argument('--aa-cut', type=int, default=-1, help='Amino acid cut position of sgRNA')
    parser.add_argument('--per-peptide', type=float, default=-1, help='Percentage of protein cut by sgRNA')
    return parser
if __name__ == '__main__':
    args = get_parser().parse_args()
    input = args.input
    output = args.output
    aa_cut = args.aa_cut
    per_peptide = args.per_peptide
    model_file_1 = 'Rule_Set_2_scoring_v1/saved_models/V3_model_nopos.pickle'
    model_file_2 = 'Rule_Set_2_scoring_v1/saved_models/V3_model_full.pickle'
    if (aa_cut == -1) or (per_peptide == -1):
        model_file = model_file_1
    else:
        model_file = model_file_2
    try:
        with open(model_file, 'rb') as f:
            model= pickle.load(f)    
    except:
        raise Exception("could not find model stored to file %s" % model_file)
    i = open(input,'r')
    o = open(output,'w')
    seq={}
    for line in i:
        if line.startswith('>'):
                name = line.replace('\n','')
        else:
                seq[name]=line.replace('\n','')
                score = model_comparison.predict(seq[name], aa_cut, per_peptide, model=model)
                o.write(name + '#%.6f'% (score) + '\n' + seq[name][4:24] + '\n')
    i.close()
    o.close()
    os.remove(input)
