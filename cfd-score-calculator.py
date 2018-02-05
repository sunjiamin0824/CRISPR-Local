#Calculates the Cutting Frequency Determination score
#Requirements: 1. Pickle file with mismatch scores in working directory
#              2. Pickle file containing PAM scores in working directory 
#Input: 1. 23mer WT sgRNA sequence
#       2. 23mer Off-target sgRNA sequence
#Output: CFD score
import pickle
import argparse
import re
import numpy as np
import os

def get_parser():
    parser = argparse.ArgumentParser(description='Calculates CFD score')
    parser.add_argument('--input', type=str, help='WT and Off-target 23mer sgRNA sequence file.')
    parser.add_argument('--output', type=str, help='CFD score result file.')
    return parser

#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

#Unpickle mismatch scores and PAM scores
def get_mm_pam_scores():
    try:
        mm_scores = pickle.load(open('mismatch_score.pkl','rb'))
        pam_scores = pickle.load(open('pam_scores.pkl','rb'))
        return (mm_scores,pam_scores)
    except: 
        raise Exception("Could not find file with mismatch scores or PAM scores")

#Calculates CFD score
def calc_cfd(wt,sg,pam):
    mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score*= mm_scores[key]
    score*=pam_scores[pam]
    return (score)

if __name__ == '__main__':
	args = get_parser().parse_args()
	mm_scores,pam_scores = get_mm_pam_scores()
	input = args.input
	output = args.output
	i = open(input,'r')
	o = open(output,'w')
	for line in i:
		line = line.replace('\n','')
		line_list=line.split('\t')
		if line_list[0]==line_list[4] and line_list[7] == 0:continue
		off = line_list[6]
		wt = line_list[2]
		pam = off[-2:]
		sg = off[:-3]
		cfd_score = calc_cfd(wt,sg,pam)
		CFD = '%.4f'%(cfd_score)
		o.write(line + '\t' + CFD + '\n')
	i.close()
	o.close()
	os.remove(input)



