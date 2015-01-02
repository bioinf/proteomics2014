#!/bin/env python2.7
# -*- coding: utf-8 -*-
"""
dssp - find secondary structure of protein by its tertiary structure

Usage:
  dssp.py -i <input_pdb_file> -o <output_tracks_file>

Options:
  -h --help                Show this screen.
  --version                Show version.
  -i <input_pdb_file>      File with optical reads.
  -o <output_tracks_file>  File for false reads.
 """

import sys

modules = ["os", "docopt", "Bio"]
exit_flag = False
for module in modules:
    try:
        __import__(module)
    except ImportError:
        exit_flag = True
        sys.stderr.write("Error: Python module " + module + " is not installed.\n")
        
if exit_flag:
    sys.exit("Some required Python modules are not installed. Exit.")

from docopt import docopt
import os
from Bio.PDB import *
from math import sqrt
import string
import itertools

# Consts
max_line_length = 50

def is_non_zero_file(fpath):  
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False


def calc_H_coords(atom_C, atom_O, atom_N):
    # NH distance is approximately 1 A (Creighton, 1993)
    # NH || CO
    C_x, C_y, C_z = atom_C.get_coord()
    O_x, O_y, O_z = atom_O.get_coord()
    N_x, N_y, N_z = atom_N.get_coord()
    OC_x = C_x - O_x
    OC_y = C_y - O_y
    OC_z = C_z - O_z
    OC_len = sqrt(OC_x**2 + OC_y**2 + OC_z**2)
    OC_unit_x = OC_x / OC_len
    OC_unit_y = OC_y / OC_len
    OC_unit_z = OC_z / OC_len
    H_x = N_x + OC_unit_x
    H_y = N_y + OC_unit_y
    H_z = N_z + OC_unit_z
    return H_x, H_y, H_z


def parall_gen():
    letters = itertools.cycle(string.ascii_uppercase)
    for letter in letters:
        yield letter


def antiparall_gen():  
    letters = itertools.cycle(string.ascii_lowercase)
    for letter in letters:
        yield letter


if __name__ == '__main__':
    arguments = docopt(__doc__, version='dssp 0.1')
    input_pdb_filename = arguments["-i"]
    output_tracks_filename = arguments["-o"]

if not is_non_zero_file(input_pdb_filename):
    sys.exit(input_pdb_filename + ": no such file or directory. Exit.")

aa_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', \
           'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', \
           'TYR', 'VAL']
aa_set = frozenset(itertools.chain(aa_list))

parser = PDBParser()
structure = parser.get_structure('luc', input_pdb_filename)
ppb = PPBuilder()
pp = ppb.build_peptides(structure)
sequence = pp[0].get_sequence()

print "Protein sequence extracted"

# Find H-bonds

print "Find H-bonds...",
residues_list = Selection.unfold_entities(structure, 'R')
residues_list_len = len(residues_list)
hbond = [[False for x in range(residues_list_len)] \
         for x in range(residues_list_len)]
res_num = 0
for i, residue in enumerate(residues_list):
    resname1 = residue.get_resname()
    if resname1 not in aa_set:
        break
    res_num += 1
    atom_C_1 = residue['C']
    atom_O_1 = residue['O']
    C1_x, C1_y, C1_z = atom_C_1.get_coord()
    O1_x, O1_y, O1_z = atom_O_1.get_coord()
    for j, residue in enumerate(residues_list):
        if j == i:
            continue
        resname2 = residue.get_resname()
        if resname2 not in aa_set:
            break
        atom_N = residue['N']
        N_x, N_y, N_z = atom_N.get_coord()
        if j < residues_list_len - 1 and \
           residues_list[j + 1].get_resname() in aa_set:
            next_residue = residues_list[j + 1]
            atom_C_2 = next_residue['C']
            atom_O_2 = next_residue['O']
            dist_ON = atom_O_1 - atom_N
            dist_CN = atom_C_1 - atom_N
            H_x, H_y, H_z = calc_H_coords(atom_C_2, atom_O_2, atom_N)
            dist_CH = sqrt((C1_x - H_x)**2 + (C1_y - H_y)**2 + (C1_z - H_z)**2)
            dist_OH = sqrt((O1_x - H_x)**2 + (O1_y - H_y)**2 + (O1_z - H_z)**2)
            f = 332 # [A * kcal / e^2]
            q1 = 0.42 # [e]
            q2 = 0.2 # [e]
            energy = f * q1 * q2 * (1/dist_ON + 1/dist_CH - 1/dist_OH - 1/dist_CN)
            if energy < -0.5: 
                hbond[i][j] = True
        else:
            break
print "finished"
          
# Describe secondary structures
# Notation from original DSSP paper is used

# 1. Turns and bridges

turns = [[' ' for i in range(res_num)] for j in range(3)]
bridges = [' ' for i in range(res_num)]

# Check n-turns (n = 3, 4, 5)

print "Find turns...",
for i in range(res_num):
    for n in [3, 4, 5]:
        if hbond[i][i + n]:
            if turns[n - 3][i] == '>' or turns[n - 3][i] == '<':
                turns[n - 3][i] = 'X'
            else:
                turns[n - 3][i] = '>'
            if turns[n - 3][i + n] == '>' or turns[n - 3][i + n] == '<':
                turns[n - 3][i + n] = 'X'
            else:
                turns[n - 3][i + n] = '<'
            for k in range(i + 1, i + n):
                if turns[n - 3][k] == ' ':
                    turns[n - 3][k] = str(n)
                    
print "finished"

# Check bridges

print "Find bridges...",
parall_generator = parall_gen()
antiparall_generator = antiparall_gen()
for i in range(res_num):
    for j in range(i + 1, res_num):
        if i > 0 and hbond[i - 1][j] and hbond[j][i + 1] or \
           j > 0 and hbond[j - 1][i] and hbond[i][j + 1]:
               letter = next(parall_generator) if i == 0 or bridges[i - 1] == ' ' \
                                         else bridges[i - 1] 
               bridges[i] = letter # 'P' stands for 'Parallel'
               bridges[j] = letter #'P'
               #print letter
        if hbond[i][j] and hbond[j][i] or \
           i > 0 and j > 0 and hbond[i - 1][j + 1] and hbond[j - 1][i + 1]:
               letter = next(antiparall_generator) if i == 0 or bridges[i - 1] == ' ' \
                                             else bridges[i - 1]
               bridges[i] = letter # 'A' stands for 'Antiparallel'
               bridges[j] = letter #'A'
               #print letter
print "finished"

# 2. Helixes and sheets

helixes = [[' ' for i in range(res_num)] for j in range(3)]

# Check n-helixes (n = 3, 4, 5)

print "Find helixes...",
indexes = [0, 0, 0]
while indexes[0] < res_num or indexes[1] < res_num or indexes[2] < res_num:
    for n in [3, 4, 5]:
        i = indexes[n - 3]
        if i > 0 and i < res_num and hbond[i - 1][i + n - 1] and hbond[i][i + n]:
            indexes[n - 3] += n
            for j in range(i, i + n): # mark helix from i to i + n - 1
                if n == 4:
                    helixes[n - 3][j] = 'H'
                elif n == 3:
                    helixes[n - 3][j] = 'G'
                else: # n == 5
                    helixes[n - 3][j] = 'I'
        else:
            indexes[n - 3] += 1
print "finished"

# Check sheets

print "Find sheets...",
sheets = [' ' for x in range(res_num)]
for i in range(res_num):
    if bridges[i] != ' ': # some bridge here
        if i == 0 and bridges[i + 1] == ' ' or \
           i == res_num - 1 and bridges[i - 1] == ' ' or \
           bridges[i - 1] == ' ' and bridges[i + 1] == ' ':
               sheets[i] = 'B' # single bridge
        else:
            sheets[i] = 'E' # 'extended' sequence of bridges
print "finished"

# Summary

print "Build secondary structure summary...",
summary = [' ' for i in range(res_num)]
for i in range(res_num):
    if sheets[i] == 'B':
        summary[i] = 'B'
    elif sheets[i] == 'E':
        summary[i] = 'E'
    elif helixes[1][i] == 'H':
        summary[i] = 'H'
    elif helixes[0][i] == 'G':
        summary[i] = 'G'
    elif helixes[2][i] == 'I':
        summary[i] = 'I'
    elif turns[0][i] != ' ' or turns[1][i] != ' ' or turns[2][i] != ' ':
        summary[i] = 'T'
print "finished"
     
# Output tracks

print "Output secondary structure to", output_tracks_filename + "...",
with open(output_tracks_filename, 'w') as dst:
    index = 0;
    while index < res_num:
        dst.write('Bridges: ' + ''.join(bridges[index:index + max_line_length]) + '\n')
        dst.write('5-turns: ' + ''.join(turns[2][index:index + max_line_length]) + '\n')
        dst.write('4-turns: ' + ''.join(turns[1][index:index + max_line_length]) + '\n')
        dst.write('3-turns: ' + ''.join(turns[0][index:index + max_line_length]) + '\n')
        dst.write('Summary: ' + ''.join(summary[index:index + max_line_length]) + '\n')
        dst.write('AA  seq: ' + ''.join(sequence[index:index + max_line_length]) + '\n')
        index += max_line_length;
        if index < res_num:
            dst.write('-'*(max_line_length + 9) + '\n');
print "finished"
