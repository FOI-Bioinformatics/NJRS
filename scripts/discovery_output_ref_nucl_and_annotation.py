#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from Bio import SeqIO
from matplotlib import pyplot as plt
import argparse

software_description = """
From an inputfile and column of reference positions, output table with reference nucleotide at that position and annotation
"""


### Parse input arguments
## setup
argparser = argparse.ArgumentParser(description=software_description)

argparser.add_argument('-i','--input',required=True,help='Path to input table')
argparser.add_argument('--ref_pos_col',required=False,type=int,default=0,help='Column of reference positions, 0-based (default: column 0)')

argparser.add_argument('-r','--reference',required=True,help='Path to reference fasta')
argparser.add_argument('-g','--gff',required=True,help='Path to reference annotation GFF file')

argparser.add_argument('--split_input_by_comma',required=False,default=None,action='store_true',help='If toggled, will assume input is a comma-separated list on first row after header')

argparser.add_argument('-o','--output',required=True,help='Path to output file')
##/
## parse input
args = argparser.parse_args()

input_table = args.input
ref_pos_col = args.ref_pos_col

reference_fasta_path = args.reference
reference_gff_path = args.gff

split_input_by_comma = args.split_input_by_comma

output_file = args.output
##/
###/


####### SCRIPT START
### Load ref seq and gff
## Fasta
print('Importing reference-sequence')
ref_seqs = {}
for record in SeqIO.parse(reference_fasta_path, "fasta"):
    print("Sequence ID:", record.id)
    print("Description:", record.description)
    
    ref_seqs[record.id] = record.seq

ref_name = list(ref_seqs.keys())[0]
ref_seq = ref_seqs[ref_name]

print(f'Ref seq identified as {ref_name} with length {len(ref_seq)}')
##/
## GFF
print('Importing reference-annotation')
gffs = []
gff_types_count = {}
with open(reference_gff_path) as f:
    for line in f:
        # skip if header
        if line[0] == '#':
            continue
        #/
        # parse line
        line = line.strip('\n')
        if not line: continue
        line = line.split('\t')
        #/
        # parse data
        rname,source_,gff_type,rstart,rend,score,strand,frame,attributes = line
        rstart = int(rstart)
        rend = int(rend)
        #/
        # parase attributes-column
        attributes_field = {}
        for chunk in attributes.split(';'):
            if chunk == '': continue # skip empty
            key,val = chunk.split('=',1)
            attributes_field[key] = val
        #/
        # skip if "genome region"
        if gff_type == 'region':
            print(f'Skipping "Region" entry with size {rend-rstart}')
            continue
        #/
        ## Compile save
        tmp_save = {'rname':rname,'rcoords':[rstart,rend],'strand':strand,'type':gff_type}
        tmp_save['ID'] = attributes_field['ID']
        
        # get a name of entry, such as "Name", "gene", etc
        ids_to_check_name = ('Name','gene','locus_tag','Note')
        tmp_save['name'] = None
        for CK in ids_to_check_name:
            if CK in attributes_field:
                tmp_save['name'] = attributes_field[CK]
                break
        if tmp_save['name'] == None:
            print(f'WARNING: Was unable to assign a name for entry\n{line}\ntested keys: {ids_to_check_name}')
            print('Moving on...')
        #/
        
        # get a description of what this entry does, such as "product"
        ids_to_check_description = ('description','product','gene_biotype','regulatory_class')
        tmp_save['description'] = None
        for CK in ids_to_check_description:
            if CK in attributes_field:
                tmp_save['description'] = attributes_field[CK]
                break
        if tmp_save['description'] == None:
            print(f'WARNING: Was unable to assign a description for entry\n{line}\ntested keys: {ids_to_check_description}')
            print('Moving on...')
        ##/
        # save
        gffs.append(tmp_save)
        
        if not gff_type in gff_types_count:         gff_types_count[gff_type] = 0
        gff_types_count[gff_type] += 1
        #/
print(f'Imported {len(gffs)} GFF-entries')
print('Number of entries per type:')
for gff_type,count in gff_types_count.items():
    print(f'\t{gff_type} \t{count}')
##/
###/

### Read input, get positions
ref_poses = set()
with open(input_table,'r') as f:
    for enum,line_raw in enumerate(f):
        # parse line
        line_raw = line_raw.strip('\n')
        if not line_raw: continue
        # skip if header
        if line_raw[0] == '#': continue
        #/
        line = line_raw.split('\t')
        #/
        ## WEDGED CODE: check if user wanted to parse a comma-separated list of positions
        if split_input_by_comma:
            line = line_raw.split(',')
            for pos in line:
                try:
                    ref_pos = int(pos)
                    ref_poses.add(ref_pos)
                except:
                    print(f'Unable to parse line with enum={enum}: {line}')
                    print('This is expected if there is a header in the file. Skipping this line and moving on...')
                    continue
            continue # if line was parsed as "comma-separated", do not attempt to process it in the normal way
        ##/
        # parse data
        try:
            ref_pos = int(line[ref_pos_col])
        except:
            print(f'Unable to parse line with enum={enum}: {line}')
            print('This is expected if there is a header in the file. Skipping this line and moving on...')
            continue
        #/
        # save
        ref_poses.add(ref_pos)
        #/
print(f'Imported {len(ref_poses)} reference positions')
###/

### Find position overlaps with GFF
print('Finding GFF overlaps to positions')
ref_poses_gff_ovlps = {}
for ref_pos in ref_poses:
    # get ovlps
    gff_ovlps = []
    for gff in gffs:
        if ref_pos >= gff['rcoords'][0] and ref_pos <= gff['rcoords'][-1]:
            gff_ovlps.append(gff)
    #/
    # compile save (as stated in the cell)
    gff_ovlps_save = {'types':[],'names':[],'descriptions':[]}
    for gff in gff_ovlps:
        gff_ovlps_save['types'].append(gff['type'])
        gff_ovlps_save['names'].append(gff['name'])
        gff_ovlps_save['descriptions'].append(gff['description'])
    #/
    # save
    ref_poses_gff_ovlps[ref_pos] = gff_ovlps_save
    #/
###/

### Output
print(f'Writing output at {output_file}')
with open(output_file,'w') as nf:
    # write header
    header = ['ref_pos','ref_nucl','types','names','descriptions']
    nf.write('\t'.join(map(str,header))+'\n')
    #/
    # write rows
    for ref_pos,data in sorted(ref_poses_gff_ovlps.items(),key=lambda x: x[0]):
        ref_nucl = ref_seq[ref_pos]
        writeArr = [ref_pos,ref_nucl,','.join(data['types']),','.join(data['names']),','.join(data['descriptions'])]
        nf.write('\t'.join(map(str,writeArr))+'\n')
    #/
###/

### Finalize
print('All done!')
###/
