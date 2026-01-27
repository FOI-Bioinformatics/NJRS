#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import argparse
from Bio import SeqIO


### Parse input arguments
## setup
argparser = argparse.ArgumentParser()

argparser.add_argument('--basecounts_dir',required=True,help='Path to clade basecounts directory')
argparser.add_argument('--msa_slim',required=True,help='Path to MSA input file')
argparser.add_argument('--msa_slim_coords',required=True,help='Path to file with coordinates of input MSA file')
argparser.add_argument('-o','--outdir',required=True,help='Path to output directory')
argparser.add_argument('-r','--ref','--reference','--ref_fasta',required=True,help='Path to REF fasta')

argparser.add_argument('--ambig_output_cutoff',required=False,type=float,default=0.2,help='Maximum fraction of datasets with "ambiguous" allowed in MSA output (default: 0.2)')
##/
## parse input
args = argparser.parse_args()

basecounts_dir_path = args.basecounts_dir
msa_slim_path = args.msa_slim
msa_slim_coords_path = args.msa_slim_coords
ref_fasta = args.ref
output_dir = args.outdir

ambig_output_cutoff = args.ambig_output_cutoff
##/
###/


####### FUNCTIONS
# file open function (gzip and non-compressed)
def open_file(filepath):
    # Check if the file is gzipped by checking its extension or magic number
    if filepath.endswith('.gz'):
        # Open gzipped file
        return gzip.open(filepath, 'rt')  # 'rt' mode for reading text
    else:
        # Open regular file
        return open(filepath, 'r')
#/
# parse rows from basecounts
def parse_basecounts(basecounts_file_path,poses_to_get=set(),return_header=False):
    print(f'Parse basecounts: {basecounts_file_path}')
    pos_basecounts = {} # pos -> rows with basecounts
    basecounts_header = None
    with open_file(basecounts_file_path) as f:
        for enum,line in enumerate(f):
            # parse line
            line = line.strip('\n')
            line = line.split('\t')
            #/
            # check if header (save it once and skip)
            if enum == 0:
                if basecounts_header == None:
                    basecounts_header = line
                continue
            #/
            # parse data
            pos = line[0]
            #/
            # check if skip pos
            if poses_to_get and not pos in poses_to_get: continue
            #/
            # save
            pos_basecounts[pos] = line
            #/
    
    # define return
    tmp_return = pos_basecounts
    if return_header == True:
        tmp_return = [pos_basecounts,basecounts_header]
    #/
    
    return tmp_return
###/
#######/

####### SCRIPT START
### parse reference sequence
print('Importing reference-sequence')
ref_seqs = {}
for record in SeqIO.parse(ref_fasta, "fasta"):
    print("Sequence ID:", record.id)
    print("Description:", record.description)
    
    ref_seqs[record.id] = record.seq

ref_name = list(ref_seqs.keys())[0]
ref_seq = ref_seqs[ref_name]

print(f'Ref seq identified as {ref_name} with length {len(ref_seq)}')
###/

### Parse MSA seq
msa_seqs = {}
for record in SeqIO.parse(msa_slim_path, "fasta"):
    print("Sequence ID:", record.id)
    print("Description:", record.description)
    
    msa_seqs[record.id] = record.seq
###/

### Import coords to parse from MSA_slim_coords
ref_coords_to_import = {} # ref_coord -> msa_idx
msa_slim_coords_header = None
with open(msa_slim_coords_path,'r') as f:
    for enum,line in enumerate(f):
        # parse line
        line = line.strip('\n')
        #/
        # skip if header
        if enum == 0:
            if msa_slim_coords_header == None:
                msa_slim_coords_header = line
            continue
        #/
        # parse coords (assumes a comma-separated list on line 2)
        if enum == 1:
            for msa_idx,coord in enumerate(line.split(',')):
                if not coord in ref_coords_to_import:
                    ref_coords_to_import[coord] = msa_idx
        #/
        # bugcheck if more lines
        if enum >= 2:
            print(f'WARNING: Found additional lines in {msa_slim_coords_path}')
            print('This file is onyl expected to have x2 lines')
        #/
print(f'Parsed N={len(ref_coords_to_import)} reference sequence positions from MSA file')
###/

### Import basecounts for msa-ref-positions
print('Begin importing basecounts')
msa_coords_basecounts = {}
for basecounts_file in os.listdir(basecounts_dir_path):
    basecounts_file_path = basecounts_dir_path+'/'+basecounts_file
    sample_name = basecounts_file.strip('.gz')
    pos_basecounts,header = parse_basecounts(basecounts_file_path,poses_to_get=ref_coords_to_import,return_header=True)
    basecounts_header = header
    
    for pos,basecounts_raw in pos_basecounts.items():
        # parse array-representation of basecounts to dict-representation
        basecounts = {}
        for idx,entry in enumerate(basecounts_raw):
            # check if entry is a number/float
            if entry.replace('.','').isdigit():
                if int(entry) == float(entry): # convert natural number
                    entry = int(entry)
                else: # convert decimal number
                    entry = float(entry)
            #/
            basecounts[basecounts_header[idx]] = entry
        #/
        # init pos
        if not pos in msa_coords_basecounts:                msa_coords_basecounts[pos] = {}
        #/
        # save sample basecounts at pos
        msa_coords_basecounts[pos][sample_name] = basecounts
        #/
print('Basecounts imported')
###/

### Compute stability (save in dict entry)
print('Classifying stability for positions (ATCG counts only; others ignored or set to ND="Not Determined")')
msa_coords_basecounts_classified = {} # pos -> sample -> classification with data and basecounts
IDE_basecounts_with_N = []
for pos,samples_basecounts in msa_coords_basecounts.items():
    
    for sample,basecounts in samples_basecounts.items():
        # get expected nucleotide at position from MSA
        msa_nucl_idx = ref_coords_to_import[pos]
        expected_nucl = msa_seqs[sample][msa_nucl_idx]
        #/
        
        ## determine stability (based on ATGC only according to kerstin [ref chat message])
        stability_classification = None
        # For positions with coverage
        coverage = basecounts['coverage']
        #coverage_metric_to_use = coverage
        
        ATGC_coverage = 0 # use ATGC counts for coverage instead of anything else
        for letter in ('A','T','G','C',):
            ATGC_coverage += basecounts[letter]
        coverage_metric_to_use = ATGC_coverage
        # IDE
        if basecounts['N'] > 0:
            IDE_basecounts_with_N.append(basecounts)
        #/
        
        support_fraction = 'NA'
        if expected_nucl != '-':
            expected_nucl_count = basecounts[expected_nucl]
            
            support_fraction = expected_nucl_count / coverage_metric_to_use
            
            # test rules <=20 cov
            if coverage_metric_to_use <= 20:
                if support_fraction < 0.90: # classify ambiguous when below this fraction of support
                    stability_classification = 'ambiguous'
                else:
                    stability_classification = 'stable'
            #/
            # test rules >20 cov
            elif coverage_metric_to_use > 20:
                if support_fraction < 0.97: # classify ambiguous when below this fraction of support
                    stability_classification = 'ambiguous'
                else:
                    stability_classification = 'stable'
            #/
        #/
        # else if no coverage
        else:
            stability_classification = 'ND'
            
            if basecounts['N'] > 0 and basecounts['N'] / basecounts['coverage'] > 0.5:
                stability_classification = 'ND'
                print('High fraction of N!')
                sys.exit()
        #/
        # bugcheck if no classification was made
        if stability_classification == None:
            print('FATAL: No classification for stability made. This is not expected, some classification should always be made')
            sys.exit()
        #/
        ##/
        # Save: For output table file with stability classification, save row
        tmp_save = {'stability_classification':stability_classification,'coverage':coverage,'ATGC_coverage':ATGC_coverage,'support_fraction':support_fraction,
                    'basecounts':basecounts}
        if not pos in msa_coords_basecounts_classified:             msa_coords_basecounts_classified[pos] = {}
        msa_coords_basecounts_classified[pos][sample] = tmp_save
        #/
###/

### Determine which positions to output as new MSA
poses_to_output_msa = set() # positions with enough "stable" SNP-datasets
poses_not_to_output_msa = set() # positions that were not determined as stable
for pos,samples_classification_data in msa_coords_basecounts_classified.items():
    # get stable/ambiguous samples (only at variant positions, i.e. consensus_is_ref should be 0 in basecount row)
    pos_stables = set()
    pos_ambiguous = set()
    pos_NDs = set()
    consensus_is_ref = set()
    
    for sample,classification_data in samples_classification_data.items():
        # check if this is a variant (if not then skip)
        basecounts = classification_data['basecounts']
        if basecounts['consensus_is_ref'] == 1:
            consensus_is_ref.add(sample)
            continue
        #/
        # get stability classification
        classification = classification_data['stability_classification']
        if classification == 'stable':
            pos_stables.add(sample)
        elif classification == 'ambiguous':
            pos_ambiguous.add(sample)
        elif classification == 'ND':
            pos_NDs.add(sample)
        #/
    #/
    
    # if there was >X% stable datasets, then save this position as "stable"
    frac_ambig = None
    if len(pos_stables) > 0: # cant determine if there are no stables at position
        frac_ambig = len(pos_ambiguous) / len(pos_stables)
        
    if frac_ambig != None and frac_ambig < ambig_output_cutoff:
        poses_to_output_msa.add(pos)
    # else, save as ambiguous or ND
    else:
        poses_not_to_output_msa.add(pos)
    #/
    
    
print(f'Positions with stables only: {len(poses_to_output_msa)}')
print(f'Positions with any ambiguous/not-determined: {len(poses_not_to_output_msa)}')
###/

### Make new MSA with stable positions only
msa_seqs_stable = {}
msa_seqs_stable_coords = [] # keep track of msa-coords <-> reference coords mapping
for pos,samples_basecounts in msa_coords_basecounts.items():
    # skip this position if not marked as stable
    if not pos in poses_to_output_msa:
        continue
    #/
    # get position of current pos in original MSA
    msa_nucl_idx = ref_coords_to_import[pos]
    #/
    # add nucleotide position for all samples
    for sample,seq in msa_seqs.items():
        sample_nucl = seq[msa_nucl_idx]
        
        # save nucl to new msa
        if not sample in msa_seqs_stable:           msa_seqs_stable[sample] = ''
        msa_seqs_stable[sample] += sample_nucl
        #/
        # keep track of msa<->reference coord map
        if not pos in msa_seqs_stable_coords:
            msa_seqs_stable_coords.append(pos)
        #/
    #/
###/

### Output
if 1:
    ## Init output dir
    if not os.path.exists(output_dir):
        print('init outdir')
        os.makedirs(output_dir)
    ##/
    
    ## Output stable MSA
    print('Write stable msa')
    with open(output_dir+'/'+'msa_stable.fasta','w') as nf:
        for sample,seq in msa_seqs_stable.items():
            seq = ''.join(seq)
            nf.write('>'+sample+'\n'+seq+'\n')
    ##/
    ## Output stable MSA coords
    print('Write stable msa coord <-> reference coord mapping file')
    with open(output_dir+'/'+'msa_stable_coords.list','w') as nf:
        nf.write(msa_slim_coords_header+'\n')
        nf.write(','.join(map(str,msa_seqs_stable_coords))+'\n')
    ##/
    ## Output position classifications
    print('Write position classification per sample')
    with open(output_dir+'/'+'position_classifications.tsv','w') as nf:
        # write header
        header = ['sample']+basecounts_header+['stability','ATGC_coverage','support_fraction','output_to_msa']
        nf.write('\t'.join(map(str,header))+'\n')
        #/
        # write rows
        for pos,samples_classification_data in msa_coords_basecounts_classified.items():
            for sample,classification_data in samples_classification_data.items():
                # compile writeArr
                writeArr = [sample]
                
                for key in basecounts_header:
                    writeArr.append(classification_data['basecounts'][key])
                
                writeArr += [classification_data['stability_classification']]
                writeArr += [classification_data['ATGC_coverage']]
                writeArr += [classification_data['support_fraction']]
                
                pos_output_to_msa = False
                if pos in poses_to_output_msa:
                    pos_output_to_msa = True
                
                writeArr += [pos_output_to_msa]
                #/
                # write
                nf.write('\t'.join(map(str,writeArr))+'\n')
                #/
        #/
    ##/
    ## Output stable coords
    print('Write stable coords list')
    with open(output_dir+'/'+'coords_output_to_msa.list','w') as nf:
        nf.write('\n'.join(map(str,sorted(list(poses_to_output_msa),key=lambda x: int(x)))))
    ##/
    ## Output ambiguous/ND coords
    print('Write ambig/ND coords list')
    with open(output_dir+'/'+'coords_not_output_to_msa.list','w') as nf:
        nf.write('\n'.join(map(str,sorted(list(poses_not_to_output_msa),key=lambda x: int(x)))))
    ##/
###/