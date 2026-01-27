#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import argparse
import ast
from Bio import SeqIO

software_description = """
Summary file for "stable MSA" pipeline
"""


### Parse input arguments
## setup
argparser = argparse.ArgumentParser(description=software_description)

argparser.add_argument('-i','--input',required=True,help='Path to DICOVERY clade results dir')
argparser.add_argument('-b','--basecounts',required=True,help='Path to directory tables with position base counts')
argparser.add_argument('-r','--ref','--reference','--ref_fasta',required=True,help='Path to REF fasta')
argparser.add_argument('-o','--output',required=True,help='Path to output file')
##/
## parse input
args = argparser.parse_args()

subclade_dir = args.input
base_counts_tables_path = args.basecounts
ref_fasta = args.ref
output_file = args.output
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
            if not pos in pos_basecounts:           pos_basecounts[pos] = []
            pos_basecounts[pos].append(line)
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

### INFO
print(f'Begin for: {subclade_dir}')
###/

### Get data from files
## MSA/msa_slim_coords.list
expected_path = subclade_dir+'/'+'MSA/msa_slim_coords.list'
MSA_msa_slim_coords = []
MSA_msa_slim_coords_header = None
with open(expected_path,'r') as f:
    for enum,line_raw in enumerate(f):
        # parse line
        line = line_raw.strip('\n')
        if line == '': continue
        line = line.split('\t')
        #/
        # parse header
        if enum == 0:
            MSA_msa_slim_coords_header = line
            continue
        #/
        # parse coordinates (only row in file)
        coords = line[0].split(',')
        #/
        # bugcheck
        if len(line) > 1:
            print('WARNING: multiple columns (tabseparated) found, this is not expected. Expected only a comma-separated list')
        if enum >= 2:
            print('WARING: found multiple rows. Expected only one row after header')
        #/
        # save
        MSA_msa_slim_coords = coords
        #/
##/

## MSA-slim (used to determine nucleotide of datasets)
expected_path = subclade_dir+'/'+'MSA/msa_slim.fasta'
MSA_slim_coords_nucls_datasets = {} # ref_coord -> nucleotide -> datasets
datasets_within_subclade = set()
for record in SeqIO.parse(expected_path, "fasta"):
    dataset = record.id
    for coord_idx,nucl in enumerate(record.seq):
        ref_coord = MSA_msa_slim_coords[coord_idx]
        
        if not ref_coord in MSA_slim_coords_nucls_datasets:           MSA_slim_coords_nucls_datasets[ref_coord] = {}
        if not nucl in MSA_slim_coords_nucls_datasets[ref_coord]:     MSA_slim_coords_nucls_datasets[ref_coord][nucl] = set()
        MSA_slim_coords_nucls_datasets[ref_coord][nucl].add(dataset)
        
        datasets_within_subclade.add(dataset)
##/

## MSA_slim positions annotated with reference info. See print below for help
expected_path = subclade_dir+'/'+'MSA/msa_slim_coords_refinfo.tsv'
# Check if file exists. It is produced by an external "misc" script. Tell user how to generate it, should is not be present
if not os.path.exists(expected_path):
    print(f'Could not locate file with reference-annotation (gff) for MSA_slim_coords. It was expected at {expected_path}')
    print('You can generate it with "MISC-script"')
    print('discovery_output_ref_nucl_and_annotation.py --reference ref.fasta --gff ref.gff --input  <subclade_workdir>/MSA/msa_slim_coords.list --split_input_by_comma --output <subclade_workdir>/MSA/MSA_slim_coords_refinfo.tsv')
    sys.exit()
#/

MSA_msa_slim_coords_refAnnotation = []
MSA_msa_slim_coords_refAnnotation_header = []
with open(expected_path,'r') as f:
    for enum,line_raw in enumerate(f):
        # parse line
        line = line_raw.strip('\n')
        if line == '': continue
        line = line.split('\t')
        #/
        # parse header
        if enum == 0:
            MSA_msa_slim_coords_refAnnotation_header = line
            continue
        #/
        # parse row
        MSA_msa_slim_coords_refAnnotation.append(line)
        #/
##/

## Traverse basecounts files
poses_datasets_consensus_counts = {} # pos -> nucl -> dataset_consensus_nucl
for enum,file_ in enumerate(os.listdir(base_counts_tables_path)):
    # INFOPRINT
    print(f'Now processing file {enum} of {len(os.listdir(base_counts_tables_path))-1}')
    #/
    # determine dataset name and file path
    dataset_name = file_.strip('.gz')
    file_path = base_counts_tables_path + '/' + file_
    dataset_basecounts = parse_basecounts(file_path,poses_to_get = MSA_msa_slim_coords)
    #/
    # save to outer
    for pos,lines in dataset_basecounts.items():
        # bugcheck
        if len(lines) > 1:
            print('WARNING: Multiple lines found. Only expected one')
            sys.exit()
        #/
        # parse line
        _,A,T,G,C,N,missing,coverage,consensus,is_ref,consensus_is_ref = lines[0]
        A,T,G,C,N,missing = map(int,[A,T,G,C,N,missing])
        coverage = int(coverage)
        #/
        # save to outer
        if not pos in poses_datasets_consensus_counts:                  poses_datasets_consensus_counts[pos] = {}
        if not consensus in poses_datasets_consensus_counts[pos]:       poses_datasets_consensus_counts[pos][consensus] = set()
        poses_datasets_consensus_counts[pos][consensus].add(dataset_name)
        #/
    #/
##/

## Position stability information
expected_path = subclade_dir+'/'+'MSA_stable/position_classifications.tsv'
MSA_stable_position_classifications = []
MSA_stable_position_classifications_header = []
with open(expected_path,'r') as f:
    for enum,line_raw in enumerate(f):
        # parse line
        line = line_raw.strip('\n')
        if line == '': continue
        line = line.split('\t')
        #/
        # parse header
        if enum == 0:
            MSA_stable_position_classifications_header = line
            continue
        #/
        # parse row
        MSA_stable_position_classifications.append(line)
        #/

# Determine overall position stability (i.e. if it was output to stable MSA)
positions_with_stable_classification = set()
positions_with_instable_classification = set()
for line in MSA_stable_position_classifications:
    pos = line[1]
    pos_output_to_stable_msa_BOOL = ast.literal_eval(line[15])
    if pos_output_to_stable_msa_BOOL == True:
        positions_with_stable_classification.add(pos)
    else:
        positions_with_instable_classification.add(pos)
print(f'Parsed N={len(positions_with_stable_classification)} positions with "Stable" overall classification and N={len(positions_with_instable_classification)} positions with instable overall classification')
#/
##/

## Congruent SNPs from stable tree
expected_path = subclade_dir+'/'+'CANSNP_candidates/SNPs.tsv'
CANSNPs_stable_notCrosschecked_SNPs = []
CANSNPs_stable_notCrosschecked_SNPs_header = []
with open(expected_path,'r') as f:
    for enum,line_raw in enumerate(f):
        # parse line
        line = line_raw.strip('\n')
        if line == '': continue
        line = line.split('\t')
        #/
        # parse header
        if enum == 0:
            CANSNPs_stable_notCrosschecked_SNPs_header = line
            continue
        #/
        # parse row
        CANSNPs_stable_notCrosschecked_SNPs.append(line)
        #/
##/

## Incongruent SNPs after stable tree
expected_path = subclade_dir+'/'+'CANSNP_candidates/SNPs_nonCongruent.tsv'
CANSNPs_stable_notCrosschecked_SNPs_nonCongruent = []
CANSNPs_stable_notCrosschecked_SNPs_nonCongruent_header = []
with open(expected_path,'r') as f:
    for enum,line_raw in enumerate(f):
        # parse line
        line = line_raw.strip('\n')
        if line == '': continue
        line = line.split('\t')
        #/
        # parse header
        if enum == 0:
            CANSNPs_stable_notCrosschecked_SNPs_nonCongruent_header = line
            continue
        #/
        # parse row
        CANSNPs_stable_notCrosschecked_SNPs_nonCongruent.append(line)
        #/
##/

## Congruent SNPs that added nodes to tee (fit-ambiguous SNPs to tree)
expected_path = subclade_dir+'/'+'CANSNP_candidates_fitAmbiguous/SNPs_fitted.tsv'

CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_fitted = []
CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_fitted_header = []
with open(expected_path,'r') as f:
    for enum,line_raw in enumerate(f):
        # parse line
        line = line_raw.strip('\n')
        if line == '': continue
        line = line.split('\t')
        #/
        # parse header
        if enum == 0:
            CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_fitted_header = line
            continue
        #/
        # parse row
        CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_fitted.append(line)
        #/
##/

## Instable SNPs after fitting ambiguous to tree
expected_path = subclade_dir+'/'+'CANSNP_candidates_fitAmbiguous/SNPs_instable.tsv'
CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_instable = []
CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_instable_header = []
with open(expected_path,'r') as f:
    for enum,line_raw in enumerate(f):
        # parse line
        line = line_raw.strip('\n')
        if line == '': continue
        line = line.split('\t')
        #/
        # parse header
        if enum == 0:
            CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_instable_header = line
            continue
        #/
        # parse row
        CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_instable.append(line)
        #/
##/

## Incongruent SNPs after fitting to tree
expected_path = subclade_dir+'/'+'CANSNP_candidates_fitAmbiguous/SNPs_nonCongruent.tsv'
CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_nonCongruent = []
CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_nonCongruent_header = []
with open(expected_path,'r') as f:
    for enum,line_raw in enumerate(f):
        # parse line
        line = line_raw.strip('\n')
        if line == '': continue
        line = line.split('\t')
        #/
        # parse header
        if enum == 0:
            CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_nonCongruent_header = line
            continue
        #/
        # parse row
        CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_nonCongruent.append(line)
        #/
##/
###/

### Output
## Compile output: Run through all ref-positions and MSA nucleotide variants to check which nucleotides exist at each "MSA slim position"
pos_nucl_data = {} # pos -> nucl -> data
arr_of_output_dicts = [] # contains dict-representation of rows
for ref_coord in MSA_slim_coords_nucls_datasets:
    for msa_nucl in MSA_slim_coords_nucls_datasets[ref_coord]:
        # skip this msa_nucl if it is the same as ref
        if msa_nucl == ref_seq[int(ref_coord)]: continue
        #/
        # skip this msa_nucl if it was a gap
        if msa_nucl == '-': continue
        #/
        # init save-dict
        tmp_save = {}
        tmp_save['ref_coord'] = ref_coord
        tmp_save['ref_nucl'] = ref_seq[int(ref_coord)]
        tmp_save['nucl'] = msa_nucl
        #/
        # get MSA counts per nucleotide at position
        nucls_samples_at_pos = {} # nucl -> samples
        for tmp_nucl in ('A','T','G','C','-'):
            tmp_samples = set()
            if tmp_nucl in MSA_slim_coords_nucls_datasets[ref_coord]:
                tmp_samples = MSA_slim_coords_nucls_datasets[ref_coord][tmp_nucl]
            
            nucls_samples_at_pos[tmp_nucl] = tmp_samples
            tmp_save['msa_names_'+tmp_nucl] = ','.join(tmp_samples)
            tmp_save['msa_'+tmp_nucl] = len(tmp_samples)
        #/
        # get number of samples that differ from ref
        nucls_samples_at_pos_diffFromRef = set()
        for tmp_nucl,samples in nucls_samples_at_pos.items():
            if tmp_nucl == '-': continue
            if tmp_nucl == ref_seq[int(ref_coord)]: continue
            nucls_samples_at_pos_diffFromRef.update(samples)
        tmp_save['msa_diffFromRef'] = len(nucls_samples_at_pos_diffFromRef)
        tmp_save['msa_diffFromRef_names'] = ','.join(nucls_samples_at_pos_diffFromRef)
        #/
        # get ref-pos info (find row for current coord)
        for row in MSA_msa_slim_coords_refAnnotation:
            ref_pos,ref_nucl,types,names,descriptions = row
            
            # skip row if it was not the current coord
            if not ref_pos == ref_coord: continue
            #/
            
            # save
            tmp_save['ann_types'] = types
            tmp_save['ann_names'] = names
            tmp_save['ann_descriptions'] = descriptions
            #/
        #/
        
        ## Get stable MSA SNPs that were congruent to stable tree (find row for current coord+nucl)
        for row in CANSNPs_stable_notCrosschecked_SNPs:
            internal_branch_id,ref,ref_pos,ref_nucl,nucl,samples = row
            
            # skip row if it was not the current coord and current alt-nucl
            if not (ref_pos == ref_coord and msa_nucl == nucl): continue
            #/
            
            # bugcheck
            if not ref_nucl == ref_seq[int(ref_coord)]:
                print('Reference nucleotide do not match! ID:9b98h9a8y259hg')
                sys.exit()
            #/
            
            # save
            samples = set(samples.split(','))
            
            tmp_save['branch_id'] = internal_branch_id
            tmp_save['nucl'] = nucl
            tmp_save['samples'] = samples
            tmp_save['congruency'] = 'congruent'
            #/
        ##/
        
        ## Get instable SNPs that were congruent to stable tree
        for row in CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_instable:
            internal_branch_id,ref,ref_pos,ref_nucl,nucl,samples = row
            
            # skip row if it was not the current coord and current alt-nucl
            if not (ref_pos == ref_coord and msa_nucl == nucl): continue
            #/
            
            # bugcheck
            if not ref_nucl == ref_seq[int(ref_coord)]:
                print('Reference nucleotide do not match! ID:9b98h9a8y259hg')
                sys.exit()
            #/
            
            # save
            samples = set(samples.split(','))
            
            tmp_save['branch_id'] = internal_branch_id
            tmp_save['nucl'] = nucl
            tmp_save['samples'] = samples
            tmp_save['congruency'] = 'congruent'
        ##/
        
        ## Get fit-tree stuff, i.e. new nodes made by fitted SNPs (find row for current coord+nucl)
        for row in CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_fitted:
            internal_branch_id,ref,ref_pos,ref_nucl,nucl,samples = row
            
            # skip row if it was not the current coord and current alt-nucl
            if not (ref_pos == ref_coord and msa_nucl == nucl): continue
            #/
            
            # bugcheck
            if not ref_nucl == ref_seq[int(ref_coord)]:
                print('Reference nucleotide do not match! ID:9b98h9a8y259hg')
                sys.exit()
            #/
            
            ## save
            samples = set(samples.split(','))
            
            # check if this was added to a previous branch and if so then append samples
            if tmp_save != None and 'samples' in tmp_save:
                tmp_save['samples'].update(samples)
            #/
            # else, make new save-dict
            else:
                print(f'init new for branch: {internal_branch_id}')
                tmp_save['branch_id'] = internal_branch_id
                tmp_save['nucl'] = nucl
                tmp_save['samples'] = samples
                tmp_save['congruency'] = 'congruent'
            #/
            ##/
        ##/
        
        ## Get incongruent SNPs (after trying to fit to tree; find row for current coord+nucl)
        for row in CANSNPs_stable_notCrosschecked_fit_ambiguous_SNPs_nonCongruent:
            ref_pos,ref_nucl,A_num_samples,T_num_samples,G_num_samples,C_num_samples,missing_num_samples,A_samples,T_samples,G_samples,C_samples,missing_samples = row
            
            # skip row if it was not the current coord and current alt-nucl
            if not (ref_pos == ref_coord and msa_nucl != ref_nucl): continue # Get ALT-nucleotide from MSA-slim when not supplied by pipeline output files
            #/
            # bugcheck
            if not ref_nucl == ref_seq[int(ref_coord)]:
                print('Reference nucleotide do not match! ID:9b98h9a8y259hg')
                sys.exit()
            #/
            
            ## save
            #
            if 'samples' in tmp_save:
                print(f'This position ({ref_pos}) has multiple variants and this variant was previously saved. Will move on')
                break
            #/
            tmp_save['branch_id'] = 'NA'
            tmp_save['nucl'] = msa_nucl
            tmp_save['samples'] = MSA_slim_coords_nucls_datasets[ref_coord][msa_nucl]
            tmp_save['congruency'] = 'incongruent'
            ##
        ##/
        
        ## Get stable MSA SNPs that were incongruent to stable tree (find row for current coord+nucl)
        if not 'samples' in tmp_save:
            for row in CANSNPs_stable_notCrosschecked_SNPs_nonCongruent:
                ref_pos,ref_nucl,A_num_samples,T_num_samples,G_num_samples,C_num_samples,missing_num_samples,A_samples,T_samples,G_samples,C_samples,missing_samples = row
                
                # skip row if it was not the current coord and current msa_nucl
                if not ref_pos == ref_coord: continue
                if msa_nucl == 'A' and not int(A_num_samples) > 0: continue
                if msa_nucl == 'T' and not int(T_num_samples) > 0: continue
                if msa_nucl == 'G' and not int(G_num_samples) > 0: continue
                if msa_nucl == 'C' and not int(C_num_samples) > 0: continue
                #/
                
                # parse samples
                samples = None
                if msa_nucl == 'A': samples = set(A_samples.split(','))
                if msa_nucl == 'T': samples = set(T_samples.split(','))
                if msa_nucl == 'G': samples = set(G_samples.split(','))
                if msa_nucl == 'C': samples = set(C_samples.split(','))
                if samples == None:
                    print(f'No samples found for this nucleotide, this is unexpected. Check it out. msa_nucl={msa_nucl}, ref_coord={ref_coord}')
                #/
                
                ## save
                tmp_save['branch_id'] = 'NA'
                tmp_save['nucl'] = msa_nucl
                tmp_save['samples'] = samples
                tmp_save['congruency'] = 'incongruent'
                ##
        ##/
        
        ## Get stability for this position
        tmp_save['stability'] = 'NA'
        if ref_coord in positions_with_stable_classification:
            tmp_save['stability'] = 'stable'
        elif ref_coord in positions_with_instable_classification:
            tmp_save['stability'] = 'instable'
        
        if tmp_save['stability'] == 'NA': sys.exit('No stability classified for position. Something seems wrong, expected "stable" or "instable"')
        ##/
        
        ## get other datasets nucleotide at this position (separate datasets that are within clade from those outside clade)
        nucls_datasets_withinSubclade = {} # A/T/G/C/G/-  -> datasets
        nucls_datasets_outsideSubclade = {} # --::--
        if ref_coord in poses_datasets_consensus_counts:
            for tmp_nucl in ('A','T','G','C','-'):
                # init save
                nucls_datasets_withinSubclade[tmp_nucl] = set()
                nucls_datasets_outsideSubclade[tmp_nucl] = set()
                #/
                # save
                if tmp_nucl in poses_datasets_consensus_counts[ref_coord]:
                    for dataset in poses_datasets_consensus_counts[ref_coord][tmp_nucl]:
                        # save dataset at "within subclade" or "outside subclade"
                        if dataset in datasets_within_subclade:
                            nucls_datasets_withinSubclade[tmp_nucl].add(dataset)
                        else:
                            nucls_datasets_outsideSubclade[tmp_nucl].add(dataset)
                        #/
                #/
        
        # save
        for tmp_nucl in ('A','T','G','C','-'):
            tmp_save['withinSubclade_'+tmp_nucl] = len(nucls_datasets_withinSubclade[tmp_nucl])
        for tmp_nucl in ('A','T','G','C','-'):
            tmp_save['outsideSubclade_'+tmp_nucl] = len(nucls_datasets_outsideSubclade[tmp_nucl])
        #/
        ##/
        
        ## Compute which datasets differ from ref (separate datasets that are within clade from those outside clade)
        diff_from_ref_datasets_withinSubclade = set()
        diff_from_ref_datasets_outsideSubclade = set()
        for tmp_nucl,datasets in nucls_datasets_withinSubclade.items():
            if not tmp_nucl in (ref_seq[int(ref_coord)], '-'):
                diff_from_ref_datasets_withinSubclade.update(nucls_datasets_withinSubclade[tmp_nucl])
        for tmp_nucl,datasets in nucls_datasets_outsideSubclade.items():
            if not tmp_nucl in (ref_seq[int(ref_coord)], '-'):
                diff_from_ref_datasets_outsideSubclade.update(nucls_datasets_outsideSubclade[tmp_nucl])
        
        # save
        tmp_save['withinSubclade_diffFromRef'] = len(diff_from_ref_datasets_withinSubclade)
        tmp_save['withinSubclade_diffFromRef_names'] = ','.join(list(diff_from_ref_datasets_withinSubclade))
        
        tmp_save['outsideSubclade_diffFromRef'] = len(diff_from_ref_datasets_outsideSubclade)
        tmp_save['outsideSubclade_diffFromRef_names'] = ','.join(list(diff_from_ref_datasets_outsideSubclade))
        #/
        # make column "presence_in_other" (set as TRUE if there are "outsideSubclade_diffFromRef" datasets)
        tmp_save['presence_in_other'] = False
        if len(diff_from_ref_datasets_outsideSubclade) > 0:
            tmp_save['presence_in_other'] = True
        #/
        ##/
        
        
        ###### FIX PIPELINE LIMITATIONS
        ## For cases when multiple variants exist (different from reference), downstream scripts in pipe do not classify incongruency. Check here if no branch was assigned to SNP and if not then mark it as incongruent
        num_variants_at_pos = 0
        for tmp_nucl in ('A','T','G','C'):
            if tmp_nucl == ref_seq[int(ref_coord)]: continue
            if tmp_save['withinSubclade_'+tmp_nucl] > 0:
                num_variants_at_pos += 1
        if num_variants_at_pos >= 2:
            if not 'branch_id' in tmp_save:
                print(f'Assigning "incongruent" for pos {ref_coord} and nucleotide {msa_nucl}; no branch ID stated')
                tmp_save['congruency'] = 'incongruent'
                tmp_save['samples'] = ','.join(list(nucls_datasets_withinSubclade[msa_nucl]))
        ##/
        
        ## For cases where a stable SNP was incongruent to the stable tree
        
        ######/
        
        ## Format if needed
        if 'samples' in tmp_save and type(tmp_save['samples']) != str:                tmp_save['samples'] = ','.join(list(tmp_save['samples']))
        ##/
        
        ## Save to outer
        arr_of_output_dicts.append(tmp_save)
        ##/
print(f'Fetched data for N={len(arr_of_output_dicts)} rows')
##/

## Write output
# define columns to output (in order)
columns_to_write = ['ref_coord',
                    'ref_nucl',
                    'ann_types',
                    'ann_names',
                    'ann_descriptions',
                    'branch_id',
                    'nucl',
                    'samples',
                    'congruency',
                    'stability',
                    #'withinSubclade_A',
                    #'withinSubclade_T',
                    #'withinSubclade_G',
                    #'withinSubclade_C',
                    #'withinSubclade_-',
                    'msa_A',
                    'msa_T',
                    'msa_G',
                    'msa_C',
                    'msa_-',
                    'outsideSubclade_A',
                    'outsideSubclade_T',
                    'outsideSubclade_G',
                    'outsideSubclade_C',
                    'outsideSubclade_-',
                    'presence_in_other',
                    #'withinSubclade_diffFromRef',
                    'msa_diffFromRef',
                    'outsideSubclade_diffFromRef',
                    #'withinSubclade_diffFromRef_names',
                    'msa_diffFromRef_names',
                    'outsideSubclade_diffFromRef_names',
                    ]
#/
# write output file
print('Writing output now')
with open(output_file,'w') as nf:
    # write header
    nf.write('\t'.join(map(str,columns_to_write))+'\n')
    #/
    # write rows
    for output_dict in arr_of_output_dicts:
        writeArr = []
        for key in columns_to_write:
            # init value as empty string
            val = ''
            #/
            # get value from data dictionary if it exists
            if key in output_dict:
                val = output_dict[key]
            #/
            # save
            writeArr.append(val)
            #/
        # write row
        nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
    #/
print(f'Output written to file: {output_file}')
#/
##/
###/
