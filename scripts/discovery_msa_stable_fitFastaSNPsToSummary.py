#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import copy
from Bio import SeqIO

software_description = """
Given a summary-file of MSA stable pipe, identify SNPs in a collection of fasta genomes pseudo-MSA (progressivemauve+x2fa pipe) and fit them into each subdir tree.

An input table is required of minimum:
    (1) fasta-name
    (2) subclade cansnp identifier
"""


### Parse input arguments
## setup
argparser = argparse.ArgumentParser(description=software_description)

argparser.add_argument('--msa',required=True,help='Path to fasta MSA')
argparser.add_argument('--table',required=True,help='Path to dataset <-> subclade table. Expects dataset name in column0 and subclade in column1')
argparser.add_argument('--summary',required=True,help='Path to summary of all subclades/groups')
argparser.add_argument('--msa_ref_name',required=False,default='FSC200',help='Name of reference in MSA')
argparser.add_argument('-o','--output',required=True,help='Path to output file')
##/
## parse input
args = argparser.parse_args()

msa_path = args.msa
table_path = args.table
summary_path = args.summary
msa_fasta_reference_name = args.msa_ref_name
output_file = args.output
##/
###/


### Import fasta-MSA
fasta_msa = {} # seqheader -> MSA-seq
for record in SeqIO.parse(msa_path, "fasta"):
    dataset = record.id
    fasta_msa[dataset] = record.seq

print(f'Import N={len(fasta_msa)} sequences including reference')
###/

### Import fasta<->subclade
fasta_subclade = {} # dataset -> subclade
fasta_subclade_header = None
with open(table_path,'r') as f:
    for enum,line in enumerate(f):
        # parse line
        line = line.strip('\n')
        line = line.split('\t')
        #/
        # skip if empty
        if set(line) == {''}: continue
        #/
        # parse header
        if enum == 0:
            fasta_subclade_header = line
            continue
        #/
        # parse row
        dataset = line[0]
        subclade = line[1]
        #/
        # check if this was already saved (not expected)
        if dataset in fasta_subclade:
            print(f'WARNING: dataset "{dataset}" laready had a subclade assigned. This indicates that there are duplicate rows in the input')
            print('Terminating!')
            sys.exit()
        #/
        # save
        fasta_subclade[dataset] = subclade
        #/
print(f'Saved subclade information for N={len(fasta_subclade)} datasets')

## restructure: subclade -> datasets
subclade_datasets = {}
for dataset,subclade in fasta_subclade.items():
    if not subclade in subclade_datasets:       subclade_datasets[subclade] = set()
    subclade_datasets[subclade].add(dataset)
print('Parsed the following number of datasets per subclade:')
for k,v in sorted(subclade_datasets.items(),key=lambda x: len(x[1]),reverse=True):
    print(f'{k:10}{len(v)}')
##/
###/

### Import summary
summary_rows = []
summary_header = None
with open(summary_path,'r') as f:
    for enum,line in enumerate(f):
        # parse line
        line = line.strip('\n')
        line = line.split('\t')
        #/
        # parse header
        if enum == 0:
            summary_header = line
            continue
        #/
        # convert row to dict
        tmp_dict = {}
        for idx,val in enumerate(line):
            col = summary_header[idx]
            tmp_dict[col] = val
        #/
        # save row
        summary_rows.append(tmp_dict)
        #/
print(f'Imported N={len(summary_rows)} rows from summary file')
###/

### Check fasta variants at each summary position
ref_coord_key = 'ref_coord'
ref_nucl_key = 'ref_nucl'
alt_nucl_key = 'nucl'
subclade_key = 'subclade'
for row_dict in summary_rows:
    # get summary values
    ref_coord = row_dict[ref_coord_key]
    ref_nucl = row_dict[ref_nucl_key]
    alt_nucl = row_dict[alt_nucl_key]
    subclade = row_dict[subclade_key]
    #/
    # FIX: if this subclade did not have any imported fasta datasets, init it as empty in "subclade_datasets"
    if not subclade in subclade_datasets:           subclade_datasets[subclade] = set()
    #/
    #/
    # BUGCHECK: Ensure that fasta MSA has the same nucleotide as reported fasta nucleotide in summary file
    msa_ref_nucl = fasta_msa[msa_fasta_reference_name][int(ref_coord)]
    if not ref_nucl == msa_ref_nucl:
        print('FATAL: Reference nucleotide in summary does not match the MSA!')
        print(f'Summary ref: {ref_nucl}')
        print(f'MSA ref: {msa_ref_nucl}')
        print('Terminating!')
        sys.exit()
    #/
    # Get nucleotide counts for datasets at this pos
    nucls_datasets = {} # nucl -> dataset
    for tmp_nucl in ('A','T','G','C','-','RYSWKHBVDN'): # "RYSWKHBVDN" is ambiguous nucleotide for any A/G (purine) or C/T (pyrimidine) G/C A/T G/T, etc, or N (any nucleotide)
        nucls_datasets[tmp_nucl] = set() # init
    
    for dataset,seq in fasta_msa.items():
        # skip if ref
        if dataset == msa_fasta_reference_name: continue
        #/
        ## parse nucl
        msa_nucl = seq[int(ref_coord)]
        
        # check if parsed msa nucl was undefined
        if not msa_nucl in nucls_datasets:
            print(f'WARNING: Found nucleotide {msa_nucl} in MSA but it was not initiated in the "nucls_datsets" variable. Will reassign it as "RYSWKHBVDN"')
            msa_nucl = 'RYSWKHBVDN'
        #/
        
        nucls_datasets[msa_nucl].add(dataset)
        ##/
    #/
    # re-compute above but without ref-nucl
    nucls_datasets_noRef = copy.deepcopy(nucls_datasets)
    del nucls_datasets_noRef[msa_ref_nucl]
    #/
    # count how many datasets have different than ref
    datasets_diff_from_ref = set()
    for msa_nucl,datasets in nucls_datasets.items():
        if not msa_nucl == msa_ref_nucl and not msa_nucl == '-':
            datasets_diff_from_ref.update(datasets)
    #/
    # count how many datasets have the current row alt_nucl or other non-ref nucl
    datasets_has_alt_nucl = nucls_datasets[alt_nucl]
    datasets_has_different_alt_nucl = datasets_diff_from_ref.difference(datasets_has_alt_nucl)
    #/
    # Separate counts for datasets inside/outside subclade that have the current row alt_nucl or other non-ref nucl
    subclade_datasets_has_alt_nucl = datasets_has_alt_nucl.intersection(subclade_datasets[subclade])
    subclade_datasets_has_different_nucl = datasets_has_different_alt_nucl.intersection(subclade_datasets[subclade])
    
    oSubclade_datasets = set()
    for tmp_subclade,tmp_datasets in subclade_datasets.items():
        if not tmp_subclade == subclade:
            oSubclade_datasets.update(tmp_datasets)
    oSubclade_datasets_has_alt_nucl = datasets_has_alt_nucl.intersection(oSubclade_datasets)
    oSubclade_datasets_has_different_nucl = datasets_has_different_alt_nucl.intersection(oSubclade_datasets)
    #/
    # Get missing value datasets by subclade/oSubclade
    subclade_datasets_missing_val = nucls_datasets['-'].intersection(subclade_datasets[subclade])
    oSubclade_datasets_missing_val = nucls_datasets['-'].intersection(oSubclade_datasets)
    #/
    # save and move on
    row_dict['subclade_datasets_has_alt_nucl'] = subclade_datasets_has_alt_nucl
    row_dict['subclade_datasets_has_different_alt_nucl'] = datasets_has_different_alt_nucl
    row_dict['subclade_datasets_missing_val'] = subclade_datasets_missing_val
    row_dict['oSubclade_datasets_has_alt_nucl'] = oSubclade_datasets_has_alt_nucl
    row_dict['oSubclade_datasets_has_different_nucl'] = oSubclade_datasets_has_different_nucl
    row_dict['oSubclade_datasets_missing_val'] = oSubclade_datasets_missing_val
    #/
###/

### Dump output
with open(output_file,'w') as nf:
    # write header
    header_out = summary_header + ['num_thisAlt_thisSubclade','num_otherAlt_thisSubclade','num_missing_thisSubclade','num_thisAlt_otherSubclade','num_otherAlt_otherSubclade','num_missing_otherSubclade','datasets_thisAlt_thisSubclade','datasets_thisAlt_otherSubclade']
    nf.write('\t'.join(map(str,header_out))+'\n')
    #/
    # write rows
    for row_dict in summary_rows:
        # init
        writeArr = []
        #/
        # fill imported summary values
        for key in summary_header:
            writeArr.append(row_dict[key])
        #/
        # add values from this script
        writeArr.append( len(row_dict['subclade_datasets_has_alt_nucl']) )
        writeArr.append( len(row_dict['subclade_datasets_has_different_alt_nucl']) )
        writeArr.append( len(row_dict['subclade_datasets_missing_val']) )
        writeArr.append( len(row_dict['oSubclade_datasets_has_alt_nucl']) )
        writeArr.append( len(row_dict['oSubclade_datasets_has_different_nucl']) )
        writeArr.append( len(row_dict['oSubclade_datasets_missing_val']) )
        
        writeArr.append( ','.join(row_dict['subclade_datasets_has_alt_nucl']) )
        writeArr.append( ','.join(row_dict['oSubclade_datasets_has_alt_nucl']) )
        #/
        # write
        nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
    #/
###/