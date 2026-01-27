#!/usr/bin/env python3

import pysam
import sys
import argparse
import os
import gzip

# Import misc functions
from functions import *
#/

### Input parse
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-i','--input',required=True,help='Path to input BAM-file')
parser.add_argument('-g','--genome',required=True,help='Path to reference-genome')
parser.add_argument('-o','--output',required=True,help='Path to output table of position base counts')

args = parser.parse_args()

bam_file_path = args.input
refgenome_path = args.genome
output_file_path = args.output
###/

### Parse refgenome
ref_seqs = importReadSeqsFasta(refgenome_path)
ref_name = list(ref_seqs.keys())[0]
ref_seq = ref_seqs[ref_name]
###/

# Parse ATCG counts in all positions
print('Begin: '+bam_file_path)
bam_fo = pysam.AlignmentFile(bam_file_path, 'rb')
position_base_counts = {}

for enum,pileup_column in enumerate(bam_fo.pileup()):
    
    # In the future, might pass some arguments to bam_fo.pileup():
    # min_base_quality | (int) – Minimum base quality. Bases below the minimum quality will not be output. The default is 13.
    # stepper='all' | "skip reads in which any of the following flags are set: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP"
    # min_mapping_quality | (int) – only use reads above a minimum mapping quality. The default is 0.
    # SEE MORE AT: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile.pileup
    
    # info printer
    if enum%10**6 == 0:
        print('Processed '+str(enum/10**6)+'M entries')
    #/
    
    position = pileup_column.pos
    depth = pileup_column.nsegments
    
    base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0,'N':0}
    other_flags = {'del':0,'refskip':0}
    for pileup_read in pileup_column.pileups:
        # check flags
        if pileup_read.is_del:
            other_flags['del'] += 1
        elif pileup_read.is_refskip:
            other_flags['refskip'] += 1
        #/
        if not pileup_read.is_del and not pileup_read.is_refskip:
            base = pileup_read.alignment.query_sequence[pileup_read.query_position]
            if base in base_counts:
                base_counts[base] += 1
            else:
                sys.exit('Found an unexpected base! '+str(base)+'\nTerminating!')
    
    position_base_counts[position] = base_counts
bam_fo.close()
#/

# Fill blanks
for pos,nucl in enumerate(ref_seq):
    if not pos in position_base_counts:
        position_base_counts[pos] = {'A': 0, 'T': 0, 'G': 0, 'C': 0,'N':0}
#/

# Output base counts to file
dirname = os.path.dirname(output_file_path)
if dirname != '' and not os.path.exists(dirname):             os.makedirs(dirname)
with gzip.open(output_file_path, 'wt') as output_file:
    header = ['pos','A','T','G','C','N','-','coverage','consensus','is_ref','consensus_is_ref']
    output_file.write('\t'.join(map(str,header))+'\n')
    for position, counts in sorted(position_base_counts.items()):
        # Get coverage
        pos_cov = sum(counts.values())
        #/
        # Check if all was zero
        no_counts = 0
        if pos_cov == 0:
            no_counts = 1
        #/
        # Get consensus
        consensus_nucl, consensus_count = sorted(counts.items(),key=lambda x: x[1], reverse=True)[0]
        
        if no_counts:   consensus_nucl = '-' # check if row had 0 counts, then dont make a consensus at position
        #/
        # check if nucleotide at position is same as reference (single-reported and consensus)
        is_ref = 0
        consensus_is_ref = 0
        if consensus_count == pos_cov and ref_seq[position] == consensus_nucl:          is_ref = 1 # single-reported nucleotide at position
        if  ref_seq[position] == consensus_nucl:                                        consensus_is_ref = 1 # consensus
        #/
        
        writeArr = [position,counts['A'],counts['T'],counts['G'],counts['C'],counts['N'],no_counts,pos_cov,consensus_nucl,is_ref,consensus_is_ref]
        output_file.write('\t'.join(map(str,writeArr))+'\n')

print('Finished: '+output_file_path)
#/
