#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse

software_description = """
Generate total summary of all subclades and include columns to locate "hotspots":
    1. Import stable-pipe summaries of all subclades
    2. Sort by coordinate
    3. Chain positions within X base pairs
    4. Require at least X number of SNPs and X number of datasets to call a cluster (e.g., require at least x2 SNPs [in any amount of datasets] and at least x2 datasets [in any amount of SNPs])
    
Make new columns:
    A. is_hotspot (BOOL; see point 3 above)
    B. num_datasets
    c. num_SNPs
"""


### Parse input arguments
## setup
argparser = argparse.ArgumentParser(description=software_description)

argparser.add_argument('-i','--input',required=True,help='Path to DICOVERY maindir. All subclades should be subfolders in this dir. Expects file "summary_stablePipe.tsv" at <discovery_maindir>/<subclade_dir>/summary_stablePipe.tsv>')
argparser.add_argument('-o','--output',required=True,help='Path to output file')

argparser.add_argument('--clust_dist',type=int,required=False,default=500,help='Cluster positions within this distance (default:500)')
##/
## parse input
args = argparser.parse_args()

discovery_maindir_path = args.input
output_file = args.output

clust_dist = args.clust_dist
##/
###/

##### HARDCODED
ref_coord_arrIdx = 0
snp_nucl_arrIdx = 6
snp_datasets_arrIdx = 7
#####/

### Import all subclade stablePipe summaries
summary_header = None
summary_rows = []
summary_rowHash_to_subdir = {} # hash of row -> subdir
for subdir in os.listdir(discovery_maindir_path):
    # skip if this was not a folder
    if os.path.isdir(discovery_maindir_path+'/'+subdir):
        expected_file_path = discovery_maindir_path+'/'+subdir+'/'+'summary_stablePipe.tsv'
        
        # check if expected path does not exist
        if not os.path.exists(expected_file_path):
            print(f'FATAL: Expected file at this path: {expected_file_path}')
            print('Terminating')
            sys.exit()
        #/
        # parse summary-file
        with open(expected_file_path,'r') as f:
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
                # save row
                summary_rows.append(line)
                #/
                
                # save subdir tracker
                row_hash = hash(str(line))
                summary_rowHash_to_subdir[row_hash] = subdir
                #/
        #/
        # INFO
        print(f'Done parsing {subdir}. Cumulative number of lines imported: {len(summary_rows)}')
        #/
    #/
###/

### Sort by ref-coordinate and identify hotspots
## Sort
summary_rows_sorted_wOrigEnum = []
for enum,row in enumerate(sorted(summary_rows,key=lambda x: int(x[ref_coord_arrIdx]))):
    summary_rows_sorted_wOrigEnum.append([enum,row])
##/
## chain
rows_chained = [] # chain rows by coord into this arr
print(f'Will sort by key "{summary_header[ref_coord_arrIdx]}" (should be reference coordinate)')
for rowEnum,row in summary_rows_sorted_wOrigEnum:
    # if no chain was initiated, then initiate chain
    tmp_init_new = {'rowEnums':[rowEnum],'rows':[row]} # also used to make new chain after testing distance to a previous chain
    if len(rows_chained) == 0:
        rows_chained.append(tmp_init_new)
        continue
    #/
    # check if append previous chain, else init new chain
    dist_to_prev = int(row[ref_coord_arrIdx]) - int(rows_chained[-1]['rows'][-1][ref_coord_arrIdx])
    if dist_to_prev < clust_dist:
        rows_chained[-1]['rowEnums'].append(rowEnum) # save original rowEnum from summary_rows
        rows_chained[-1]['rows'].append(row)
    else:
        rows_chained.append(tmp_init_new)
        continue
    #/
##/
## INFO
INFO_num_clusts = len(rows_chained)
INFO_num_clusts_nonSingletons = 0
INFO_nonSingletons_rows = 0
for clust in rows_chained:
    if len(clust['rows']) > 1:
        INFO_num_clusts_nonSingletons += 1
        INFO_nonSingletons_rows += len(clust['rows'])
print(f'From N={len(summary_rows)} rows, formed N={INFO_num_clusts} clusters of which N={INFO_num_clusts_nonSingletons} clusts are not singletons with a total of N={INFO_nonSingletons_rows} rows')
##/
###/

### Write output
with open(output_file,'w') as nf:
    # write header
    out_header = ['clustEnum','num_datasets','num_snps','>=2datasets_>=2snps','subclade'] + summary_header
    nf.write('\t'.join(out_header)+'\n')
    #/
    # write rows
    for clustEnum,clust in enumerate(rows_chained):
        for idx,row in enumerate(clust['rows']):
            # get rowEnum
            rowEnum = clust['rowEnums'][idx]
            #/
            # init
            writeArr = []
            #/
            # clust-info
            writeArr.append(f'clust{clustEnum}')
            
            clust_datasets = set()
            clust_snps = set()
            for tmp_row in clust['rows']:
                datasets = tmp_row[snp_datasets_arrIdx].split(',')
                snp = tmp_row[ref_coord_arrIdx]+':'+tmp_row[snp_nucl_arrIdx]
                clust_datasets.update(datasets)
                clust_snps.add(snp)
            writeArr.append(len(clust_datasets))
            writeArr.append(len(clust_snps))
            #/
            # check if clust has X datasets and X SNPs
            clust_passes_thresholds = False
            if len(clust_datasets) >= 2 and len(clust_snps) >= 2:
                clust_passes_thresholds = True
            writeArr.append(clust_passes_thresholds)
            #/
            # add subclade
            row_hash = hash(str(row))
            writeArr.append(summary_rowHash_to_subdir[row_hash])
            #/
            
            # append rest of row
            writeArr += row
            #/
            # write
            nf.write('\t'.join(map(str,writeArr))+'\n')
            #/
    #/
###/