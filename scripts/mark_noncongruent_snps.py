#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from Bio import Phylo, SeqIO
from matplotlib import pyplot as plt
import argparse

software_description = """
Produce trees with marked non-congruent SNP datasets from SNP summary file and raw tree
"""


### Parse input arguments
## setup
argparser = argparse.ArgumentParser(description=software_description)

argparser.add_argument('--tree',required=True,help='Path to tree-file (computer tree with original branch names, e.g. "CANSNPs_stable_notCrosschecked_fit_ambiguous/tree_nodifyed_raw.nwk") ')
argparser.add_argument('--msa_slim',required=True,help='Path to MSA (e.g. msa_slim.fasta)')
argparser.add_argument('--msa_slim_coords',required=True,help='Path to MSA-to-ref-coord mapping file (e.g. msa_slim_coords.list)')
argparser.add_argument('--snp_summary',required=True,help='Path to SNP summary file, will parse noncongruent SNPs (e.g. <wd>/summary_stablePipe.tsv)')
argparser.add_argument('-o','--output',required=True,help='Path to output directory')
##/
## parse input
args = argparser.parse_args()

tree_file_path = args.tree
msa_file_path = args.msa_slim
msa_coords_file_path = args.msa_slim_coords
SNP_file_path = args.snp_summary
output_dir = args.output
##/
###/


####### SCRIPT START

### import MSA and MSA-to-ref-coord map
## MSA seq
msa_seqs = {}
for record in SeqIO.parse(msa_file_path, "fasta"):
    msa_seqs[record.id] = record.seq
##/
## MSA coords
refCoord_msaIdx = {} # ref_coord -> MSA idx
msa_slim_coords_header = None
with open(msa_coords_file_path,'r') as f:
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
                refCoord_msaIdx[coord] = msa_idx
        #/
        # bugcheck if more lines
        if enum >= 2:
            print(f'WARNING: Found additional lines in {msa_coords_file_path}')
            print('This file is onyl expected to have x2 lines')
        #/
print(f'Parsed N={len(refCoord_msaIdx)} reference sequence positions from MSA file')
##/
###/

### Import tree and compute all nodes per dataset
# Import tree
tree = Phylo.read(tree_file_path, "newick")
#/
# set root
tree.root_with_outgroup('reference_sequence') # set reference-sequence as outgroup
#/
# IDE-print
if 0 and 'IDE-print tree in terminal':
    print('Imported tree:')
    Phylo.draw_ascii(tree)
#/
# get all leaf-nodes in tree
datasets_in_tree = set()
for dataset_node in tree.get_terminals():
    dataset = dataset_node.name
    if dataset == 'reference_sequence': continue # skip ref
    datasets_in_tree.add(dataset)
#/
# For each branch-node in tree, compute all childs
datasets_nodes = {} # dataset -> nodes
for dataset in datasets_in_tree:
    dataset_ancestor_path = tree.get_path(dataset)
    for node in dataset_ancestor_path:
        if node.name == dataset: continue #skip if node is self (this entry should appears to always be last but I implement this in a general way regardless)
        if not node.name in datasets_nodes:             datasets_nodes[node.name] = set()
        datasets_nodes[node.name].add(dataset)
#/
# Remove branch-labels from graph (de-clutter the plot)
for branch_node in tree.get_nonterminals():
    branch_node.name = ''
#/
###/

### Import SNP file, parse noncongruent SNPs
print(f'Importing non-congruent SNPs from file: {SNP_file_path}')
header = None
pos_nucls_datasets = {} # ref_pos -> nucl -> dataset
pos_ref_nucl = {} # ref_pos -> nucl
pos_stabilities = {} # ref_pos -> stability classification
with open(SNP_file_path,'r') as f:
    for enum,line in enumerate(f):
        # parse line
        line = line.strip('\n')
        line = line.split('\t')
        #/
        # skip header (save and move on)
        if enum == 0:
            header = line
            continue
        #/
        # check if position was incongruent
        congruency = line[8]
        
        if not congruency == 'incongruent': continue
        #/
        # parse data
        ref_pos = line[0]
        ref_nucl = line[1]
        stability = line[9]
        #/
        
        # get dataset nucleotide counts from MSA, given ref pos
        msa_idx = refCoord_msaIdx[ref_pos]
        nucls_datasets = {}
        for nucl in ('A','T','G','C','-'):
            tmp_datasets = set()
            for dataset,seq in msa_seqs.items():
                tmp_nucl = seq[msa_idx]
                if tmp_nucl == nucl:
                    tmp_datasets.add(dataset)
            nucls_datasets[nucl] = tmp_datasets
        #/
        # save
        pos_nucls_datasets[ref_pos] = nucls_datasets
        pos_ref_nucl[ref_pos] = ref_nucl
        pos_stabilities[ref_pos] = stability
        #/
print(f'Imported {len(pos_nucls_datasets)} non-congruent SNPs')
###/

### Output
# init outdir
print(f'Begin output at {output_dir}')
if not os.path.exists(output_dir):          os.makedirs(output_dir)
#/

## Make plot for each non-congruent SNP, highlighting the datasets in tree
nucl_colors = {'A':'red','T':'orange','G':'green','C':'blue','-':'black'}
for pos,nucl_datasets in pos_nucls_datasets.items():
    # INFOPRINTER
    print(f'Now plotting for position {pos}')
    #/
    # init figure
    fig, ax = plt.subplots(figsize=(10, 6))
    #/
    # determine scale to use for plotting (text and elements)
    cur_font_size = plt.rcParams['font.size']
    plot_scaler = 1
    if len(datasets_in_tree) > 20:
        plot_scaler = 0.9
    elif len(datasets_in_tree) > 30:
        plot_scaler = 0.8
    elif len(datasets_in_tree) > 40:
        plot_scaler = 0.6
    elif len(datasets_in_tree) > 50:
        plot_scaler = 0.4
    #/
    # plot tree with biopython-Phylo
    Phylo.draw(tree, do_show=False, axes=ax)
    #/
    # get position of each dataset text label
    datasets_textlabel_pos = {} # dataset -> pos
    for child in ax.get_children():
        if isinstance(child, plt.Text):
            label_text = child.get_text().lstrip(' ') # strip bullshit space infront of every text label
            # skip label if it is a branch-node or reference_sequence or is empty
            if label_text.startswith('Inner') == True or label_text == 'reference_sequence' or label_text == '': continue
            #/
            # get pos and save
            position = child.get_position()
            datasets_textlabel_pos[label_text] = position
            #/
    #/
    # get all x_poses and sort (highest X-value last, for placement of dots)
    x_vals_sorted = sorted([x for x,y in datasets_textlabel_pos.values()])
    #/
    # re-position textlabels so that they (1) align and (2) position behind dot and (3) change font size
    for child in ax.get_children():
        if isinstance(child, plt.Text):
            old_pos = child.get_position()
            child.set_position([x_vals_sorted[-1]+x_vals_sorted[-1]*0.04,old_pos[1]])
            child.set_fontsize(cur_font_size*(plot_scaler**8))
    #/
    # Plot dots for SNP datasets (A=Red,T=Orange,G=Green,C=Blue,Missing=Black)
    for dataset,text_coords in datasets_textlabel_pos.items():
        # plot dot based on nucleotide
        for nucl,datasets in nucl_datasets.items():
            # dont plot ref-nucls
            if nucl == pos_ref_nucl[pos]: continue
            #/
            # dont plot dataset if it didnt have current nucl
            if not dataset in datasets: continue
            #/
            # plot dot at dataset label
            ax.scatter(x_vals_sorted[-1]+x_vals_sorted[-1]*0.015,text_coords[1],color=nucl_colors[nucl],s=100*(plot_scaler**4))
            #/
        #/
    #/
    # Plot legend
    legend_dummys = []
    for nucl,color in nucl_colors.items():
        tmp_ = ax.scatter(0,0,color=nucl_colors[nucl],s=100*(plot_scaler**4),label=nucl)
        legend_dummys.append(tmp_)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=5,fontsize=cur_font_size*(plot_scaler**4),markerscale=plot_scaler)
    for legend_dummy in legend_dummys:
        legend_dummy.set_visible(False)
    #/
    # Cleanup plot
    ax.spines['right'].set_visible(False) # remove right plot border
    ax.spines['top'].set_visible(False) # remove top plot border
    ax.set_xlabel(None) # remove ylabel
    ax.set_ylabel(None) # remove ylabel
    ax.yaxis.set_ticks([]) # remove yticks
    ax.set_title(rf'$\bf{{snp:}}$ {pos}, $\bf{{ref:}}$ {pos_ref_nucl[pos]}, $\bf{{stability:}}$ {pos_stabilities[pos]}',fontsize=cur_font_size*plot_scaler)
    ax.title.set_position([0.5,1.1])
    
    plt.tight_layout()
    #plt.show() # IDE only
    #/
    # Save plot
    fig.savefig(output_dir+'/'+str(pos)+'.pdf')
    #/
    # close plot
    plt.close()
    #/
###/