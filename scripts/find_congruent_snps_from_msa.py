#!/usr/bin/env python3

import argparse
from Bio import Phylo,AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from matplotlib import pyplot as plt

# use tree functions from external file
import parse_pileups_functions
#/
# Import misc functions
from functions import *
#/

### Input parse
## Setup parser
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

grp_IO = parser.add_argument_group('Input and output')
grp_IO.add_argument('-i','--input',required=True,help='Path to MSA-file (e.g. msa_slim.fasta)')
grp_IO.add_argument('-c','--coordinates',required=True,help='Path to MSA coordinates-file (e.g. msa msa_slim_coords.list)')
grp_IO.add_argument('-g','--genome',required=True,help='Path to reference-genome')
grp_IO.add_argument('-o','--output',required=True,help='Path to output-directory')

grp_misc = parser.add_argument_group('Miscellaneous')
grp_misc.add_argument('--include_leafs',required=False,action='store_true',help='If specified, will include SNPs that exist in one sample (default: require at least two samples)')
##/

# Parse commandline
args = parser.parse_args()
#/
# Parse variables from command
msa_file_path = args.input
msa_nucl_refCoords_file_path = args.coordinates
refgenome_path = args.genome
output_dir = args.output

use_snps_with_one_sample = args.include_leafs
#/
## Validate
##/
###/

### Init outdir and write paramters used
mkdir(output_dir)
with open(output_dir+'/'+'args.tsv','w') as nf:
    for arg_name, arg_value in vars(args).items():
        nf.write(arg_name+'\t'+str(arg_value)+'\n')
###/

### Import MSA and make neighbour joining
print('Importing MSA')
# Load MSA
msa_aln = AlignIO.read(msa_file_path, 'fasta') # WARNING: I will modify the MSA (include reference-sequence and add preceding dummy-nucleotide. DO NOT use this downstream)
#/
## Add reference-sequence to MSA to use reference as root in the tree
# Import reference fasta
print('Importing reference-sequence')
ref_seqs = importReadSeqsFasta(refgenome_path)
ref_name = list(ref_seqs.keys())[0]
ref_seq = ref_seqs[ref_name]
#/
# import "slim" msa nucleotide positions in reference-genome
nucls_coords = None
with open(msa_nucl_refCoords_file_path,'r') as f:
    for line in f:
        if line[0] == '#': continue # skip header-line
        nucls_coords = list(map(int,line.strip('\n').split(',')))
#/
# construct reference MSA-sequence
print('Constructing MSA-sequence for reference')
ref_msa_seq = ''
for coord in nucls_coords:
    ref_msa_seq += ref_seq[coord]
#/
# add reference msa sequence to msa
if 'add_reference_to_msa':
    new_entry_record = SeqRecord(Seq(ref_msa_seq), id='reference_sequence', name='reference_sequence', description='reference_sequence')
    msa_aln.append(new_entry_record)
#/
# Add preceding dummy-nucleotide: manipulate dataset and reference sequences in the MSA to always separate reference from datasets
print('Adding preceding dummy-nucleotides to ensure the reference-sequence is always separated from the datasets')
for entry in msa_aln:
    if entry.id != 'reference_sequence':
        entry.seq = 'G' + entry.seq
    else:
        entry.seq = 'A' + entry.seq
#/
##/
# Make Neighbor-Joining tree (nj according to "Build a NJ Tree:" at https://biopython.org/docs/1.75/api/Bio.Phylo.TreeConstruction.html)
print('Constructing tree')
constructor = DistanceTreeConstructor()
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(msa_aln)
tree = constructor.nj(dm)
tree.root_with_outgroup('reference_sequence') # set reference-sequence as outgroup
#/

## Check if user wants to use SNPs with one sample. Then add a connector-node before the leaf-node for the sample so that downstream code can process it as a branch-node.
if use_snps_with_one_sample:
    print('Will allow SNPs with single-sample support. Adding "connector-nodes" before leafs in the tree')
    ## Get highest node-branch ID (I.e. "InnerX"). Will be used to name the connector-nodes
    highest_node_branch_enum = 0
    for node in tree.find_clades():
        if not node.is_terminal():
            if type(node.name) == str and node.name.find('Inner') != -1:
                node_branch_enum = int(node.name.replace('Inner',''))
                if node_branch_enum > highest_node_branch_enum:
                    highest_node_branch_enum = node_branch_enum
    ##/
    ## For each leaf-node, add a "connector-node"
    new_node_branch_enum = highest_node_branch_enum + 1 # start adding new nodes at this enumerate
    for enum,leaf_node in enumerate(tree.get_terminals()):
        # Get the parent of the leaf-node, skip the leaf-node if it does not have a parent (i.e. the case for "Reference-sequence" which is at the root)
        path_to_root = tree.get_path(leaf_node)
        if len(path_to_root) == 1: continue # if there are no parents to the current node, then skip it
        #/
        # Create the "connector-node""
        leaf_node_branch_length = leaf_node.branch_length
        connector_node = Phylo.BaseTree.Clade(name='Inner'+str(new_node_branch_enum),branch_length=0)
        new_node_branch_enum += 1 # iterate on each node add
        #/
        # Add the leaf-node as a child of the "connector-node"
        connector_node.clades.append(leaf_node)
        #/
        # replace the leaf-node with the "connector-node" at the parent.
        leaf_parent = path_to_root[-2] # get second-to-last entry. last is the leaf.
        leaf_parent.clades[leaf_parent.clades.index(leaf_node)] = connector_node
        #/
    ##/
##/

## Display the raw tree in figure
if 1 and 'plot tree?':
    fig, ax = plt.subplots(figsize=(15,15))
    Phylo.draw(tree, axes=ax, do_show=False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.title('NJ-Tree')
    plt.tight_layout()
    plt.savefig(output_dir+'/'+'tree_before.pdf')
##/

## Write tree
Phylo.write(tree, output_dir+'/'+'tree_raw.nwk', 'newick')
##/

## For each leaf-node, get ancestor path (branch-nodes), and for each ancestor, save the leaf-node
print('Computing all childs at each branch-node')
tree_nodes_childs = {}
for seq_entry in msa_aln:
    sample = seq_entry.id
    sample_ancestor_path = tree.get_path(sample)
    for node in sample_ancestor_path:
        if node.name == sample: continue #skip if node is self (this entry should appears to always be last but I implement this in a general way regardless)
        init(tree_nodes_childs,node.name,set())
        tree_nodes_childs[node.name].add(sample)
##/

## Parse sample nucleotide at each position
# import "slim" msa sequences
print('Importing MSA sequences')
msa_seqs = importReadSeqsFasta(msa_file_path)
#/
# Determine sample nucleotides per position
pos_nucl_samples = {} # pos -> nucl -> samples
for sample,seq in msa_seqs.items():
    for nucl_idx,nucl in enumerate(seq):
        ref_pos = nucls_coords[nucl_idx]
        init(pos_nucl_samples,ref_pos,{})
        init(pos_nucl_samples[ref_pos],nucl,set())
        pos_nucl_samples[ref_pos][nucl].add(sample)
#/
##/

## Find congruent SNPs using the tree: traverse SNPs and find which are congruent (i.e., which SNPs exist in all samples and only at those samples at the current node in the tree)
print('Determining congruent SNPs')
congruent_branches_snps = {} # branch_id -> SNP entries
noncongruent_snps = {} # ref_pos -> nucl -> sample
for ref_pos,nucl_samples in pos_nucl_samples.items():
    # keep track of noncongruent samples->nucls
    pos_noncongruent_snps = {} # nucl -> samples
    #/
    for nucl,samples in nucl_samples.items():
        # check if ignore this nucleotide: is not a SNP (it is same as ref)
        ref_nucl = ref_seq[ref_pos]
        if nucl == ref_nucl:
            continue
        #/
        # check if ignore this nucleotide: is "missing"
        if nucl == '-': continue
        #/
        # check if congruent (only if not a missing data)
        pos_had_congruent = False
        if nucl != '-':
            congruent_branch_node = parse_pileups_functions.check_sample_tree_conflict(samples,tree_nodes_childs,return_congruent_node=True)
            if congruent_branch_node:
                init(congruent_branches_snps,congruent_branch_node,[])
                congruent_branches_snps[congruent_branch_node].append([ref_pos,nucl,samples])
                pos_had_congruent = True
        #/
        # check if it was a noncongruent snp, then add nucl->samples
        if not pos_had_congruent:
            pos_noncongruent_snps[nucl] = samples
        #/
    # check if this position had noncongruent snp, then save to outer
    if len(pos_noncongruent_snps) > 0:
        # get total number of non-congruent samples (can be i.e. x1 for nucl1 and x1 for nucl2, it is not always x2+ samples for a single nucl)
        noncongruent_samples = set()
        for nucl,samples in pos_noncongruent_snps.items():
            noncongruent_samples.update(samples)
        #/
        # save to outer (do not consider single-samples as "noncongruent". A noncongruent positions must always have at least x2 samples for there to be a conflict. A N=1 case should only happen when --include_leafs is not applied)
        if len(noncongruent_samples) >= 2:
            noncongruent_snps[ref_pos] = nucl_samples
        #/
    #/
##/

### Change node labels (empty if no congruent SNP and otherwise to the congruent SNPs)
# Change branch-node name
for node in tree.find_clades():
    if not node.is_terminal():
        if node.name in congruent_branches_snps:
            node.name = str(len(congruent_branches_snps[node.name]))
        else:
            node.name = ''
#/
## Display the after-tree in figure
if 1 and 'plot tree?':
    fig, ax = plt.subplots(figsize=(15,15))
    Phylo.draw(tree, axes=ax, do_show=False)
    
    ## reposition branch-node labels
    x_min, x_max = ax.get_xlim()
    x_diff = x_max-x_min
    y_min, y_max = ax.get_ylim()
    y_diff = y_max-y_min
    for text_obj in ax.texts:
        if text_obj.get_text().strip().isdigit(): # only edit if the text is at a LABELLED branch-node (expect a digit)
            x, y = text_obj.get_position()
            text_obj.set_position([x-x_diff*0.015,y+y_diff*0.0075]) # move slightly
            text_obj.set_color('red')
            
    ##/
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.title('NJ-Tree')
    plt.tight_layout()
    plt.savefig(output_dir+'/'+'tree_after.pdf')
##/
###/

## Write tree (Branch-names set only at congruent branch-nodes [=number of SNPs at node])
Phylo.write(tree, output_dir+'/'+'tree_clean.nwk', 'newick')
##/

### Write congruent SNPs
with open(output_dir+'/'+'SNPs.tsv','w') as nf:
    # write header
    header = ['internal_branch_id','ref','ref_pos','ref_nucl','nucl','samples']
    nf.write('\t'.join(map(str,header))+'\n')
    #/
    for branch_node,entries in sorted(congruent_branches_snps.items(),key=lambda x: int(x[0].replace('Inner','')), reverse=True): # Try to traverse congruent branch-nodes sorted by their digit. Expected format: "Inner<num>"
        for ref_pos,nucl,samples in entries:
            ref_nucl = ref_seq[ref_pos]
            
            writeArr = [branch_node,refgenome_path,ref_pos,ref_nucl,nucl,','.join(sorted(samples))]
            nf.write('\t'.join(map(str,writeArr))+'\n')
###/

### Write noncongruent SNPs
with open(output_dir+'/'+'SNPs_nonCongruent.tsv','w') as nf:
    # write header
    header = ['ref_pos','ref_nucl','A_num_samples','T_num_samples','G_num_samples','C_num_samples','missing_num_samples','A_samples','T_samples','G_samples','C_samples','missing_samples']
    nf.write('\t'.join(map(str,header))+'\n')
    #/
    # write rows
    for ref_pos,nucls_samples in noncongruent_snps.items():
        # get ref nucl
        ref_nucl = ref_seq[ref_pos]
        #/
        # get samples for each nucleotide   
        A_samples = [] if not 'A' in nucls_samples else list(nucls_samples['A'])
        T_samples = [] if not 'T' in nucls_samples else list(nucls_samples['T'])
        G_samples = [] if not 'G' in nucls_samples else list(nucls_samples['G'])
        C_samples = [] if not 'C' in nucls_samples else list(nucls_samples['C'])
        missing_samples = [] if not '-' in nucls_samples else list(nucls_samples['-'])
        #/
        # determine nucl counts
        A_num_samples = len(A_samples)
        T_num_samples = len(T_samples)
        G_num_samples = len(G_samples)
        C_num_samples = len(C_samples)
        missing_num_samples = len(missing_samples)
        #/
        # write row
        writeArr = [ref_pos,ref_nucl,A_num_samples,T_num_samples,G_num_samples,C_num_samples,missing_num_samples,','.join(A_samples),','.join(T_samples),','.join(G_samples),','.join(C_samples),','.join(missing_samples)]
        nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
    #/
###/