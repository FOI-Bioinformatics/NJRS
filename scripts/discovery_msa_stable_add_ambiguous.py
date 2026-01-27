#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from matplotlib import pyplot as plt
from Bio import SeqIO, Phylo
import copy
import itertools



### Parse input arguments
## setup
argparser = argparse.ArgumentParser()

argparser.add_argument('--tree',required=True,help='Input tree built from stable MSA')
argparser.add_argument('--congruent_snp_file',required=True,help='SNP file with congruent SNPs assigned at each branch')
argparser.add_argument('--msa_raw',required=True,help='Tree MSA before applying filter for stable positions')
argparser.add_argument('--msa_raw_coords',required=True,help='File of reference coordinates for each position in raw MSA')
argparser.add_argument('--ambig_snp_coords',required=True,help='List of SNP positions that were not included in stable MSA')
argparser.add_argument('-r','--ref','--reference','--ref_fasta',required=True,help='Path to REF fasta')
argparser.add_argument('-o','--output','--outdir',required=True,help='Path to output dir (will put multiple files)')
##/
## parse input
args = argparser.parse_args()

tree_path = args.tree
congruent_snp_path = args.congruent_snp_file
msa_raw_path = args.msa_raw
msa_raw_coords_path = args.msa_raw_coords
ambig_snp_coords_path = args.ambig_snp_coords
ref_fasta = args.ref
output_dir = args.output
##/
###/


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

### Import MSA
## MSA seq
msa_seqs = {}
for record in SeqIO.parse(msa_raw_path, "fasta"):
    print("Sequence ID:", record.id)
    print("Description:", record.description)
    
    msa_seqs[record.id] = record.seq
##/
## MSA coords
ref_coords = {} # ref_coord -> MSA idx
msa_slim_coords_header = None
with open(msa_raw_coords_path,'r') as f:
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
                ref_coords[coord] = msa_idx
        #/
        # bugcheck if more lines
        if enum >= 2:
            print(f'WARNING: Found additional lines in {msa_raw_coords_path}')
            print('This file is only expected to have x2 lines')
        #/
print(f'Parsed N={len(ref_coords)} reference sequence positions from MSA file')
##/
###/

### Import ambiguous SNP coords
ambig_snp_coords = set()
with open(ambig_snp_coords_path,'r') as f:
    for line in f:
        line = line.strip('\n')
        if line == '': continue # skip empty line
        _ = int(line) # try to int it to make sure error is triggered if the file format is wrong. The file is expected to be a list of coordinates
        coord = line
        ambig_snp_coords.add(coord)
print(f'Parsed N={len(ambig_snp_coords)} coordinates with ambiguous SNPs')
###/

### Determine which nucleotides in input MSA that were not used in the stable tree generation
ambig_snp_MSA_nucls = {} # pos -> MSA_nucl -> datasets
for ref_coord,msa_idx in ref_coords.items():
    # skip if this ref_coord was not part of ambiguous set
    if not ref_coord in ambig_snp_coords: continue
    #/
    # save nucl->samples at this MSA position for the ambiguous nucl
    for sample,msa_seq in msa_seqs.items():
        msa_nucl = msa_seq[msa_idx]
        if not ref_coord in ambig_snp_MSA_nucls:                    ambig_snp_MSA_nucls[ref_coord] = {}
        if not msa_nucl in ambig_snp_MSA_nucls[ref_coord]:          ambig_snp_MSA_nucls[ref_coord][msa_nucl] = set()
        ambig_snp_MSA_nucls[ref_coord][msa_nucl].add(sample)
    #/
print(f'Parsed the nucleotide for each sample at each position with ambiguous coord, successful positions N={len(ambig_snp_MSA_nucls)}')
###/

### Remove values for "ref-nucl" in ambig_snp_MSA_nucls
ambig_snp_MSA_nucls_refNuclDiscarded = {}
for ref_coord,msa_nucl_samples in ambig_snp_MSA_nucls.items():
    for msa_nucl,samples in msa_nucl_samples.items():
        # only save if this nucl was not the ref nucl or a gap
        if not msa_nucl == ref_seq[int(ref_coord)] and not msa_nucl == '-':
            if not ref_coord in ambig_snp_MSA_nucls_refNuclDiscarded:                    ambig_snp_MSA_nucls_refNuclDiscarded[ref_coord] = {}
            ambig_snp_MSA_nucls_refNuclDiscarded[ref_coord][msa_nucl] = samples
        #/
###/


### Import tree and congruent SNPs. Determine branches for ambiguous SNPs
## import tree and set root
tree = Phylo.read(tree_path, "newick")
tree.root_with_outgroup('reference_sequence') # set reference-sequence as outgroup
##/
## Import stable congruent SNPs at each branch
branch_stable_and_congruent_datasets = {} # branch_node -> datasets
with open(congruent_snp_path,'r') as f:
    congruent_branch_ambig_snps_HEADER = None
    for enum,line in enumerate(f):
        # parse line
        line = line.strip('\n')
        line = line.split('\t')
        #/
        # parse header
        if enum == 0:
            congruent_branch_ambig_snps_HEADER = line
            continue
        #/
        # parse row
        tree_branch_id,ref_file_path,ref_pos,ref_nucl,nucl,samples = line
        #/
        # save
        if not tree_branch_id in branch_stable_and_congruent_datasets:         branch_stable_and_congruent_datasets[tree_branch_id] = []
        branch_stable_and_congruent_datasets[tree_branch_id].append(line)
        #/
print(f'Parsed N={len(branch_stable_and_congruent_datasets)} branches with congruent SNPs')
##/
## get all leaf-nodes in tree
datasets_in_tree = set()
for dataset_node in tree.get_terminals():
    dataset = dataset_node.name
    if dataset == 'reference_sequence': continue # skip ref
    datasets_in_tree.add(dataset)
print(f'Found N={len(datasets_in_tree)} datasets in the tree')
##/
## For each branch-node in tree, compute all childs
nodes_datasets = {} # node->datasets
for dataset in datasets_in_tree:
    dataset_ancestor_path = tree.get_path(dataset)
    for node in dataset_ancestor_path:
        if node.name == dataset: continue #skip if node is self (this entry should appears to always be last but I implement this in a general way regardless)
        if not node.name in nodes_datasets:             nodes_datasets[node.name] = set()
        nodes_datasets[node.name].add(dataset)
print(f'Assigned dataset childs to N={len(nodes_datasets)} nodes')
##/
## for each branch-node in tree, get ambiguous snps
congruent_branch_ambig_snps = {} # branch_ID -> snp_pos -> msa_nucl -> samples
ambig_snps_assignment_status = {} # keep track of which ref_pos->msa_nucl were congruent in tree
for branch_node,node_datasets in nodes_datasets.items():
    for ref_coord,msa_nucl_samples in ambig_snp_MSA_nucls_refNuclDiscarded.items():
        for msa_nucl,samples in msa_nucl_samples.items():
            # init this ref_pos->nucl in assignment_status
            if not ref_coord in ambig_snps_assignment_status:               ambig_snps_assignment_status[ref_coord] = {}
            if not msa_nucl in ambig_snps_assignment_status[ref_coord]:     ambig_snps_assignment_status[ref_coord][msa_nucl] = {'samples':samples,'status':'incongruent'} # init just once, will overwrite after testing all nodes in tree if it was found congruent
            #/
            # check if congruent
            if samples == node_datasets:
                if not branch_node in congruent_branch_ambig_snps:                   congruent_branch_ambig_snps[branch_node] = {}
                if not ref_coord in congruent_branch_ambig_snps[branch_node]:        congruent_branch_ambig_snps[branch_node][ref_coord] = {}
                congruent_branch_ambig_snps[branch_node][ref_coord][msa_nucl] = samples
                
                ambig_snps_assignment_status[ref_coord][msa_nucl]['status'] = 'congruent' # overwrite incongruent
            #/
print(f'Found N={len(congruent_branch_ambig_snps)} branch-nodes with congruent SNPs')
##/
## Sometimes the stable tree has branch-nodes that do not have a stable SNP to them (i.e. 0support). Keep track of instable SNPs that support a branch-node in the stable tree.
branch_instable_and_congruent_datasets_noStablesAtNode = set()
for branch_node in congruent_branch_ambig_snps:
    if not branch_node in branch_stable_and_congruent_datasets:
        branch_instable_and_congruent_datasets_noStablesAtNode.add(branch_node)
##/
## Parse non-congruent SNPs
congruent_snps_samples = {}
noncongruent_snps_samples = {}
INFO_num_assigned_congruent = 0
INFO_num_assigned_noncongruent = 0
for ref_coord,msa_nucl_samples in ambig_snps_assignment_status.items():
    for msa_nucl,samples_status in ambig_snps_assignment_status[ref_coord].items():
        # check if congruent
        pos_nucl_isCongruent = None
        if samples_status['status'] == 'congruent':
            pos_nucl_isCongruent = True
        if samples_status['status'] == 'incongruent':
            pos_nucl_isCongruent = False
        
        if pos_nucl_isCongruent == None: sys.exit('BUG: congruent was left as NONE, should always be True/False. Check it out')
        #/
        # save samples as congruent/incongruent
        if pos_nucl_isCongruent == True:
            if not ref_coord in congruent_snps_samples:               congruent_snps_samples[ref_coord] = {}
            if not msa_nucl in congruent_snps_samples[ref_coord]:     congruent_snps_samples[ref_coord][msa_nucl] = {}
            congruent_snps_samples[ref_coord][msa_nucl] = samples_status['samples']
            INFO_num_assigned_congruent += 1
        else:
            if not ref_coord in noncongruent_snps_samples:               noncongruent_snps_samples[ref_coord] = {}
            if not msa_nucl in noncongruent_snps_samples[ref_coord]:     noncongruent_snps_samples[ref_coord][msa_nucl] = {}
            noncongruent_snps_samples[ref_coord][msa_nucl] = samples_status['samples']
            INFO_num_assigned_noncongruent += 1
        #/
print(f'Number assigned congruent: {INFO_num_assigned_congruent}')
print(f'Number assigned non-congruent: {INFO_num_assigned_noncongruent}')
##/
###/

### Flatten tree and try to fit ambig SNPs where possible. Add a node for such SNPs
## Flatten tree
tree_flatten = copy.deepcopy(tree)
branch_trimmed = None
while branch_trimmed == None or branch_trimmed == True:
    branch_trimmed = False # toggle to true if we trimmed a branch
    for branch_node in tree_flatten.get_nonterminals():
        # get branch name
        branch_node_id = branch_node.name
        #/
        # skip branch if it does not have "Inner" as branch identifier
        if branch_node_id == None or branch_node_id.find('Inner') == -1: continue
        #/
        # If this node did not have congruent SNPs, then remove it and put all its childrens at the parent (P -> cur_node -> C ===> P -> C)
        if not branch_node_id in branch_stable_and_congruent_datasets and not branch_node_id in branch_instable_and_congruent_datasets_noStablesAtNode:
            # get path to root and skip this branch_node if it does not have a parent
            branch_node_path = tree_flatten.get_path(branch_node_id)
            if branch_node_path == None or len(branch_node_path) < 2: continue
            #/
            # INFO print
            print(f'Begin deletion of {branch_node_id}')
            #/
            # get parent from path
            parent_node = branch_node_path[-2] # Last element is current node, next-to-last element is the parent of current node
            #/
            # add childrens from current node to parent
            for child in branch_node.clades:
                parent_node.clades.append(child)
            #/
            # remove self from parent
            for rm_idx,child_of_parent in enumerate(parent_node.clades):
                if child_of_parent.name == branch_node.name:
                    del parent_node.clades[rm_idx]
                    break
            #/
            # remove childrens from self
            branch_node.clades = []
            #/
            # Set branch_trimmed toggle to trigger re-traversal of tree
            branch_trimmed = True
            #/
            # INFO-print
            print(f'Removed node {branch_node_id}')
            #/
if 0 and 'IDE':
    Phylo.draw(tree, do_show=False)
    Phylo.draw(tree_flatten, do_show=False)
##/

## Fit ambig SNPs, try to add nodes where they are congruent
def getTreeNode(tree_to_traverse=None,node_name=''):
    target_node = None
    for clade in tree_to_traverse.find_clades():
        if clade.name == node_name:
            target_node = clade
            break
    return target_node

def find_branches_to_nodify(tree_to_traverse=None,flatTree_levels_nodes_datasets=None,flatTree_nodes_datasets=None):
    return_data = {} # will return a "branch":branch_to_nodify and "branchify_with":[childs to put under new branch]
    for root_dist in sorted(flatTree_levels_nodes_datasets):
        for branch_node,node_datasets in flatTree_levels_nodes_datasets[root_dist].items():
            # get tree-object for current node
            branch_node_obj = getTreeNode(tree_to_traverse=tree_to_traverse,node_name=branch_node)
            #/
            # Check if current branch node has more than x2 childs (i.e. we might be able to add a subnode with congruent datasets)
            branch_node_childs = []
            for child_node in branch_node_obj.clades:
                if not child_node.name in datasets_in_tree: # skip if it is a leaf (dataset)
                    branch_node_childs.append(child_node)
            
            if len(branch_node_childs) > 2:
                ## Test all possible combinations (biggest group first)
                combination_sizes = list(range(2, len(branch_node_childs)+1)) # +1 since zero-based
                for combination_size in combination_sizes[::-1]:
                    combinations = list(itertools.combinations(branch_node_childs, combination_size))
                    for combination in combinations:
                        # get datasets at combination
                        combination_datasets = set()
                        combination_branch_names = []
                        for tmp_branch in combination:
                            tmp_branch_name = tmp_branch.name
                            
                            tmp_branch_datasets = flatTree_nodes_datasets[tmp_branch_name]
                            
                            # save1
                            combination_datasets.update(tmp_branch_datasets)
                            #/
                            # save2
                            combination_branch_names.append(tmp_branch_name)
                            #/
                        #/
                        # check if this combination matches an ambig snp
                        for ref_coord,msa_nucl_samples in ambig_snp_MSA_nucls_refNuclDiscarded.items():
                            for msa_nucl,samples in ambig_snp_MSA_nucls_refNuclDiscarded[ref_coord].items():
                                
                                # Check if this MSA-nucleotide was already assigned in the non-flattened tree: then skip it
                                if ref_coord in ambig_snps_assignment_status:
                                    if msa_nucl in ambig_snps_assignment_status[ref_coord]:
                                        if ambig_snps_assignment_status[ref_coord][msa_nucl]['status'] == 'congruent':
                                            continue
                                #/

                                # check if congruent
                                if samples == combination_datasets:
                                    # Return this combination
                                    return_data['branch'] = branch_node
                                    return_data['branchify_with'] = combination_branch_names
                                    return return_data
                                    #/
                                #/
                        #/
                ##/
            #/
    # will return empty if we didnt return in loops
    return return_data
    #/

tree_nodifyed = copy.deepcopy(tree_flatten)

## PRE: need to put a branch-node dummy infront of any datasets that doesnt have it (only those with singleton/unique SNP got it from previous script)
print('Adding "connector-nodes" before leafs in the tree that do not already have it.')
# Get highest node-branch ID (I.e. "InnerX"). Will be used to name the connector-nodes
highest_node_branch_enum = 0
for node in tree_nodifyed.find_clades():
    if not node.is_terminal():
        if type(node.name) == str and node.name.find('Inner') != -1:
            node_branch_enum = int(node.name.replace('Inner',''))
            if node_branch_enum > highest_node_branch_enum:
                highest_node_branch_enum = node_branch_enum
#/
# For each leaf-node, add a "connector-node"
new_node_branch_enum = highest_node_branch_enum + 1 # start adding new nodes at this enumerate
for enum,leaf_node in enumerate(tree_nodifyed.get_terminals()):
    # Get the parent of the leaf-node, skip the leaf-node if it does not have a parent (i.e. the case for "Reference-sequence" which is at the root)
    path_to_root = tree_nodifyed.get_path(leaf_node)
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
#/
##/

ambig_snp_fitted = None
new_nodes_names = set() # keep track of new nodes added (1-based names)
IDE_iters = 0
IDE_maxiters = 100 # maximum number of times to iterate over tree and branchify it
branches_expanded = [] # keep track of which nodes were already expanded. Do not expand them again
while ambig_snp_fitted == None or ambig_snp_fitted == True:
    # set toggle state (toggle to true if we made a new node with ambig snps a branch)
    ambig_snp_fitted = False
    #/
    ## Traverse tree levels and assign datasets at each level->branch
    flatTree_levels_nodes_datasets = {} # treeLevel (distance to root) -> node -> datasets
    flatTree_nodes_datasets = {} # node -> datasets
    for dataset in datasets_in_tree:
        dataset_ancestor_path = tree_nodifyed.get_path(dataset)
        for root_dist,node in enumerate(dataset_ancestor_path):
            if node.name == dataset: continue #skip if node is self (this entry should appears to always be last but I implement this in a general way regardless)
            
            # save1
            if not root_dist in flatTree_levels_nodes_datasets:                     flatTree_levels_nodes_datasets[root_dist] = {}
            if not node.name in flatTree_levels_nodes_datasets[root_dist]:          flatTree_levels_nodes_datasets[root_dist][node.name] = set()
            flatTree_levels_nodes_datasets[root_dist][node.name].add(dataset)
            #/
            
            # save2
            if not node.name in flatTree_nodes_datasets:            flatTree_nodes_datasets[node.name] = set()
            flatTree_nodes_datasets[node.name].add(dataset)
            #/
            
    print(f'Assigned N={len(flatTree_levels_nodes_datasets)} levels an N={len(flatTree_nodes_datasets)} branch nodes')
    ##/
    ## Check if there is any set of branches that can be make into a new node (i.e. congruent datasets at that new branch)
    # get a node to branchify
    branches_to_nodify = find_branches_to_nodify(tree_to_traverse=tree_nodifyed,
                                                 flatTree_levels_nodes_datasets=flatTree_levels_nodes_datasets,
                                                 flatTree_nodes_datasets=flatTree_nodes_datasets)
    #/
    # make branch
    if branches_to_nodify:
        # WARNING: THIS MAY NOT WORK AS INTENDED! check so that this branch was not previously expanded
        if branches_to_nodify['branchify_with'] in branches_expanded:
            print('This branch was already nodifyed, skipping')
            #print('WARNING: THIS MAY NOT WORK AS INTENDED!')
            continue
        branches_expanded.append(branches_to_nodify['branchify_with'])
        #/
        
        print(f'Found branch to nodify: {branches_to_nodify}')
        # make branch (P -> C ===> P -> new_branch -> C)
        # get parent/child names
        parent_branch_name = branches_to_nodify['branch']
        child_names = branches_to_nodify['branchify_with']
        #/
        # get parent/child objects
        parent_branch = getTreeNode(tree_to_traverse=tree_nodifyed,node_name=parent_branch_name)
        childs = []
        for child_name in child_names:
            childs.append(getTreeNode(tree_to_traverse=tree_nodifyed,node_name=child_name))
        #/
        # remove childs from parent node
        child_removal_done = None
        while child_removal_done == None or child_removal_done == False:
            child_removal_done = True # assume we will not remove any child
            for rm_idx,child_of_parent in enumerate(parent_branch.clades):
                if child_of_parent.name in child_names:
                    del parent_branch.clades[rm_idx]
                    child_removal_done = False # set toggle when removing a child
                    break # trigger iter reset
        #/
        # make new node
        new_node_name = f'Nodify:{len(new_nodes_names)+1}' # 1-based
        new_nodes_names.add(new_node_name) # add to tracker
        
        new_node_obj = Phylo.BaseTree.Clade(name=new_node_name,branch_length=0.1337) # dummy-set a branch length. will update later
        new_node_obj.clades = childs
        parent_branch.clades.append( new_node_obj )
        #/
        # set toggle state
        ambig_snp_fitted = True
        #/
    #/
    ##/
    
    IDE_iters += 1
    if IDE_iters > IDE_maxiters:
        print('maxiter reached, breaking')
        break

if 0 and 'IDE':
    fig, ax = plt.subplots(figsize=(15,15))
    Phylo.draw(tree_nodifyed, axes=ax, do_show=False)
##/

## After adding nodes to tree where possible, find congruent SNPs
nodifyTree_congruent_branch_ambig_snps = {} # branch_ID -> snp_pos -> msa_nucl -> samples
nodifyTree_ambig_snps_assignment_status = {} # keep track of which ref_pos->msa_nucl were congruent in tree
for branch_node,node_datasets in flatTree_nodes_datasets.items():
    for ref_coord,msa_nucl_samples in ambig_snp_MSA_nucls_refNuclDiscarded.items():
        for msa_nucl,samples in ambig_snp_MSA_nucls_refNuclDiscarded[ref_coord].items():
            
            # Check if this MSA-nucleotide was already assigned in the non-flattened tree: then skip it
            if ref_coord in ambig_snps_assignment_status:
                if msa_nucl in ambig_snps_assignment_status[ref_coord]:
                    if ambig_snps_assignment_status[ref_coord][msa_nucl]['status'] == 'congruent':
                        continue
            #/
            # Check if this MSA-nucleotide was already assigned in a previous iteration: then skip it
            #/
            
            # init this ref_pos->nucl in assignment_status
            if not ref_coord in nodifyTree_ambig_snps_assignment_status:               nodifyTree_ambig_snps_assignment_status[ref_coord] = {}
            if not msa_nucl in nodifyTree_ambig_snps_assignment_status[ref_coord]:     nodifyTree_ambig_snps_assignment_status[ref_coord][msa_nucl] = {'samples':samples,'status':'incongruent'} # init just once, will overwrite after testing all nodes in tree if it was found congruent
            #/
            # check if congruent
            if samples == node_datasets:
                if not branch_node in nodifyTree_congruent_branch_ambig_snps:                   nodifyTree_congruent_branch_ambig_snps[branch_node] = {}
                if not ref_coord in nodifyTree_congruent_branch_ambig_snps[branch_node]:        nodifyTree_congruent_branch_ambig_snps[branch_node][ref_coord] = {}
                nodifyTree_congruent_branch_ambig_snps[branch_node][ref_coord][msa_nucl] = samples
                
                nodifyTree_ambig_snps_assignment_status[ref_coord][msa_nucl]['status'] = 'congruent' # overwrite incongruent
            #/
print(f'Found N={len(nodifyTree_congruent_branch_ambig_snps)} branch-nodes with congruent SNPs')
##/
###/
### Begin output
tree_to_output = copy.deepcopy(tree_nodifyed) # do a deepcopy backup before modifying label names
## Init output dir
if not os.path.exists(output_dir):          os.makedirs(output_dir)
##/

## Change node labels (empty if no congruent SNP and otherwise to the congruent SNPs)
# Change branch-node name
for node in tree_to_output.find_clades():
    if not node.is_terminal():
        node_name = node.name
        # determine name on format "totSnps (stable_congruent/ambig_congruent/ambig_nodify_congruent)"
        stable_congruent = 0
        if node_name in branch_stable_and_congruent_datasets:
            stable_congruent = len(branch_stable_and_congruent_datasets[node_name])
            
        ambig_congruent = 0
        if node_name in congruent_branch_ambig_snps:
            ambig_congruent = len(congruent_branch_ambig_snps[node_name])
        
        ambig_nodify_congruent = 0
        if node_name in nodifyTree_congruent_branch_ambig_snps:
            ambig_nodify_congruent = len(nodifyTree_congruent_branch_ambig_snps[node_name])
        
        tot_congruent = stable_congruent + ambig_congruent + ambig_nodify_congruent
        #/
        # set branch length to total number of SNPs
        print(f'Set branch length {node_name}={tot_congruent}')
        node.branch_length = tot_congruent
        #/
        if tot_congruent > 0:
            node.name = f'{tot_congruent} ({stable_congruent}/{ambig_congruent}/{ambig_nodify_congruent})'
        else:
            node.name = ''
#/
# Display the after-tree in figure
if 1 and 'plot tree?':
    def plotTree(tree_to_plot,saveLoc=None):
        fig, ax = plt.subplots(figsize=(15,15))
        Phylo.draw(tree_to_plot, axes=ax, do_show=False)
        
        ## reposition branch-node labels
        x_min, x_max = ax.get_xlim()
        x_diff = x_max-x_min
        y_min, y_max = ax.get_ylim()
        y_diff = y_max-y_min
        for text_obj in ax.texts:
            if text_obj.get_text().strip().split()[0].isdigit(): # only edit if the text is at a LABELLED branch-node (expect a digit)
                x, y = text_obj.get_position()
                text_obj.set_position([x-x_diff*0.015,y+y_diff*0.0075]) # move slightly
                text_obj.set_color('red')
                
        ##/
        ax.spines['right'].set_visible(False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        plt.title('NJ-Tree')
        plt.tight_layout()
        if saveLoc != None:
            plt.savefig(saveLoc)
    
    # save tree before label mod
    plotTree(tree_nodifyed,saveLoc = output_dir+'/'+'tree_before.pdf')
    #/
    # save tree with label mod
    plotTree(tree_to_output,saveLoc = output_dir+'/'+'tree_after.pdf')
    #/
#/
##/
## Write tree (Branch-names set only at congruent branch-nodes [=number of SNPs at node])
Phylo.write(tree_flatten, output_dir+'/'+'tree_flatten.nwk', 'newick')
Phylo.write(tree_nodifyed,output_dir+'/'+'tree_nodifyed_raw.nwk', 'newick')
Phylo.write(tree_to_output, output_dir+'/'+'tree_nodifyed.nwk', 'newick')
##/

## Write congruent SNPs (before modifying tree)
with open(output_dir+'/'+'SNPs_instable.tsv','w') as nf:
    # write header
    header = ['internal_branch_id','ref','ref_pos','ref_nucl','nucl','samples']
    nf.write('\t'.join(map(str,header))+'\n')
    #/
    #
    for branch_node,ref_coord_msa_nucl_samples in sorted(congruent_branch_ambig_snps.items(),key=lambda x: int(x[0].replace('Inner','')), reverse=True): # Try to traverse congruent branch-nodes sorted by their digit. Expected format: "Inner<num>"
        for ref_coord,nucl_samples in ref_coord_msa_nucl_samples.items():
            for nucl,samples in nucl_samples.items(): # nucl here refers to MSA-nucl, do not confuse with ref nucl
                ref_nucl = ref_seq[int(ref_coord)]
                
                writeArr = [branch_node,ref_fasta,ref_coord,ref_nucl,nucl,','.join(sorted(samples))]
                nf.write('\t'.join(map(str,writeArr))+'\n')
    #/
##/

## Write fitted SNPs
with open(output_dir+'/'+'SNPs_fitted.tsv','w') as nf:
    # write header
    header = ['internal_branch_id','ref','ref_pos','ref_nucl','nucl','samples']
    nf.write('\t'.join(map(str,header))+'\n')
    #/
    #
    for branch_node,ref_coord_msa_nucl_samples in sorted(nodifyTree_congruent_branch_ambig_snps.items(),key=lambda x: int(x[0].replace('Nodify:','').replace('Inner','')), reverse=True): # Try to traverse congruent branch-nodes sorted by their digit. Expected format: "Inner<num>"
        for ref_coord,nucl_samples in ref_coord_msa_nucl_samples.items():
            for nucl,samples in nucl_samples.items(): # nucl here refers to MSA-nucl, do not confuse with ref nucl
                ref_nucl = ref_seq[int(ref_coord)]
                
                writeArr = [branch_node,ref_fasta,ref_coord,ref_nucl,nucl,','.join(sorted(samples))]
                nf.write('\t'.join(map(str,writeArr))+'\n')
    #/
##/

## Write incongruent SNPs (not fitted to tree)
# determine which SNPs are still incongruent after running this script
noncongruent_snps_samples_after_fitting = {} # ref_coord -> nucl -> samples
for ref_coord,nucls_samples in noncongruent_snps_samples.items():
    new_nucls_samples = {}
    for nucl,samples in nucls_samples.items():
        for sample in samples:
            # check if this position and samples were fitted
            sample_was_fitted = False
            if ref_coord in nodifyTree_ambig_snps_assignment_status:
                if nucl in nodifyTree_ambig_snps_assignment_status[ref_coord]:
                    if sample in nodifyTree_ambig_snps_assignment_status[ref_coord][nucl]['samples'] and nodifyTree_ambig_snps_assignment_status[ref_coord][nucl]['status'] == 'congruent':
                        sample_was_fitted = True
            #/
            # save
            if not sample_was_fitted:
                if not ref_coord in noncongruent_snps_samples_after_fitting:                    noncongruent_snps_samples_after_fitting[ref_coord] = {}
                if not nucl in noncongruent_snps_samples_after_fitting[ref_coord]:              noncongruent_snps_samples_after_fitting[ref_coord][nucl] = set()
                noncongruent_snps_samples_after_fitting[ref_coord][nucl].add(sample)
            #/
#/
# write
with open(output_dir+'/'+'SNPs_nonCongruent.tsv','w') as nf:
    # write header
    header = ['ref_pos','ref_nucl','A_num_samples','T_num_samples','G_num_samples','C_num_samples','missing_num_samples','A_samples','T_samples','G_samples','C_samples','missing_samples']
    nf.write('\t'.join(map(str,header))+'\n')
    #/
    # write rows
    for ref_coord,nucls_samples in noncongruent_snps_samples_after_fitting.items():
        # get ref nucl
        ref_nucl = ref_seq[int(ref_coord)]
        #/
        # get samples for each nucleotide   
        A_samples = [] if not 'A' in nucls_samples else list(nucls_samples['A'])
        T_samples = [] if not 'T' in nucls_samples else list(nucls_samples['T'])
        G_samples = [] if not 'G' in nucls_samples else list(nucls_samples['G'])
        C_samples = [] if not 'C' in nucls_samples else list(nucls_samples['C'])
        missing_samples = [] if not '-' in nucls_samples else list(nucls_samples['-'])
        #/
        ## add refcounts and gaps
        msa_idx = ref_coords[ref_coord]
        for dataset,msa_seq in msa_seqs.items():
            dataset_nucl = msa_seq[msa_idx]
            
            if dataset_nucl == 'A' and not dataset in A_samples:            A_samples.append(dataset)
            if dataset_nucl == 'T' and not dataset in T_samples:            T_samples.append(dataset)
            if dataset_nucl == 'G' and not dataset in G_samples:            G_samples.append(dataset)
            if dataset_nucl == 'C' and not dataset in C_samples:            C_samples.append(dataset)
            if dataset_nucl == '-' and not dataset in missing_samples:      missing_samples.append(dataset)
        #/
        ##/
        # determine nucl counts
        A_num_samples = len(A_samples)
        T_num_samples = len(T_samples)
        G_num_samples = len(G_samples)
        C_num_samples = len(C_samples)
        missing_num_samples = len(missing_samples)
        #/
        # write row
        writeArr = [ref_coord,ref_nucl,A_num_samples,T_num_samples,G_num_samples,C_num_samples,missing_num_samples,','.join(A_samples),','.join(T_samples),','.join(G_samples),','.join(C_samples),','.join(missing_samples)]
        nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
    #/
#/
##/