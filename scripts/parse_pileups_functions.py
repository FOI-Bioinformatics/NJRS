
"""
This is a draft-script of function to find if an input of samples conflict their relations in a provided tree
"""

import os
import sys
from matplotlib import pyplot as plt
import re
import ete3

# Import misc functions
from functions import *
#/

### Make tree and get all leaf-nodes for each branch-node
def construct_tree_get_nodes_childs(newick_tree,tree_root_setName=None):
    ## Import tree
    tree = ete3.Tree(newick_tree.format('newick'))
    if tree_root_setName:
        root_node = tree.search_nodes(name=tree_root_setName)[0] # if this rows fail it is an indication that the tree root name was not found
        tree.set_outgroup(root_node)
    ##/
    
    ## For each node in the tree, calculate all samples it holds
    # pre-flight: assign a name for un-named nodes
    leaf_node_enum = 0
    branch_node_enum = 0
    for node in tree.traverse():
        # for leaves that do not have a name, add an arbitrary one. However, we expect all nodes to havea a name, matching the SNP dataset
        if node.is_leaf():
            if not node.name:
                node.name = 'leaf:'+str(leaf_node_enum)
                leaf_node_enum += 1
                print('WARNING: leaf node had no name! Asssigning arbirary name: '+str(node.name))
        #/
        # add branch-name for branch-nodes with empty names
        if not node.is_leaf():
            if not node.name:
                node.name = 'branch:'+str(branch_node_enum)
                branch_node_enum += 1
        #/
    #/
    # for each branch-node in the tree, save all of its leaf-nodes
    tree_nodes_childs = {}
    for node in tree.traverse():
        # skip row if it is not a branch-node
        if node.is_leaf(): continue
        #/
        # get all leaves for current branch-node
        leaf_nodes = set()
        for desc_node in node.get_descendants():
            if desc_node.is_leaf():
                leaf_nodes.add(desc_node.name)
        #/
        # save
        tree_nodes_childs[node.name] = leaf_nodes
        #/
    #/
    return tree_nodes_childs,tree
    ##/


### Check conflicts relative to the tree
## Traverse SNPs and find conflicts relative to tree
def check_sample_tree_conflict(samples_with_var,tree_nodes_childs,return_congruent_node=False):
    """
    Given a set of samples, determine if the samples are neighbours in the tree (or if they have a conflict accroding to the tree)
    Names are matched by exact-matching
    Optionally, if the input samples are found to be congruent in the tree, their node is returned by "return_congruent_node=True"
    """
    # check input
    if not type(samples_with_var) == set:
        try:
            samples_with_var = set(samples_with_var)
        except:
            sys.exit('This function expects a set as input (or an input that can be converted to a set')
    #/
    # check branches in the tree
    branch_node_found = False
    congruent_branch_node = None
    for branch_node,leaf_nodes in tree_nodes_childs.items():
        if leaf_nodes == samples_with_var:
            branch_node_found = True
            congruent_branch_node = branch_node
            break #break on first
    #/
    # Check if we want to return congruent node
    if return_congruent_node:
        return congruent_branch_node
    #/
    # determine if they havea conflict. If we found a branch with these samples as childs (and no other childs on that branch), then there is no conflict
    samples_have_tree_conflict = False
    if not branch_node_found:
        samples_have_tree_conflict = True
    #/
    # check if only one samples was input. Then there is nothing to check, there is no conflict
    if len(samples_with_var) == 1:
        samples_have_tree_conflict = False
    #/
    return samples_have_tree_conflict
##/
