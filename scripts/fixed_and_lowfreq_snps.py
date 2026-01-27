#!/usr/bin/env python3

from matplotlib import pyplot as plt
import sys
import os
import argparse
import pickle
import gzip
import math

# tree-related functions (identify samples with e.g. shared snp that are not relatives in the tree)
import parse_pileups_functions
#/
# Import misc functions
from functions import *
#/

### Input parse
## Setup parser
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

grp_IO = parser.add_argument_group('Input and output')
grp_IO.add_argument('-i','--input',required=True,help='Path to directory of tables with position base counts')
grp_IO.add_argument('-o','--output',required=True,help='Path to output-directory')
grp_IO.add_argument('-g','--genome',required=True,help='Path to reference-genome')

grp_classification = parser.add_argument_group('Classification and parsing')
grp_classification.add_argument('--snp_min_count',type=int,required=False,default=1,help='Minimum count of reads supportive of SNPs to import position and test it across datasets (default:1)')
grp_classification.add_argument('--pos_min_cov',type=int,required=False,default=1,help='Minimum coverage of SNP position to consider a SNP. Positions below the input value are imported and classified as "low coverage" (default:1)')
grp_classification.add_argument('--lowcount_thresh',type=float,default=2,help='Consider SNPs with counts below input value as "low-count" (default:2)')
grp_classification.add_argument('--lowfreq_frac_thresh',type=float,default=0.2,help='Consider SNPs with fractional-count below input value as "low-frequent" (default:0.2)')
grp_classification.add_argument('--fixed_min_count',type=int,default=5,help='Minimum number of counts to consider a fixed SNP (default:5)')
grp_classification.add_argument('--fixed_min_frac',type=float,default=0.7,help='Minimum fraction of counts to consider a fixed SNP (default:0.7)')

grp_msa = parser.add_argument_group('MSA')
grp_msa.add_argument('--ignore_same',action='store_true',required=False,default=False,help='If specified, will not output nucleotide-positions where all datasets have the same value (default: not set)')
grp_msa.add_argument('--no_N',action='store_true',required=False,default=False,help='If specified, will remove nucleotide-positions with conflicts (default: set to N)')
#grp_msa.add_argument('--no_noise',action='store_true',required=False,default=False,help='[Not implemented] If specified, will skip fixed nucleotides where "unstable" classifications exist for the same sample. Example: sample={nucl_fixed=G(61),nucl_lowcount=A(2)} will be discarded (default: output all fixed nucleotide positions)')
grp_msa.add_argument('--allow_lowcov_fraction',type=float,required=False,default=None,help='For MSA output, include positions where up to this fraction of datasets (rounded up) have no/low coverage. No/low coverage datasets will get a "-" output [may or may not be compatible will complete pipe] (default:Not set)')
##/

# Parse commandline
args = parser.parse_args()
#/
# Parse variables from command
base_counts_tables_path = args.input
output_dir = args.output
refgenome_path = args.genome

snp_min_count = args.snp_min_count
pos_min_cov = args.pos_min_cov

lowcount_thresh = args.lowcount_thresh
lowfreq_frac_thresh = args.lowfreq_frac_thresh

fixed_min_count = args.fixed_min_count
fixed_min_frac = args.fixed_min_frac

msa_ignore_same = args.ignore_same
msa_no_N_out = args.no_N
#msa_no_noise = args.no_noise
msa_allow_nocov_lowcov_fraction = args.allow_lowcov_fraction
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

### Parse refgenome
print('Importing reference-sequence',flush=True)
ref_seqs = importReadSeqsFasta(refgenome_path)
ref_name = list(ref_seqs.keys())[0]
ref_seq = ref_seqs[ref_name]
###/

### Parse basecounts
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


## Pre-flight: scan all input and find which bases are modified in any sample. Then we will import only those positions, to save memory.
print('Pre-flight: parsing input-files to determine which positions to import',flush=True)
samples_found = set()
poses_to_import = {} # pos -> samples with any variant
skip_reasons_count = {}
for fileEnum,file_ in enumerate(os.listdir(base_counts_tables_path)):
    if fileEnum % max(int(len(os.listdir(base_counts_tables_path))/10),1) == 0: print('Processed ' + str(fileEnum) + ' files',flush=True)
    sample = os.path.splitext(file_)[0]
    samples_found.add(sample)
    with open_file(base_counts_tables_path+'/'+file_) as f:
        for enum,line in enumerate(f):
            if enum == 0: continue # skip header line
            line = line.strip('\n')
            line = line.split('\t')
            
            # parse line
            pos,num_A,num_T,num_G,num_C,num_N,is_empty,coverage,consensus,is_ref,consensus_is_ref = line
            pos = int(pos)
            num_A = int(num_A)
            num_T = int(num_T)
            num_G = int(num_G)
            num_C = int(num_C)
            num_N = int(num_N)
            coverage = int(coverage)
            #/
            
            # check if there is any variation at position
            has_value = []
            if num_A >= snp_min_count: has_value.append('A')
            if num_T >= snp_min_count: has_value.append('T')
            if num_G >= snp_min_count: has_value.append('G')
            if num_C >= snp_min_count: has_value.append('C')
            if num_N >= snp_min_count: has_value.append('N')
            #/
            
            ## Check if skip pos
            skip_pos = True # assume true, prove otherwise below
            
            # check if coverage is 0 or below pos_min_cov, then keep position (do not want to exclude positions with lacking data)
            if coverage == 0 or coverage < pos_min_cov:
                skip_pos = False
            #/
            # check if multiple values exist
            if len(has_value) > 1: # multiple values means at least one is different from reference
                skip_pos = False
            
            elif len(has_value) == 1: # if one value, then test so it is not hte same as reference
                # keep position only if the count is NOT the reference nucleotide
                if not ref_seq[pos] in has_value:
                    skip_pos = False
                #/
            #/
            ##/
            # Add position as "to be imported" if successful
            if not skip_pos:
                init(poses_to_import,pos,[])
                poses_to_import[pos].append(sample)
            #/
print('Finalized pre-flight processing of '+str(fileEnum+1)+' files',flush=True)
##/
## Import data at positions selected above
pos_data = {} # pos -> data {A/T/G/C -> sample -> counts} , {coverage} , ref_nucl
for fileEnum,file_ in enumerate(os.listdir(base_counts_tables_path)):
    if fileEnum % max(int(len(os.listdir(base_counts_tables_path))/10),1) == 0: print('Processed ' + str(fileEnum) + ' files',flush=True)
    sample = os.path.splitext(file_)[0]
    with open_file(base_counts_tables_path+'/'+file_) as f:
        for enum,line in enumerate(f):
            if enum == 0: continue # skip header line
            line = line.strip('\n')
            line = line.split('\t')
            
            # parse line
            pos,num_A,num_T,num_G,num_C,num_N,is_empty,coverage,consensus,is_ref,consensus_is_ref = line
            pos = int(pos)
            #/
            # Check if position should be parsed
            if not pos in poses_to_import: continue
            #/
            # convert some variables to int
            num_A,num_T,num_G,num_C,num_N,is_empty,coverage,is_ref,consensus_is_ref = map(int,[num_A,num_T,num_G,num_C,num_N,is_empty,coverage,is_ref,consensus_is_ref])
            #/
            
            ## Save data at pos
            # save nucleotide counts
            for nucl,count in ( ['A',num_A],['T',num_T],['G',num_G],['C',num_C],['N',num_N] ):
                SK = 'counts'
                init(pos_data,pos,{})
                init(pos_data[pos],SK,{})
                init(pos_data[pos][SK],nucl,{})
                pos_data[pos][SK][nucl][sample] = count
            #/
            # save total counts per sample (including N)
            for nucl in pos_data[pos]['counts']:
                count = pos_data[pos]['counts'][nucl][sample]
                SK = 'total_counts'
                init(pos_data,pos,{})
                init(pos_data[pos],SK,{})
                init(pos_data[pos][SK],sample,0)
                pos_data[pos][SK][sample] += count
            #/
            # save nucleotide fractions
            for nucl in pos_data[pos]['counts']:
                count = pos_data[pos]['counts'][nucl][sample]
                SK = 'samples_fractions'
                init(pos_data,pos,{})
                init(pos_data[pos],SK,{})
                init(pos_data[pos][SK],sample,{})
                
                frac = 0 # if totalcounts is 0, then we used default as 0
                if pos_data[pos]['total_counts'][sample] > 0:
                    frac = count / pos_data[pos]['total_counts'][sample]
                pos_data[pos][SK][sample][nucl] = frac
            #/
            # save position coverage
            SK = 'coverage'
            init(pos_data,pos,{})
            init(pos_data[pos],SK,{})
            pos_data[pos][SK][sample] = coverage
            #/
            # save position reference nucleotide
            SK = 'ref_nucl'
            init(pos_data,pos,{})
            pos_data[pos][SK] = ref_seq[pos]
            #/
            ##/
            
print('Finalized processing of '+str(fileEnum)+' files',flush=True)
##/
###/

### Sort imported data: (1) "stable and fixed SNPs" [across all samples], (2) "low-frequency and mixed-support SNPs" [in any sample], (3) "no/low coverage SNPs" [in any sample]
pos_samples_nucls_classified = {} # pos -> sample -> classification -> nucleotide+data
pos_classification_summary = {} # pos -> classification -> sample -> nucleotide+data
for pos,data in pos_data.items():
    ## Classify all "shaky-types-of-classifications"
    # Classify lowcount, lowfreq, fixed, mixed, and shaky
    samples_nucls_classified = {}
    for sample,nucl_fractions in data['samples_fractions'].items():
        ### Method: put all nucl-fractions > 0 in an array. Move away values that are classified as lowfreq. Then mova away those classified as fixed. If the remainder is at least of size two, then these are "mixed"
        ## init (parse fraction-counts greater than 0; i.e. skip if no counts)
        nucls_fracs = []
        for nucl,frac in nucl_fractions.items():
            if frac > 0:
                nucls_fracs.append([nucl,frac])
        ##/
        ## pull lowcounts
        nucls_lowcounts = []
        while True:
            resetIteration = False
            for enum,(nucl,frac) in enumerate(nucls_fracs):
                count = frac*data['total_counts'][sample]
                if count < lowcount_thresh:
                    nucls_lowcounts.append([nucl,frac,count])
                    del nucls_fracs[enum]
                    resetIteration = True
                    break
            
            if not resetIteration:
                break
        ##/
        ## pull lowfreqs
        nucls_lowfreq = []
        while True:
            resetIteration = False
            for enum,(nucl,frac) in enumerate(nucls_fracs):
                if frac < lowfreq_frac_thresh:
                    nucls_lowfreq.append([nucl,frac])
                    del nucls_fracs[enum]
                    resetIteration = True
                    break
            
            if not resetIteration:
                break
        ##/
        ## pull fixed (SAME AS REF)
        nucls_fixed_isRef = []
        while True:
            resetIteration = False
            for enum,(nucl,frac) in enumerate(nucls_fracs):
                if not nucl == data['ref_nucl']: continue # test only nucls that are same as ref
                count = frac*data['total_counts'][sample]
                if frac >= fixed_min_frac and count >= fixed_min_count:
                    nucls_fixed_isRef.append([nucl,frac])
                    del nucls_fracs[enum]
                    resetIteration = True
                    break
            
            if not resetIteration:
                break
        ##/
        ## pull fixed (OTHER THAN REF)
        nucls_fixed = []
        while True:
            resetIteration = False
            for enum,(nucl,frac) in enumerate(nucls_fracs):
                count = frac*data['total_counts'][sample]
                if frac >= fixed_min_frac and count >= fixed_min_count:
                    nucls_fixed.append([nucl,frac])
                    del nucls_fracs[enum]
                    resetIteration = True
                    break
            
            if not resetIteration:
                break
        ##/
        ## pull mixed/orphan-shaky-nucl
        nucls_mixed = []
        nucls_shaky = []
        if len(nucls_fracs) >= 2:
            nucls_mixed = nucls_fracs
            #/
        else:
            nucls_shaky = nucls_fracs
        ##/
        ###/
        
        # Classify nocov/lowcov
        lowcov = False
        cov = data['coverage'][sample]
        if cov == 0 or cov < pos_min_cov:
            lowcov = True
        #/
        
        # Save to outer
        samples_nucls_classified[sample] = {'nucls_lowcounts':nucls_lowcounts,'nucls_lowfreq':nucls_lowfreq,
                                            'nucls_fixed_isRef':nucls_fixed_isRef,'nucls_fixed':nucls_fixed,
                                            'nucls_mixed':nucls_mixed,'nucls_shaky':nucls_shaky,
                                            'lowcov':lowcov}
        #/
    #/
    ##/
    ## Save classifications to outer
    pos_samples_nucls_classified[pos] = samples_nucls_classified
    ##/
    ## Save sample classification summary at position to outer
    # compute summary
    classifications_per_sample = {} # classi -> sample
    for sample,classis in samples_nucls_classified.items():
        for classi,entries in classis.items():
            if entries:
                init(classifications_per_sample,classi,{})
                classifications_per_sample[classi][sample] = entries
    #/
    # save to outer
    pos_classification_summary[pos] = classifications_per_sample
    #/
    ##/
###/


### Compile summary of position classifications
positions_classifications_summary = {} # classification -> count
positions_classifications_summary_allReasonsAsJoinedKey = {}
for pos,ref_nucl in enumerate(ref_seq):
    ## check if position did not have a classification
    if not pos in pos_classification_summary:
        SK = 'no_classification'
        # save1
        if not SK in positions_classifications_summary:             positions_classifications_summary[SK] = 0
        positions_classifications_summary[SK] += 1
        #/
        # save2
        if not SK in positions_classifications_summary_allReasonsAsJoinedKey:             positions_classifications_summary_allReasonsAsJoinedKey[SK] = 0
        positions_classifications_summary_allReasonsAsJoinedKey[SK] += 1
        #/
        continue
    ##/
    ## Get classifications at pos and save
    classification_summary = pos_classification_summary[pos]
    # save1
    for classification in classification_summary:
        if not classification in positions_classifications_summary:     positions_classifications_summary[classification] = 0
        positions_classifications_summary[classification] += 1
    #/
    # save2
    classifications_present_ID = '||'.join(sorted(classification_summary.keys()))
    if not classifications_present_ID in positions_classifications_summary_allReasonsAsJoinedKey:         positions_classifications_summary_allReasonsAsJoinedKey[classifications_present_ID] = 0
    positions_classifications_summary_allReasonsAsJoinedKey[classifications_present_ID] += 1
    #/
    ##/
##/
###/

### Write position classifications as (1) full table, (2) summarized table [only classified positions]
if 1 and 'OUTPUT?':
    ## init fileopens
    nf_all = open(output_dir+'/'+'positions_all.tsv','w') # full table [including classified positions]
    nf_summarized = open(output_dir+'/'+'positions_classified.tsv','w') # summarized table [only classified position]
    nfs = [nf_all,nf_summarized]
    ##/
    ## Header
    # define header
    header = ['pos','ref_nucl','datatable'] # "base info". datatable => "position taken from pos_classification_summary" or "position taken from pos_data"
    nucl_fixed_header = ['fixed_consensus','fixed_A','fixed_T','fixed_G','fixed_C','fixed_consensus_samples','fixed_A_samples','fixed_T_samples','fixed_G_samples','fixed_C_samples']
    header += nucl_fixed_header
    classi_keys_to_check_ordered = ['nucls_fixed','nucls_fixed_isRef','lowcov','nucls_lowcounts','nucls_lowfreq','nucls_shaky','nucls_mixed']
    header += classi_keys_to_check_ordered
    #/
    # write (both files)
    for nf in nfs:
        nf.write('\t'.join(map(str,header))+'\n')
    #/
    ##/
    ## Data
    for pos,ref_nucl in enumerate(ref_seq):
        ## Check if current position has classifications
        if pos in pos_classification_summary:
            classification_summary = pos_classification_summary[pos]
            # bugcheck: check so ref_nucleotide is same in imported ref-file and annotated ref_nucl in pos_data
            if not ref_nucl == pos_data[pos]['ref_nucl']:
                print('FATAL: reference nucleotide in imported fasta does not match reported ref_nucl in input basecount files!')
                sys.exit()
            #/
            # add "base info"
            writeArr = [ pos,ref_nucl,'pos_classification_summary' ]
            #/
            # add "nucl_fixed_header" info
            nucl_samples = {'A':set(),'T':set(),'G':set(),'C':set()}
            CK = 'nucls_fixed'
            if CK in classification_summary:
                for sample,nucl_data in classification_summary[CK].items():
                    for nucl,frac in nucl_data:
                        if nucl in nucl_samples:
                            nucl_samples[nucl].add(sample)
            #@ write consensus of nucleotide-numbers
            consensus_nucl,consensus_samples = sorted(nucl_samples.items(),key=lambda x: len(x[1]), reverse=True)[0]
            consensus_num = len(consensus_samples)
            consensus_samples_writeReady = 0 # default
            if consensus_num > len(samples_found)/2: # if fixed nucleotide exist in more than half of imported samples, then consider it consensus
                writeArr.append(consensus_nucl)
                consensus_samples_writeReady = ', '.join(sorted(consensus_samples))
            else:
                writeArr.append(0)
            #@/
            #@ write nucleotide numbers
            for nucl,samples in nucl_samples.items():
                writeArr.append(len(samples))
            #@/
            #@ write consensus samples
            writeArr.append(consensus_samples_writeReady)
            #@/
            #@ write nucleotide samples
            for nucl,samples in nucl_samples.items():
                samples_writeReady = 0
                if samples:
                    samples_writeReady = ', '.join(sorted(samples))
                writeArr.append(samples_writeReady)
            #@/
            #/
            # add classi info
            for classi in classi_keys_to_check_ordered:
                if classi in classification_summary:
                    writeArr.append(len(classification_summary[classi]))
                else:
                    writeArr.append(0)
            #/
            # write (both files)
            for nf in nfs:
                nf.write('\t'.join(map(str,writeArr))+'\n')
            #/
            # prevent downstream code
            continue
            #/
        ##/
        ## Check if we want to output non-classified spots too ()
        # add "base info"
        writeArr = [ pos,ref_nucl,'ref_fasta' ]
        #/
        # write 0-entries for all nucl_fixed_header
        for i in nucl_fixed_header:
            writeArr.append(0)
        #/
        # write 0-entries for all classifications
        for classi in classi_keys_to_check_ordered:
            writeArr.append(0)
        #/
        # write (ONLY TO "all"-file)
        nf_all.write('\t'.join(map(str,writeArr))+'\n')
        #/
        # prevent downstream code
        continue
        #/
        ##/
    ##/
    ## Close files
    for nf in nfs:
        nf.close()
    ##/
###/

###### SECTION: Consensus sequences of "stable" positions (e.g. for MSA)
### Parse "stable" positions only (skip unstable positions = "if any sample has unstable classifications")
pos_fixed_nucls_samples_vals = {} # pos -> sample -> fixed_nucleotide (if present) otherwise reference-nucleotide
skipped_positions = set() # keep track of positions that are skipped due to unstable classifications
positions_status = {}
pos_lowcov_dataset_rescues = {} # pos -> dataset_with_lowcov [will add "-" for these in MSA later]
for pos,ref_nucl in enumerate(ref_seq):
    ## Check if current position has classifications
    # save INFO
    SK = 'num_positions'
    if not SK in positions_status:     positions_status[SK] = 0
    positions_status[SK] += 1
    #/
    if pos in pos_classification_summary:
        # Get classifications
        classification_summary = pos_classification_summary[pos]
        #/
        
        ## Check if skip position
        # save INFO
        SK = 'tested_positions'
        if not SK in positions_status:     positions_status[SK] = 0
        positions_status[SK] += 1
        #/
        
        skip_pos = False
        # Check if position has any samples with nucls_fixed (otherwise it is not an interesting position)
        if not 'nucls_fixed' in classification_summary:
            skip_pos = True
            # save INFO
            SK = 'no_fixed_nonRef_nucl'
            if not SK in positions_status:     positions_status[SK] = 0
            positions_status[SK] += 1
            #/
        #/
        # Check if multiple nucls_fixed are reported across samples
        fixed_nucls = {} # nucl -> sample
        if 'nucls_fixed' in classification_summary:
            for sample,nucls_data in classification_summary['nucls_fixed'].items():
                # bugcheck: we do not yet support multiple fixed nuccleotides reported at a single sample
                if len(nucls_data) != 1:
                    print('FATAL: found multiple fixed variants at sample. this is not yet supported. Terminating!',flush=True)
                    sys.exit()
                #/
                nucl = nucls_data[0][0]
                init(fixed_nucls,nucl,set())
                fixed_nucls[nucl].add(sample)
        
        if len(fixed_nucls) > 1:
            print('INFO: found multiple fixed nucleotides at position '+str(pos))
            for nucl,samples in fixed_nucls.items():
                print(nucl+'\t'+str(len(samples)))
            
            if 0 and 'HALT ON MULTIPLE FIXED NUCLEOTIDES ACROSS SAMPLES':
                skip_pos = True
                print('Skipping position',flush=True)
                # save INFO
                SK = 'fixed_multiple_variants'
                if not SK in positions_status:     positions_status[SK] = 0
                positions_status[SK] += 1
                #/
            else:
                print('Position is in use as of writing this.')
        #/
        ## Check if skip position (any sample has unstable classification here WITHOUT having a fixed_nucl [fixed+unstable = OK, unstable = notOK])
        # get samples with nucls_fixed or nucls_fixed_isRef
        samples_with_any_fixed_type = set()
        for classification in ('nucls_fixed','nucls_fixed_isRef',):
            if classification in classification_summary:
                samples_data = classification_summary[classification]
                for sample in samples_data:
                    samples_with_any_fixed_type.add(sample)
        #/
        # check if there are unstable classifications where not all of those "unstable" stamples also have a fixed classification at the position
        for classification,samples_data in classification_summary.items():
            if not classification in ('nucls_fixed','nucls_fixed_isRef',): # do not test the stable classifications ==> we test the unstable classifications
                if not set(samples_data).intersection(samples_with_any_fixed_type) == set(samples_data): # if we have an unstable classification for X samples, then check so that not all of those X samples also have a stable(fixed) classification
                    skip_pos = True
                    # save INFO
                    SK = 'unstable_at_fixed'
                    if not SK in positions_status:     positions_status[SK] = 0
                    positions_status[SK] += 1
                    #/
        #/
        # [rescue] check if user wants to allow X number of lowcov/nocov to be present. Consider these as with the "ref" nucleotide
        if msa_allow_nocov_lowcov_fraction != None:
            # get datasets with nocov/lowcov and dirty-add them to checkpoint "samples_with_any_fixed_type"
            datasets_with_lowcov_or_nocov = set()
            for classi_to_check in ('nucls_lowcounts','lowcov',):
                if classi_to_check in classification_summary:
                    for dataset in classification_summary[classi_to_check]:
                        # Only add the dataset if it does not have a fixed or shaky nucl. then it cannot be a nucl_lowcounts
                        skip_dataset = False
                        for CK in ['nucls_fixed','nucls_fixed_isRef','nucls_lowfreq','nucls_shaky','nucls_mixed']:
                            if CK in classification_summary and dataset in classification_summary[CK]:
                                skip_dataset = True
                        #/
                        # save
                        if not skip_dataset:
                            datasets_with_lowcov_or_nocov.add(dataset)
                        #/
            #/
            # check if the number of datasets with nocov/lowcov is within user-specified range
            too_many_lowcov_nocov_num = math.ceil(len(samples_found)*msa_allow_nocov_lowcov_fraction)
            if len(datasets_with_lowcov_or_nocov) > 0 and len(datasets_with_lowcov_or_nocov) <= too_many_lowcov_nocov_num:
                if samples_with_any_fixed_type.union(datasets_with_lowcov_or_nocov) == samples_found:
                    skip_pos = False # do not skip this pos anymore
                    fixed_nucls['-'] = datasets_with_lowcov_or_nocov # save lowcov/nocov samples under "-" at fixed_nucls as rescue
                    
                    # save to outer
                    pos_lowcov_dataset_rescues[pos] = datasets_with_lowcov_or_nocov
                    #/
                    # save INFO
                    SK = 'lowcov_nocov_rescued'
                    if not SK in positions_status:     positions_status[SK] = 0
                    positions_status[SK] += 1
                    #/
            #/
        #/
        ##/
        # Check if skip
        if skip_pos:
            skipped_positions.add(pos)
            # save INFO
            SK = 'skipped_positions'
            if not SK in positions_status:     positions_status[SK] = 0
            positions_status[SK] += 1
            #/
            continue
        #/
        ##/
        
        ## Assign nucleotides for each sample at this position
        samples_added = set()
        # add samples with fixed nucl
        for nucl,samples in fixed_nucls.items():
            for sample in samples:
                init(pos_fixed_nucls_samples_vals,pos,{})
                pos_fixed_nucls_samples_vals[pos][sample] = nucl
                samples_added.add(sample)
        #/
        # add remainder with reference-nucleotide
        for sample in samples_found:
            if not sample in samples_added:
                pos_fixed_nucls_samples_vals[pos][sample] = ref_nucl
                samples_added.add(sample)
        #/
        # save INFO
        SK = 'pass1_positions'
        if not SK in positions_status:     positions_status[SK] = 0
        positions_status[SK] += 1
        #/
        ##/
###/

### Make pseudo-sequences per sample of stable positions
samples_pseudoseqs = {} # sample -> pseudo-seq
samples_pseudoseqs_poses = [] # sample -> ref_coord per position in pseudo-seqs
INFO_num_is_same_nucl = 0 # count how many positions have same values across all datasets
INFO_num_is_same_nucl_isRef = 0 # count how many positions have same values across all datasets AND is the same nucleotide as the reference
pos_fixed_nucls_samples_vals_allTheSame = {} # store positions where all samples had the same nucl
for pos,samples_nucl in pos_fixed_nucls_samples_vals.items():
    ## Get number of nucls reported at pos
    sample_nucls_counts = {}
    for sample,nucl in samples_nucl.items():
        init(sample_nucls_counts,nucl,0)
        sample_nucls_counts[nucl] += 1
        
    if len(sample_nucls_counts) == 1:
        INFO_num_is_same_nucl += 1
    if len(sample_nucls_counts) == 1 and list(sample_nucls_counts)[0] == ref_seq[pos]:
        INFO_num_is_same_nucl_isRef += 1
    ##/
    ## Check if skip positions where all samples have the same value (dont count "missing" as a value, i.e. missing+allsame is considered as allsame)
    if len(set(sample_nucls_counts).difference('-')) == 1 and msa_ignore_same:
        # save INFO
        SK = 'all_samples_same'
        if not SK in positions_status:     positions_status[SK] = 0
        positions_status[SK] += 1
        #/
        # save2 (output to "rescue" position list)
        pos_fixed_nucls_samples_vals_allTheSame[pos] = samples_nucl
        #/
        continue
    ##/
    
    for sample,nucl in samples_nucl.items():
        init(samples_pseudoseqs,sample,[])
        samples_pseudoseqs[sample].append(nucl)
    
    # save INFO
    SK = 'pass2_positions'
    if not SK in positions_status:     positions_status[SK] = 0
    positions_status[SK] += 1
    #/
    
    samples_pseudoseqs_poses.append(pos)

pseudoseqs_samples = {} # pseudo-seq -> sample (find how many samples have identical sequences)
pseudoseqs_lens = 0 # pseudo-seq -> seq_len
for sample,seq in samples_pseudoseqs.items():
    seq = ''.join(seq)
    init(pseudoseqs_samples,seq,[])
    pseudoseqs_samples[seq].append(sample)
    
    pseudoseqs_lens = len(seq)

# write info to file
with open(output_dir+'/'+'msa_seq_info.fasta','w') as nf:
    for writeData in ( ['msa_ignore_same',msa_ignore_same] , ['INFO_num_is_same_nucl',INFO_num_is_same_nucl] , ['INFO_num_is_same_nucl_isRef',INFO_num_is_same_nucl_isRef] ):
        nf.write('\t'.join(map(str,writeData))+'\n')
#/
print('Amongst samples, found '+str(len(pseudoseqs_samples))+' sequence variants of length '+str(pseudoseqs_lens),flush=True)
print('Positions where all datasets have the same nucleotide: '+str(INFO_num_is_same_nucl)+', same nucleotide and is same as reference: '+str(INFO_num_is_same_nucl_isRef))
###/

### Output some stats
## positions_classifications_summary
with open(output_dir+'/'+'pos_classification_summary.tsv','w') as nf:
    # write header
    header = ['classification','abundance']
    nf.write('\t'.join(map(str,header))+'\n')
    #/
    # write data
    for classification,abundance in sorted(positions_classifications_summary.items(),key=lambda x: x[1], reverse=True):
        writeArr = [classification,abundance]
        nf.write('\t'.join(map(str,writeArr))+'\n')
    #/
##/

## positions_classifications_summary_allReasonsAsJoinedKey
with open(output_dir+'/'+'pos_classification_detailed.tsv','w') as nf:
    # write header
    header = ['classifications','abundance']
    nf.write('\t'.join(map(str,header))+'\n')
    #/
    # write data
    for classifications,abundance in sorted(positions_classifications_summary_allReasonsAsJoinedKey.items(),key=lambda x: x[1], reverse=True):
        writeArr = [classifications,abundance]
        nf.write('\t'.join(map(str,writeArr))+'\n')
    #/
##/

## positions_status
with open(output_dir+'/'+'pos_classification_final.tsv','w') as nf:
    # write header
    header = ['description','abundance']
    nf.write('\t'.join(map(str,header))+'\n')
    #/
    # write data
    for description,abundance in sorted(positions_status.items(),key=lambda x: x[1], reverse=True):
        writeArr = [description,abundance]
        nf.write('\t'.join(map(str,writeArr))+'\n')
    #/
##/
###/

### Write slim seqs
if 1 and 'write pseudo seqs?':
    # Write seq
    with open(output_dir+'/'+'msa_slim.fasta','w') as nf:
        for sample,seq in samples_pseudoseqs.items():
            seq = ''.join(seq)
            nf.write('>'+sample+'\n'+seq+'\n')
    #/
    # Write positions of SNPs in reference
    with open(output_dir+'/'+'msa_slim_coords.list','w') as nf:
        nf.write('# Msa nucleotide positions and ref_nucl in reference: '+refgenome_path+'\n')
        nf.write(','.join(map(str,samples_pseudoseqs_poses))+'\n')
    #/
###/
### Make real seqs per sample of stable positions
samples_seqs = {}
for pos,ref_nucl in enumerate(ref_seq):
    for sample in samples_found:
        # use ref_nucl as default
        nucl = ref_nucl
        #/
        # Update nucleotide if we have it in "pos_fixed_nucls_samples_vals"
        if pos in pos_fixed_nucls_samples_vals:
            nucl = pos_fixed_nucls_samples_vals[pos][sample]
        #/
        # Update nucleotide if we have it in "skipped_positions" (unstable nucls)
        if pos in skipped_positions:
            nucl = 'N'
        #/
        # check if discard "N" positions (do not output)
        if msa_no_N_out:
            continue
        #/
        # save
        init(samples_seqs,sample,[])
        samples_seqs[sample].append(nucl)
        #/
###/
### Write full seqs
if 1 and 'write full seqs?':
    with open(output_dir+'/'+'msa_full.fasta','w') as nf:
        for sample,seq in samples_seqs.items():
            seq = ''.join(seq)
            nf.write('>'+sample+'\n'+seq+'\n')
###/

### Write rescued "missing data" positions
with open(output_dir+'/pos_lowcov_dataset_rescues.tsv','w') as nf:
    # write header
    header = ['position','all_datasets_same_nucl','datasets']
    nf.write('\t'.join(map(str,header))+'\n')
    #/
    # write rows
    for pos,datasets in sorted(pos_lowcov_dataset_rescues.items(),key=lambda x: x[0]):
        # check if this pos was discarded due to all datasets having the same nucl in the MSA
        all_datasets_same_nucl = False
        if pos in pos_fixed_nucls_samples_vals_allTheSame:
            all_datasets_same_nucl = True
        #/
        writeArr= [pos,all_datasets_same_nucl,','.join(datasets)]
        nf.write('\t'.join(map(str,writeArr))+'\n')
    #/
###/
######/