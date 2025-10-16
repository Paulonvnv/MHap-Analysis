#!/usr/bin/env python
"""Align ASVs to target amplicon reference and report variants as CIGAR strings
"""
import argparse
import sys
import os
import subprocess
import re

from Bio import SeqIO, AlignIO

# default location of amplicon sequences and DUST mask info for the gt-seq panel
AMPLICON_DATABASE="/gsap/garage-protistvector/ampseq_data/AmpSeQC/amplicons_noprimers.fasta"

verbose = False # set to true to report more messages

# parse amplicon dust mask info
def parse_dustmasker(mask_info):
    """Parse DUST accloc format mask info"""
    if not mask_info:
        return
    mask = {}
    with open(mask_info) as f:
        for line in f:
            line = line.strip().split("\t")
            gene = line[0].split(":")[0][1:]
            if gene not in mask:
                mask[gene] = set()
            start = int(line[1])+1 # mask info is 0-based, but we want 1-based
            end = int(line[2])+2 # +1 for 1-based and +1 to include last pos in range
            mask[gene].update(list(range(start, end))) # add all pos in between start and end
    if not mask:
        print("ERROR: No mask data loaded! Is the file the correct format?", file=sys.stderr)
        sys.exit(1)
    return mask


# parse amplicon database
def parse_amp_db(fasta_file=AMPLICON_DATABASE):
    """Load sequences from fasta file of amplicons"""
    amplicons = {}
    for seq in SeqIO.parse(fasta_file, "fasta"):
        amplicons[seq.id] = seq
    return amplicons


# parse asv to amplicon table
def parse_asv_table(file, min_reads=0, min_samples=0, max_snv_dist=-1, max_indel_dist=-1, include_failed=False, exclude_bimeras=False):
    """Parse DADA2 ASV table format"""
    bins = {}
    with open(file) as f:
        f.readline()
        for line in f:
            line = line.strip().split("\t")
            # total reads
            if int(line[1]) < min_reads: 
                continue # skip if too few total reads
            # total samples
            if int(line[2]) < min_samples: 
                continue # skip if in too few samples
            # minimum SNV distance
            if max_snv_dist >= 0 and int(line[5]) > max_snv_dist:
                continue # skip if snv distance > threshold
            # minimum indel distance
            if max_indel_dist >= 0 and int(line[6]) > max_indel_dist:
                continue # skip if indel distance > threshold
            # check for failing the snv_filter and indel_filter
            if not include_failed and (line[-3] == "FAIL" or line[-2] == "FAIL"):
                continue # failed post-DADA2 filters
            # check for dada2 calling asv a bimera
            if exclude_bimeras and line[-1] == "TRUE":
                continue # skip if dada2 called bimera
            ASV = line[0] # (e.g. ASV123)
            amplicon = line[4] # target gene/amplicon
            if amplicon not in bins:
                bins[amplicon] = []
            bins[amplicon].append(ASV)
    return bins # bins is dict of amplicon -> list of ASVs assigned to amplicon


# parse ASV fasta file
def get_asv_seqs(file):
    """Load ASV sequences from fasta file"""
    return {seq.id: seq for seq in SeqIO.parse(file, "fasta")}


# write amplicon fasta files
def wrte_amplicon_fastas(asvs, bins, amplicons, outdir="ASVs"):
    """Write one fasta file per amplicon, containing reference sequence and assigned ASVs"""
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    for amplicon in bins:
        if amplicon not in amplicons:
            print(f"WARNING: {amplicon} target not found in amplicon sequence database", file=sys.stderr)
            continue
        with open(os.path.join(outdir, f"{amplicon}.fasta"), "w") as w:
            SeqIO.write(amplicons[amplicon], w, "fasta")
            SeqIO.write([asvs[asv] for asv in bins[amplicon]], w, "fasta")


# run muscle for each amplicon
def run_muscle(bins, outdir="ASVs"):
    """Iterate through amplicons, aligning each one with MUSCLE"""
    for amplicon in bins:
        fasta = os.path.join(outdir, f"{amplicon}.fasta")
        if not os.path.isfile(fasta):
            print(f"ERROR: Could not find {fasta}", file=sys.stderr)
            continue
        msa = os.path.join(outdir, f"{amplicon}.msa")
        subprocess.run(["muscle", "-in", fasta, "-out", msa, "-cluster", "neighborjoining"], capture_output=True) # CHANGE
        # With the argument -cluster we can apply different clustering algorithms for the progressive step,
        # Here I am using neighborjoining which seems to be the most suitable for placing indels correctly (before upgmb was used as default, other possible option is upgma)
        # Also we can use arguments like -gapopen (default -400) or -gapextend (default 0) to improve the alignment.
        # It would be good to allow flexibility in the usage of these parameters as they might change between different alignments.

# get coords of homopolymer runs
def _get_homopolymer_runs(seq, min_length=5):
    """Detect and report homopolymer runs of minimum length"""
    runs = set()
    prev = None
    run = 1
    start = None
    last_non_gap = None
    for i in range(len(seq)):
        if seq[i] == "-":
            continue
        if seq[i] == prev:
            if not start:
                if i > 1 and seq[i-2] == '-':
                    # gap at start of run
                    j = i - 2
                    while j >= 0:
                        if seq[j] != "-":
                            start = j+1 # start is the start of the gap
                            break
                        j -= 1
                    else:
                        start = 0
                else:
                    start = last_non_gap
            run += 1
        else:
            if run >= min_length:
                runs.update(list(range(start, i)))
            run = 1
            start = None
        prev = seq[i]
        last_non_gap = i
    
    return runs

# parse muscle alignment
def parse_alignment(alignment, mask={}, min_homopolymer_length=5, amplicon=None):
    """
	Parse amplicon alignment file, converting ASV to CIGAR string
	
	Description:
	This function parses an amplicon alignment file and converts ASV (Amplicon Sequence Variant) alignments to CIGAR (Compact Idiosyncratic Gapped Alignment Report) strings. It processes the alignment, identifies indels, and handles masking of positions. The function returns a dictionary mapping ASV IDs to their respective CIGAR strings.

	Parameters:

	alignment (str): The path to the amplicon alignment file.
	mask (dict, optional): A dictionary containing masking information for specific amplicon sequences. Default is an empty dictionary ({}).
	min_homopolymer_length (int, optional): The minimum length of a homopolymer run to be reported. Default is 5.
	amplicon (str, optional): The ID of the amplicon reference sequence. Default is None.
	Returns:

	asv_to_cigar (dict): A dictionary mapping ASV IDs to their corresponding CIGAR strings.
	"""
    aln = AlignIO.read(alignment, "fasta")
    # sort such that amplicon reference is first in alignment
    aln.sort(key = lambda record: (record.id != amplicon, record.id))
    anchor = aln[0]
    if anchor.id != amplicon:
        print(f"ERROR: No anchor gene for {alignment}", file=sys.stderr)
        # don't parse if amplicon reference not in alignment (this shouldn't happen)
        return

    if min_homopolymer_length > 1:
        # detect homopolymer runs in reference sequence
        homopolymer_runs = _get_homopolymer_runs(aln[0], min_length=min_homopolymer_length)

    if len(anchor.seq.lstrip("-")) != aln.get_alignment_length():
        print(f"WARNING: {os.path.basename(alignment)} extends beyond 5' end of reference sequence!", file=sys.stderr)
    elif len(anchor.seq.rstrip("-")) != aln.get_alignment_length():
        print(f"WARNING: {os.path.basename(alignment)} extends beyond 3' end of reference sequence!", file=sys.stderr)

    masked = mask.get(aln[0].id, None)
    type_of_indel = "NULL"
    asv_to_cigar = {}
    for seq in aln[1:]: # for each sequence but the reference or anchor
        #print(seq.name)
        
        if len(anchor.seq.lstrip("-")) < len(anchor):
            
            pos = 1 # start at position 1 in anchor (reference) sequence
            cigar = {} #""  # cigar string to dictionary output, start empty
            indel = False # indicate alignment column in an indel
            masking = False # indicate alignment column is being masked
            for i in range(aln.get_alignment_length()): # For each position in the alignment
                
                if len(cigar) < i + 1:
                    cigar[i] = '' # dictionary slot start empty
                    
                # if anchor pos masked, or next base in anchor is masked and anchor position is a gap
                if masked and (pos in masked or (pos+1 in masked and anchor[i] == '-')):
                    if verbose and seq.id == aln[1].id:
                        if not masking:
                            print(f"INFO: Skipping masked positions starting at {pos} in {os.path.basename(alignment)}", file=sys.stderr)
                            if anchor[i] == '-':
                                print(f"INFO: Gap in alignment at start of masked region!", file=sys.stderr)
                            masking = True
                        elif pos not in masked:
                            print(f"INFO: Ending masked positions at {pos-1} in {os.path.basename(alignment)}", file=sys.stderr)
                            masking = False
                elif min_homopolymer_length > 1 and i in homopolymer_runs:
                    if verbose and seq.id == aln[1].id:
                        if i and i-1 not in homopolymer_runs:
                            print(f"INFO: Skipping homopolymer run (poly-{anchor[i]}) beginning at position {pos} in {os.path.basename(alignment)}", file=sys.stderr)
                        elif i+1 not in homopolymer_runs:
                            print(f"INFO: End of homopolymer run (poly-{anchor[i]}) at position {pos} in {os.path.basename(alignment)}", file=sys.stderr)
                elif seq[i] != anchor[i]: # If there is a mutation
                    if anchor[i] == "-" or seq[i] == "-": # If mutation is an INDEL
                        if not indel: # if previous position was not an indel
                            indel = True
                        
                            # DEFINE TYPE OF INDEL
                            if anchor[i] == "-" and seq[i] != "-": # A possible insertion (it could be a an insertion or and insertion plus a deletion)
                                
                                anchor_indel_ends = i

                                # Define the end of the gap in the anchor
                                if anchor_indel_ends + 1 < len(anchor.seq.rstrip("-")):
                                    while anchor[anchor_indel_ends + 1] == "-" and anchor_indel_ends < len(anchor.seq.rstrip("-")): 
                                        if anchor_indel_ends + 1 < len(anchor.seq.rstrip("-")):
                                            anchor_indel_ends += 1
                                        if anchor_indel_ends == len(anchor) - 1:
                                            break
                            
                                    if (anchor[anchor_indel_ends + 1] != "-" and seq[anchor_indel_ends + 1] == "-"): # if the positon in the sequence of interest after the end in the anchor is a gap
                                        type_of_indel = "indel"
                                        
                                    elif seq[anchor_indel_ends + 1] != "-": # if the positon in the sequence of interest after the end in the anchor is not a gap
                                        type_of_indel = "insertion"
                                else: 
                                    type_of_indel = "insertion"                        
                                    
                            elif anchor[i] != "-" and seq[i] == "-": # A possible deletion
                                
                                
                                if i == len(seq.seq.rstrip('-')):
                                    
                                    type_of_indel = "deletion"
                                    
                                else:
                                    seq_indel_ends = i

                                    # Define the end of the gap in the sequence of interest
                                    if seq_indel_ends + 1 < len(anchor):
                                        while seq[seq_indel_ends + 1] == "-" and seq_indel_ends + 1 < len(anchor):
                                            if seq_indel_ends + 1 < len(anchor):
                                                seq_indel_ends += 1
                                            if seq_indel_ends == len(anchor) - 1:
                                                break
                                        if seq_indel_ends == len(anchor) - 1:
                                            if (anchor[seq_indel_ends] == "-" and seq[seq_indel_ends] != "-"): # if the positon in the anchor after the end in the sequence is a gap
                                                type_of_indel = "indel"
                                            else:
                                                type_of_indel = "deletion"
                                        elif (anchor[seq_indel_ends + 1] == "-" and seq[seq_indel_ends + 1] != "-"): # if the positon in the anchor after the end in the sequence is a gap
                                            type_of_indel = "indel"
                                        elif ((anchor[seq_indel_ends + 1] != "-" and seq[seq_indel_ends + 1] != "-") or # if the positon in the anchor after the end in the sequence is not a gap
                                            (anchor[seq_indel_ends + 1] == "-" and seq[seq_indel_ends + 1] == "-" and seq_indel_ends == len(anchor.seq.rstrip("-")))): # or is not the end of the alignment
                                            type_of_indel = "deletion"
                                    else:
                                        type_of_indel = "deletion"
                                        
                            
                            # START WRITING the CIGAR string for INDELs
                            if type_of_indel == "insertion": # If the indel is an insertion
                                cigar[pos - 1] = cigar[pos - 1] + f"{pos - 1}I=" # We are referencing the position BEFORE the INSERTION starts, ei. if insertion starts at position 13, we will write position 12, the nucleotide present in the sequence of interest in position 12, and then all nucleotides the are part of the insertion
                                if pos - 1 > 0:
                                    
                                    j = 1
                                    while seq[i-j] == "-":
                                        j += 1
                                    cigar[pos - 1] = cigar[pos - 1] + f"{seq[i - j]}"
                                    
                                    #for j in range(1,len(anchor)-i + 1): # localize the nucleotide in the sample before INSERTION starts
                                    #    if anchor[i-j] != "-":
                                    #        cigar[pos - 1] = cigar[pos - 1] + f"{seq[i - j]}" # Write the nucleotide before insertion starts to the CIGAR string
                                    #        break
                                        
                                if seq[i] != "-":
                                    cigar[pos - 1] = cigar[pos - 1] + f"{seq[i]}" # Write the first inserted nucleotide
                                
                            elif type_of_indel == "deletion": # if the indel is a deletion
                                cigar[pos] = cigar[pos] + f"{pos}D=" # In contrast to insertions, here we are going to report the position where the DELETION starts
                                if anchor[i] != "-":
                                    cigar[pos] = cigar[pos] + f"{anchor[i]}" # Write the first deleted nucleuotide
                        
                            elif type_of_indel == "indel": # if the indel includes both an insertion and a deletion
                            
                                insertion_pos = pos - 1 
                                cigar[insertion_pos] = cigar[insertion_pos] + f"{insertion_pos}I=" # Write the position before insertion starts
                                if insertion_pos > 0:
                                    
                                    j = 1
                                    while seq[i-j] == "-":
                                        j += 1
                                    cigar[pos - 1] = cigar[pos - 1] + f"{seq[i - j]}"
                                    
                                    #for j in range(1,len(anchor)-i + 1): # localize the nucleotide in the sample before INSERTION starts
                                    #    if anchor[i-j] != "-":
                                    #        cigar[insertion_pos] = cigar[insertion_pos] + f"{seq[i - j]}" # Write the nucleotide before insertion starts to the CIGAR string
                                    #        break
                                if seq[i] != "-":
                                    cigar[insertion_pos] = cigar[insertion_pos] + f"{seq[i]}" # Write the first inserted nucleotide
                                
                                deletion_pos = pos
                                if deletion_pos >= 1 and i == 0:
                                    print("INDEL happens before the start of the reference sequence")
                                    cigar[deletion_pos] = f"{deletion_pos}D=" # Report the position where the DELETION starts

                                else:
                                    cigar[deletion_pos] = cigar[deletion_pos] + f"{deletion_pos}D=" # Report the position where the DELETION starts
                                    if anchor[i] != "-":
                                        cigar[deletion_pos] = cigar[deletion_pos] + f"{anchor[i]}" # Write the first deleted nucleuotide
                    
                        # CONTINUE WRITING CIGAR string for indels     
                        elif type_of_indel == 'insertion':
                            if seq[i] != "-" and anchor[i] == '-':
                                cigar[pos - 1] = cigar[pos - 1] + f"{seq[i]}" # Write the following inserted nucleotides
                        elif type_of_indel == "deletion":
                            if anchor[i] != "-":
                                cigar[pos] = cigar[pos] + f"{anchor[i]}" # Write the following deleted nucleotides
                            
                        elif type_of_indel == 'indel':
                            if (seq[i] != "-" and anchor[i] == "-"):
                                cigar[insertion_pos] = cigar[insertion_pos] + f"{seq[i]}" # Write the following inserted nucleotides
                            if anchor[i] != "-" and seq[i] == "-":
                                cigar[deletion_pos] = cigar[deletion_pos] + f"{anchor[i]}" # Write the following deleted nucleotides

                    else: # If mutation is a SUBSTITUTION
                        
                        if pos <= aln.get_alignment_length() - str(anchor.seq).count('-') and i <= len(anchor.seq.rstrip("-")): # Count substitutions up to the end of the length of the reference, for mutations in insertions in the 3' flanking region they are going to be reported in the insertion
                            
                            if i == len(anchor) - 1: # if it is the last position
                                cigar[pos] = cigar[pos] + f"{pos}{seq[i]}" # Write substitution in the current position
                            elif anchor[i + 1] != "-" and seq[i + 1] != "-": # if there are no INDELs in the next position
                                cigar[pos] = cigar[pos] + f"{pos}{seq[i]}" # Write substitution in the current position
                            elif anchor[i + 1] == "-" and seq[i + 1] == "-": # If there are INDELS in the next position of both the Reference and the Sequence of interest
                                # Measure where the deletion ends for the reference and the sequence of interest
                                anchor_indel_ends = i
                                seq_indel_ends = i
                                
                                if anchor_indel_ends < len(anchor) - 1:
                                    while (anchor[anchor_indel_ends + 1] == "-" or seq[seq_indel_ends + 1] == "-") and (anchor_indel_ends < len(anchor) - 1):
                                        
                                        if anchor[anchor_indel_ends + 1] == "-":
                                            anchor_indel_ends += 1
                                        
                                        if seq[seq_indel_ends + 1] == "-":
                                            seq_indel_ends += 1
                                        if seq_indel_ends == len(anchor) - 1 or anchor_indel_ends == len(anchor) - 1:
                                            break
                                            
                                    if anchor_indel_ends <= seq_indel_ends: # if the deletion length is the same or greater than the reference
                                        cigar[pos] = cigar[pos] + f"{pos}{seq[i]}"
                                elif anchor_indel_ends == len(anchor) - 1:
                                    cigar[pos] = cigar[pos] + f"{pos}{seq[i]}"
                            elif anchor[i + 1] != "-" and seq[i + 1] == "-": # if there is a gap only in the sequence of interest
                                seq_indel_ends = i + 1
                                if seq_indel_ends < len(anchor) - 1:
                                    while seq[seq_indel_ends + 1] == "-" and seq_indel_ends + 1 < len(anchor.seq.rstrip("-")) :
                                        if seq_indel_ends + 1 < len(anchor.seq.rstrip("-")):
                                            seq_indel_ends += 1
                                        if seq_indel_ends == len(anchor) - 1:
                                            break
                                    if seq_indel_ends == len(anchor.seq.rstrip("-")) - 1:
                                        if (anchor[seq_indel_ends] == "-" and seq[seq_indel_ends] != "-"):
                                            print("Substitution will be written within the inserttion")
                                        else:
                                            cigar[pos] = cigar[pos] + f"{pos}{seq[i]}"
                                    elif (anchor[seq_indel_ends + 1] == "-" and seq[seq_indel_ends + 1] != "-"):
                                        print("Substitution will be written within the inserttion")
                                    elif ((anchor[seq_indel_ends + 1] != "-" and seq[seq_indel_ends + 1] != "-") or 
                                        (anchor[seq_indel_ends + 1] == "-" and seq[seq_indel_ends + 1] == "-" and seq_indel_ends == len(anchor.seq.rstrip("-")))):
                                
                                        cigar[pos] = cigar[pos] + f"{pos}{seq[i]}"
                                
                                else:
                                    cigar[pos] = cigar[pos] + f"{pos}{seq[i]}"
                        indel = False
                elif seq[i] == "-" and anchor[i] == "-" and indel: # I the gap is present in both the anchor and the sequence, extend the indel
                    indel = True
                    
                else:
                    indel = False
                if anchor[i] != '-': # Move to next position in the anchor if there is no gap in it
                    pos += 1
            
            cigar = dict(sorted(cigar.items())) # Sort all positions in CIGAR string
            cigar = ''.join(cigar.values()) # # Join all positions in a single string
            if cigar == '':
                cigar = "."
                
            asv_to_cigar[seq.id] = cigar # Store CIGAR string in the final output
            
        elif len(anchor.seq.lstrip("-")) == len(anchor):
            
            pos = 1 # start at position 1 in anchor (reference) sequence
            cigar = {} #""  # cigar string to dictionary output, start empty
            indel = False # indicate alignment column in an indel
            masking = False # indicate alignment column is being masked
            for i in range(aln.get_alignment_length()): # For each position in the alignment
                
                if len(cigar) < i + 1:
                    cigar[i + 1] = '' # dictionary slot start empty
            
                # if anchor pos masked, or next base in anchor is masked and anchor position is a gap
                if masked and (pos in masked or (pos+1 in masked and anchor[i] == '-')):
                    if verbose and seq.id == aln[1].id:
                        if not masking:
                            print(f"INFO: Skipping masked positions starting at {pos} in {os.path.basename(alignment)}", file=sys.stderr)
                            if anchor[i] == '-':
                                print(f"INFO: Gap in alignment at start of masked region!", file=sys.stderr)
                            masking = True
                        elif pos not in masked:
                            print(f"INFO: Ending masked positions at {pos-1} in {os.path.basename(alignment)}", file=sys.stderr)
                            masking = False
                elif min_homopolymer_length > 1 and i in homopolymer_runs:
                    if verbose and seq.id == aln[1].id:
                        if i and i-1 not in homopolymer_runs:
                            print(f"INFO: Skipping homopolymer run (poly-{anchor[i]}) beginning at position {pos} in {os.path.basename(alignment)}", file=sys.stderr)
                        elif i+1 not in homopolymer_runs:
                            print(f"INFO: End of homopolymer run (poly-{anchor[i]}) at position {pos} in {os.path.basename(alignment)}", file=sys.stderr)
                elif seq[i] != anchor[i]: # If there is a mutation
                    if anchor[i] == "-" or seq[i] == "-": # If mutation is an INDEL
                        if not indel: # if previous position was not an indel
                            indel = True
                        
                            # DEFINE TYPE OF INDEL
                            if anchor[i] == "-" and seq[i] != "-": # A possible insertion (it could be a an insertion or and insertion plus a deletion)
                                anchor_indel_ends = i
                                
                                # Define the end of the gap in the anchor
                                if anchor_indel_ends + 1 < len(anchor.seq.rstrip("-")):
                                    
                                    while anchor[anchor_indel_ends + 1] == "-" and anchor_indel_ends < len(anchor.seq.rstrip("-")): 
                                        if anchor_indel_ends + 1 < len(anchor.seq.rstrip("-")):
                                            anchor_indel_ends += 1
                                        if anchor_indel_ends == len(anchor) - 1:
                                            break
                                            
                                    if (anchor[anchor_indel_ends + 1] != "-" and seq[anchor_indel_ends + 1] == "-"): # if the positon in the sequence of interest after the end in the anchor is a gap
                                        type_of_indel = "indel"
                                
                                    elif seq[anchor_indel_ends + 1] != "-": # if the positon in the sequence of interest after the end in the anchor is not a gap
                                        type_of_indel = "insertion"
                                        
                                else: 
                                    type_of_indel = "insertion"                           
                                    
                            elif anchor[i] != "-" and seq[i] == "-": # A possible deletion
                                
                                if i == len(seq.seq.rstrip('-')):
                                    
                                    type_of_indel = "deletion"
                                else:
                                    seq_indel_ends = i

                                    # Define the end of the gap in the sequence of interest
                                    if seq_indel_ends + 1 < len(anchor):
                                        
                                        while seq[seq_indel_ends + 1] == "-" and seq_indel_ends + 1 < len(anchor):
                                            if seq_indel_ends + 1 < len(anchor):
                                                seq_indel_ends += 1
                                            if seq_indel_ends == len(anchor) - 1:
                                                break
                                                
                                        if seq_indel_ends == len(anchor) - 1:
                                            if (anchor[seq_indel_ends] == "-" and seq[seq_indel_ends] != "-"): # if the positon in the anchor after the end in the sequence is a gap
                                                type_of_indel = "indel"
                                                
                                            else:
                                                type_of_indel = "deletion"
                                                
                                        elif (anchor[seq_indel_ends + 1] == "-" and seq[seq_indel_ends + 1] != "-"): # if the positon in the anchor after the end in the sequence is a gap
                                            type_of_indel = "indel"
                                            
                                        elif ((anchor[seq_indel_ends + 1] != "-" and seq[seq_indel_ends + 1] != "-") or # if the positon in the anchor after the end in the sequence is not a gap
                                            (anchor[seq_indel_ends + 1] == "-" and seq[seq_indel_ends + 1] == "-" and seq_indel_ends == len(anchor.seq.rstrip("-")))): # or is not the end of the alignment
                                            type_of_indel = "deletion"
                                            
                                    else:
                                        type_of_indel = "deletion"
                                
                        
                            # START WRITING the CIGAR string for INDELs
                            
                            if type_of_indel == "insertion": # If the indel is an insertion
                                cigar[pos - 1] = cigar[pos - 1] + f"{pos - 1}I=" # We are referencing the position BEFORE the INSERTION starts, ei. if insertion starts at position 13, we will write position 12, the nucleotide present in the sequence of interest in position 12, and then all nucleotides the are part of the insertion
                            
                                if pos - 1 > 0:
                                    
                                    j = 1
                                    while seq[i-j] == "-":
                                        j += 1
                                    cigar[pos - 1] = cigar[pos - 1] + f"{seq[i - j]}"
                                    
                                    #for j in range(1,len(anchor) - i + 1): # localize the nucleotide in the sample before INSERTION starts
                                        
                                    #    if anchor[i-j] != "-":
                                    #        cigar[pos - 1] = cigar[pos - 1] + f"{seq[i - j]}" # Write the nucleotide before insertion starts to the CIGAR string
                                    #        break
                                if seq[i] != "-":
                                    cigar[pos - 1] = cigar[pos - 1] + f"{seq[i]}" # Write the first inserted nucleotide
                                
                            elif type_of_indel == "deletion": # if the indel is a deletion
                                cigar[pos] = cigar[pos] + f"{pos}D=" # In contrast to insertions, here we are going to report the position where the DELETION starts
                                if anchor[i] != "-":
                                    cigar[pos] = cigar[pos] + f"{anchor[i]}" # Write the first deleted nucleuotide
                        
                            elif type_of_indel == "indel": # if the indel includes both an insertion and a deletion
                                
                                insertion_pos = pos - 1
                                if insertion_pos == 0:
                                    cigar[insertion_pos] = f"{insertion_pos}I=" # Write the position before insertion starts
                                else:
                                    cigar[insertion_pos] = cigar[insertion_pos] + f"{insertion_pos}I=" # Write the position before insertion starts
                                if insertion_pos > 0:
                                    
                                    j = 1
                                    while seq[i-j] == "-":
                                        j += 1
                                    cigar[pos - 1] = cigar[pos - 1] + f"{seq[i - j]}"
                                    
                                    #for j in range(1,len(anchor)-i + 1): # localize the nucleotide in the sample before INSERTION starts
                                    #    if anchor[i-j] != "-":
                                    #        cigar[insertion_pos] = cigar[insertion_pos] + f"{seq[i - j]}" # Write the nucleotide before insertion starts to the CIGAR string
                                    #        break
                                if seq[i] != "-":
                                    cigar[insertion_pos] = cigar[insertion_pos] + f"{seq[i]}" # Write the first inserted nucleotide
                                
                                deletion_pos = pos
                                cigar[deletion_pos] = cigar[deletion_pos] + f"{deletion_pos}D=" # Report the position where the DELETION starts
                                if anchor[i] != "-":
                                    cigar[deletion_pos] = cigar[deletion_pos] + f"{anchor[i]}" # Write the first deleted nucleuotide
                    
                        # CONTINUE WRITING CIGAR string for indels     
                        elif type_of_indel == 'insertion':
                            if seq[i] != "-" and anchor[i] == '-':
                                cigar[pos - 1] = cigar[pos - 1] + f"{seq[i]}" # Write the following inserted nucleotides
                        elif type_of_indel == "deletion":
                            if anchor[i] != "-":
                                cigar[pos] = cigar[pos] + f"{anchor[i]}" # Write the following deleted nucleotides
                            
                        elif type_of_indel == 'indel':
                            if (seq[i] != "-" and anchor[i] == "-"):
                                cigar[insertion_pos] = cigar[insertion_pos] + f"{seq[i]}" # Write the following inserted nucleotides
                            if (anchor[i] != "-" and seq[i] == "-"):
                                cigar[deletion_pos] = cigar[deletion_pos] + f"{anchor[i]}" # Write the following deleted nucleotides

                    else: # If mutation is a SUBSTITUTION
                        
                        if pos <= aln.get_alignment_length() - str(anchor.seq).count('-') and i <= len(anchor.seq.rstrip("-")): # Count substitutions up to the end of the length of the reference, for mutations in insertions in the 3' flanking region they are going to be reported in the insertion
                            
                            if i == len(anchor) - 1: # if it is the last position
                                cigar[pos] = cigar[pos] + f"{pos}{seq[i]}" # Write substitution in the current position
                            elif anchor[i + 1] != "-" and seq[i + 1] != "-": # if there are no INDELs in the next position
                                cigar[pos] = cigar[pos] + f"{pos}{seq[i]}" # Write substitution in the current position
                            elif anchor[i + 1] == "-" and seq[i + 1] == "-": # If there are INDELS in the next position of both the Reference and the Sequence of interest
                                # Measure where the deletion ends for the reference and the sequence of interest
                                
                                anchor_indel_ends = i
                                seq_indel_ends = i
                                if anchor_indel_ends < len(anchor) - 1:
                                    
                                    while (anchor[anchor_indel_ends + 1] == "-" or seq[seq_indel_ends + 1] == "-") and (anchor_indel_ends < len(anchor) - 1):
                                        
                                        if anchor[anchor_indel_ends + 1] == "-":
                                            anchor_indel_ends += 1
                                        
                                        if seq[seq_indel_ends + 1] == "-":
                                            seq_indel_ends += 1
                                            
                                        if seq_indel_ends == len(anchor) - 1 or anchor_indel_ends == len(anchor) - 1:
                                            break
                                    
                                            
                                    if anchor_indel_ends <= seq_indel_ends: # if the deletion length is the same or greater than the reference
                                        cigar[pos] = cigar[pos] + f"{pos}{seq[i]}"
                                elif anchor_indel_ends == len(anchor) - 1:
                                    cigar[pos] = cigar[pos] + f"{pos}{seq[i]}"
                            elif anchor[i + 1] != "-" and seq[i + 1] == "-": # if there is a gap only in the sequence of interest
                                seq_indel_ends = i + 1
                                
            
                                if seq_indel_ends < len(anchor) - 1:
                                    while seq[seq_indel_ends + 1] == "-" and seq_indel_ends + 1 < len(anchor.seq.rstrip("-")) :
                                        if seq_indel_ends + 1 < len(anchor.seq.rstrip("-")):
                                            seq_indel_ends += 1
                                        if seq_indel_ends == len(anchor) - 1:
                                            break
                                            
                                    if seq_indel_ends == len(anchor.seq.rstrip("-")) - 1:
                                        if (anchor[seq_indel_ends] == "-" and seq[seq_indel_ends] != "-"):
                                            print("Substitution will be written within the inserttion")
                                        else:
                                            cigar[pos] = cigar[pos] + f"{pos}{seq[i]}"
                                    elif (anchor[seq_indel_ends + 1] == "-" and seq[seq_indel_ends + 1] != "-"):
                                        print("Substitution will be written within the inserttion")
                                    elif ((anchor[seq_indel_ends + 1] != "-" and seq[seq_indel_ends + 1] != "-") or 
                                        (anchor[seq_indel_ends + 1] == "-" and seq[seq_indel_ends + 1] == "-" and seq_indel_ends == len(anchor.seq.rstrip("-")))):
                                
                                        cigar[pos] = cigar[pos] + f"{pos}{seq[i]}"
                                
                                else:
                                    cigar[pos] = cigar[pos] + f"{pos}{seq[i]}"
                                
                        indel = False
                elif seq[i] == "-" and anchor[i] == "-" and indel: # I the gap is present in both the anchor and the sequence, extend the indel
                    indel = True
                else:
                    indel = False
                if anchor[i] != '-': # Move to next position in the anchor if there is no gap in it
                    pos += 1
        
            cigar = dict(sorted(cigar.items())) # Sort all positions in CIGAR string
            cigar = ''.join(cigar.values()) # # Join all positions in a single string
        
            if cigar == '':
                cigar = "."
            
            asv_to_cigar[seq.id] = cigar # Store CIGAR string in the final output
        
    return asv_to_cigar

# get variants per amplicon per position
def parse_alignments(bins, mask={}, min_homopolymer_length=5, outdir="ASVs"):
    """Parse multi-sequence alignment fasta file from MUSCLE"""
    cigars = {}
    for amplicon in sorted(bins):
        print(amplicon)
        msa = os.path.join(outdir, f"{amplicon}.msa")
        if not os.path.isfile(msa):
            print(f"ERROR: Could not find {msa}", file=sys.stderr)
            continue
        # store CIGAR strings per amplicon in dict
        cigars[amplicon] = parse_alignment(msa, mask=mask, min_homopolymer_length=min_homopolymer_length, amplicon=amplicon)
    
    return cigars


# write table of asv -> amplicon/cigar
def write_cigar_strings(cigars, out):
    """Write conversion table from ASV to CIGAR string"""
    number = re.compile(r"\d+")
    with open(out, 'w') as w:
        # write tab file with ASV, amplicon target, and CIGAR string
        w.write("ASV\tAmplicon\tCIGAR\n")
        for amplicon in sorted(cigars):
            # sort on ASV number
            for ASV in sorted(cigars[amplicon], key = lambda x: int(number.search(x).group())):
                w.write(f"{ASV}\t{amplicon}\t{cigars[amplicon][ASV]}\n")


def convert_seqtab(file, cigars, out):
    """Parse seqtab file, converting ASVs to CIGAR strings"""
    # get dict of ASVs -> amplicon/CIGAR
    asv_to_cigar = {}
    variants = set()
    for amplicon in sorted(cigars):
        for ASV in sorted(cigars[amplicon]):
            variant = f"{amplicon},{cigars[amplicon][ASV]}"
            asv_to_cigar[ASV] = variant
            variants.add(variant)
    
    if not variants:
        print("ERROR: No haplotypes to convert!", file=sys.stderr)
        return
    
    total_reads = {}
    # parse seqtab file
    with open(file) as f:
        seqtab = {}
        f.readline()
        for line in f:
            line = line.strip().split("\t")
            sample = line[0]
            seqtab[sample] = {}
            # iterate through each ASV (i.e. ASV1, ASV2, ... ASVn)
            for i, count in enumerate(line[1:]):
                asv = f"ASV{i+1}" # don't use actual sequence
                variant = asv_to_cigar.get(asv)
                if not variant:
                    continue # ASV was filtered out
                # sum ASVs per sample that are the same variant
                count = int(count)
                if variant not in total_reads:
                    total_reads[variant] = 0
                total_reads[variant] += count
                if variant not in seqtab[sample]:
                    seqtab[sample][variant] = 0
                seqtab[sample][variant] += count
        
        if not seqtab:
            print("ERROR: No seqtab data to write!", file=sys.stderr)
            return

        # write output file (sort variants first)
        variants = sorted(list(variants), key=lambda variant: total_reads.get(variant, 0), reverse=True)
        with open(out, "w") as w:
            # write header
            w.write("sample\t" + "\t".join(variants) + "\n")
            # write one sample per line
            for sample in sorted(seqtab):
                w.write(f"{sample}\t" + "\t".join([f"{seqtab[sample].get(variant, 0)}" for variant in variants]) + "\n")
        
            return True


parser = argparse.ArgumentParser(usage="%(prog)s [options] fasta table alignments out",
                                 description="Convert ASVs from DADA2 pipeline to pseudo-CIGAR strings.",
                                 epilog="(C)2021 Broad Institute")
parser.add_argument("fasta", help="Fasta file of ASV sequences from DADA2 pipeline")
parser.add_argument("table", help="ASV table from DADA2 pipeline")
parser.add_argument("seqtab", help="DADA2 seqtab tsv file")
parser.add_argument("out", help="Output seqtab tsv file with amplicon/variant counts")
parser.add_argument("--asv_to_cigar", help="Output file for ASV -> CIGAR string table")
parser.add_argument("-a", "--alignments", default="alignments", help="Directory to store ASV alignment files (default: alignments)")
parser.add_argument("-p", "--polyN", type=int, default=5, help="Mask homopolymer runs length >= polyN (default: 5; disabled < 2)")
parser.add_argument("-r", "--min_reads", type=int, default=0, help="Minimum total reads to include ASV (default: 0, disabled)")
parser.add_argument("-n", "--min_samples", type=int, default=0, help="Minimum samples to include ASV (default: 0, disabled)")
parser.add_argument("-f", "--include_failed", action="store_true", default=False, help="INCLUDE ASVs that failed post-DADA2 filters (default: False)")
parser.add_argument("-b", "--exclude_bimeras", action="store_true", default=False, help="EXCLUDE ASVs that DADA2 flagged as bimeras (default: False)")
parser.add_argument("-s", "--max_snv_dist", type=int, default=-1, help="Maximum SNV distance to include ASV (default: -1, disabled)")
parser.add_argument("-i", "--max_indel_dist", type=int, default=-1, help="Maximum indel distance to include ASV (default: -1, disabled)")
parser.add_argument("-d", "--amp_db", default=AMPLICON_DATABASE, help=f"Amplicon sequence fasta file (default: {AMPLICON_DATABASE})")
parser.add_argument("-m", "--amp_mask", default=None, help=f"Amplicon low complexity mask info (default: None, disabled)")
parser.add_argument("-v", "--verbose", default=False, action='store_true', help="Increase verbosity")
args = parser.parse_args()

if args.verbose:
    verbose = True

print(f"INFO: Loading {args.amp_db}", file=sys.stderr)
amplicons = parse_amp_db(args.amp_db)
if not amplicons:
    print(f"ERROR: No amplicons in {args.amp_db}", file=sys.stderr)
    sys.exit(1)


if args.amp_mask:
    print(f"INFO: Loading {args.amp_mask}", file=sys.stderr)
    mask = parse_dustmasker(args.amp_mask)
else:
    print(f"INFO: No mask data specified.", file=sys.stderr)
    mask = {}

print(f"INFO: Loading {args.fasta}")
asvs = get_asv_seqs(args.fasta)
if not asvs:
    print(f"ERROR: No ASV sequences in {args.fasta}", file=sys.stderr)
    sys.exit(1)

print(f"INFO: Parsing {args.table} with total reads >= {args.min_reads}, samples >= {args.min_samples}, snv_dist <= {args.max_snv_dist}, indel_dist <= {args.max_indel_dist}", file=sys.stderr)

if args.include_failed:
    print("WARNING: Including ASVs that failed post-DADA2 filters! This is not recommended.", file=sys.stderr)
else:
    print("INFO: Excluding ASVs that failed post-DADA2 filters.", file=sys.stderr)

if args.exclude_bimeras:
    print("INFO: Excluding ASVs that DADA2 marked as bimeras.", file=sys.stderr)

bins = parse_asv_table(args.table, min_reads=args.min_reads, min_samples=args.min_samples, max_snv_dist=args.max_snv_dist, max_indel_dist=args.max_indel_dist, include_failed=args.include_failed, exclude_bimeras=args.exclude_bimeras)
if not bins:
    print(f"ERROR: No useable data in {args.table}", file=sys.stderr)
    sys.exit(1)

outdir = args.alignments
print(f"INFO: Writing amplicon fasta files to {outdir}", file=sys.stderr)
if not os.path.isdir(outdir):
    os.mkdir(outdir)
wrte_amplicon_fastas(asvs, bins, amplicons, outdir=outdir)

print("INFO: Running MUSCLE aligner on amplicon fasta files. Please wait...", file=sys.stderr)
run_muscle(bins, outdir=outdir)

print("INFO: Parsing alignments to CIGAR strings", file=sys.stderr)
cigars = parse_alignments(bins, mask=mask, min_homopolymer_length=args.polyN, outdir=outdir)
if not cigars:
    print("ERROR: could not determine CIGAR strings", file=sys.stderr)
    sys.exit(1)

if args.asv_to_cigar:
    write_cigar_strings(cigars, args.asv_to_cigar)
    print(f"INFO: Wrote ASV->CIGAR table to {args.asv_to_cigar}", file=sys.stderr)


print(f"INFO: Converting DADA2 seqtab file {args.seqtab} to {args.out}", file=sys.stderr)
if convert_seqtab(args.seqtab, cigars, args.out):
    print("INFO: Completed successfully!", file=sys.stderr)

