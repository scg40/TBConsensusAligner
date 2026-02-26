import argparse
import os
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from pprint import pprint
import gzip
'''

This script creates a SNP alignment in multi fasta format from a collection of consensus genomes
one of which should be the outgroup if intended to be used later for phylogenetic analyses.

The user can choose, how many undefined states ('-' or 'N') for a polymorphic position are accepted with
the -g argument. By default, if a position has more than 90% undefined states, it will be excluded (-g 0.9).
For a polymorphic site with Ns or - or both to be included in the final SNP alignment, the site has to be poly-
morphic in terms of at least 2 different bases, plus Ns or -s or both. 
A site with just one base and Ns and or -s is not considered polymorphic.

'''
  

def flatten(t):
    return [item for sublist in t for item in sublist]


# Function to open both zipped and unzipped vcf files
def open_vcf(path):
    path = Path(path)
    # Check the first two bytes for the gzip magic number
    with open(path, 'rb') as f:
        is_gzip = f.read(2) == b'\x1f\x8b'
    
    if is_gzip:
        return gzip.open(path, "rt")
    else:
        return open(path, "r")

# Create a dictionary of fasta sequences where the sample name is the key and the value the string of bases
def fastas2dict(args):
    sequences_dict = {}
    consensus_dir = Path(args.output_dir)
    consensus_files = sorted(consensus_dir.glob("*.fasta"))

    for fasta in consensus_files:
        with fasta.open() as fh:
            header = None
            seq_parts = []
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    header = line[1:]
                else:
                    seq_parts.append(line)
            
            if header:
                # Join only once at the very end for this sample
                sequences_dict[header] = "".join(seq_parts)
                print(f"Loaded {header}: {len(sequences_dict[header])} bases")
    
    return sequences_dict




# Check if all sequences (values, string of characters) in a dictionary are the same length, return True if so
def check_length(sequences_dict):

    sequence_lengths = [len(seq) for seq in sequences_dict.values()]

    return all(x == sequence_lengths[0] for x in sequence_lengths)


# Create a 2D array from dictionary of sequences with sequence names at pos 0 and first genomic position at python pos 1
def dict2array(sequences_dict):


    print("DEBUG dict2array:")
    print("  Number of sequences:", len(sequences_dict))
    print("  Unique sequence lengths:", set(len(seq) for seq in sequences_dict.values()))

    if check_length(sequences_dict) == True:
        
        array = [[name] + list(seq) for name, seq in sequences_dict.items()]
        
        

        return array

    else: 
        print("Not all sequences are of equal length")
        return None

# Get positions from outgroup VCF corresponding to genomic positions of polymorphic sites found in alignment
def get_outgroup_vcf_var_pos(args, var_positions):
    outgroup_vcf_dict = {}
    
    # 1. Clean the outgroup name for the FASTA header
    outgroup_name = args.outgroup_name if args.outgroup_name else Path(args.outgroup_vcf).stem
    outgroup_name = outgroup_name.replace(".vcf", "").replace(".gz", "")
    
    # 2. Populate the dictionary FIRST so debug prints work
    with open_vcf(args.outgroup_vcf) as outgroup_vcf:
        for var in outgroup_vcf:
            if var.startswith('#') or not var.strip():
                continue

            # Strict tab splitting and stripping to avoid hidden characters in Galaxy
            fields = var.strip().split('\t')

            if len(fields) < 5:
                continue

            # Force POS to be a clean string
            pos = fields[1].strip()
            ref = fields[3].strip()
            alt = fields[4].strip()

            outgroup_vcf_dict[pos] = [ref, alt]

    # 3. Handle the polymorphic coordinates
    # Ensure var_positions is a list of clean strings to match the dict keys
    unique_var_positions = [str(p) for p in dict.fromkeys(var_positions)]

    # DEBUG SECTION
    print(f"DEBUG: Input positions count: {len(unique_var_positions)}")
    print(f"DEBUG: Outgroup VCF dict entries: {len(outgroup_vcf_dict)}")
    if unique_var_positions:
        print(f"DEBUG: First search key: '{unique_var_positions[0]}'")
    
    outgroup_sequence = list()
    for pos in unique_var_positions:
        # Check dictionary using the clean string key
        key = str(pos)

        if key in outgroup_vcf_dict:
            ref_val = outgroup_vcf_dict[pos][0]
            alt_val = outgroup_vcf_dict[pos][1]

            if alt_val == '.':   # No ALT allele -> use REF
                outgroup_sequence.append(ref_val)
            
            elif alt_val != '.': # ALT allele present
                if len(alt_val) > 1 or len(ref_val) > 1: # Indel
                    outgroup_sequence.append('N')
                else:   # SNP
                    outgroup_sequence.append(alt_val)
            print(f"MATCH FOUND: Position {key} is in Outgroup VCF. Base: {outgroup_vcf_dict[key]}")
        else:
            # Position not in VCF means it's the Reference base (or missing)
            # In your case, you chose "-" for missing
            print(f"MATCH FAILED: Position {key} not found in Outgroup VCF keys.")
            outgroup_sequence.append("-")

    # Add the header at the start
    result = [outgroup_name] + outgroup_sequence
    
    print(f"DEBUG: Finalll outgroup row length: {len(result)}")
    return result


# Main script to extract polymorphic positions from the alignment [aligned_sequences] 
def main(args):
    UNDEF_STAT = args.undefined_states  # Use the arg from your parser
    SEQ_DICT = fastas2dict(args)
    
    seq_names = list(SEQ_DICT.keys())
    num_seqs = len(seq_names)
    genome_length = len(SEQ_DICT[seq_names[0]])

    # Initialize pol_pos with sequence names as the first element of each row
    pol_pos = [[name] for name in seq_names]
    polymorphic_indices = []

    print(f"Analyzing {num_seqs} sequences across {genome_length} bp...")

    # Iterate through each genomic position one by one (Memory Efficient)
    for col_index in range(genome_length):
        # Build the column on the fly
        # .replace('X', 'N') handles ambiguity
        column = [SEQ_DICT[name][col_index].replace('X', 'N') for name in seq_names]
        
        # Check for polymorphism (more than 1 unique state)
        if len(set(column)) > 1:
            
            # Logic for adding column to alignment
            should_add = False
            
            if 'N' not in column and '-' not in column:
                should_add = True
            else:
                count_undef = column.count('N') + column.count('-')
                if (count_undef / num_seqs) <= UNDEF_STAT:
                    unique_states = set(column)
                    # Your specific conditions for missing data sites
                    if ('N' in column and '-' not in column) and len(unique_states) > 2:
                        should_add = True
                    elif ('-' in column and 'N' not in column) and len(unique_states) > 2:
                        should_add = True
                    elif ('N' in column and '-' in column) and len(unique_states) > 3:
                        should_add = True

            if should_add:
                for i, base in enumerate(column):
                    pol_pos[i].append(base)
                # Store the 1-based genomic index
                polymorphic_indices.append(col_index + 1)

    print("\npol_pos preview (first 5 sequences, first 10 entries):")
    for row in pol_pos[:5]:
        print(row[:10])

    # Append the outgroup sequence if provided via VCF
    if args.outgroup_vcf:
        # Use a clean, unique, and sorted list of indices
        unique_indices = sorted(list(set(polymorphic_indices)))
    
        # Get the outgroup sequence row
        outgroup_line = get_outgroup_vcf_var_pos(args, var_positions=unique_indices)

        # Safety check: Row length must match (Name + Bases)
        if len(outgroup_line) != len(pol_pos[0]):
            print(f"CRITICAL ERROR: Outgroup row length ({len(outgroup_line)}) "
                  f"mismatch with Sample row length ({len(pol_pos[0])})")
        
        pol_pos.append(outgroup_line)

    print("Number of polymorphic sites:", len(pol_pos[0]) - 1)   
    print("Number of polymorphic_indices:", len(polymorphic_indices))
   
    # Write the output file
    # Write the output file
    output_alignment_path = os.path.join(args.output_dir, 'snp_alignment.fasta')
    with open(output_alignment_path, 'w') as output_file:
        for i, row in enumerate(pol_pos):
            sequence_name = row[0]
            # Convert everything to string just in case a None or list slipped in
            sequence = "".join(str(base) for base in row[1:])

            # Standard FASTA format: >Header\nSequence\n
            output_file.write(f">{sequence_name}\n")
            output_file.write(f"{sequence}\n")
    
    print(f"Successfully wrote {len(pol_pos)} sequences to {output_alignment_path}")










        








    


