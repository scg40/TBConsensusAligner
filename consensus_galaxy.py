import argparse
import re
import os
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
import gzip
'''
This script creates consensus fastas from a collection of VCF files and the MTBC ancestor reference genome.
SNPs with a frequency >90% are encoded with the alternative base from the VCF.
SNPs with a frequency between 10% and 90% are encoded with the ambiguity base.
SNPs with a frequency <10% are encoded with the ancestral base from the reference genome.

Large deletions (coverage 0) are encoded by dashes (-).
Small deletions with frequency >90% are encoded by dashes (-).
Small deletions with frequency between 10% and 90% are encoded with a 'N'.
Small deletions with frequency <10% are encoded with the ancestral base from the reference genome.

Small insertions (alternative bases longer than reference base) are encoded with the ancestral state from the reference genome.

Sites to exclude specified in the bed files are encoded with a 'N'.
Variants that do not have the 'PASS' quality filter are encoded with a 'N'.
Sites that are not in the VCF and are covered by less than 5 reads (via depth file) are encoded with a 'N'.

'''

# Extract the basename to use in output file naming and consensus file filling
def get_basename(vcf_path):
    # Get the filename from the path (e.g., 'input_vcfs/Sample_01.recal.vcf.gz')
    filename = os.path.basename(vcf_path)
    
    # Use regex to "look ahead" for .vcf and take everything before it.
    match = re.split(r'\.vcf', filename, flags=re.IGNORECASE)
    
    return match[0]

def flatten(t):
    return [item for sublist in t for item in sublist]


# Split a string of characters into a list of the individual characters (list comprehension) useful to examine each base of f.e. ACTG individually
def split(word):
    return [char for char in word]


# Function to open either gzipped or unzipped vcf files
def open_vcf(path):

    path = Path(path)

    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    else:
        return open(path, "r")
    

# Function to get the length of the reference sequence in case of testing with a smaller reference genome
def get_reference_length(args):
    record = SeqIO.read(str(args.reference), "fasta")    # Create object with key info of the sequence
    return len(record.seq)

# Parse the BED files to get a dictionary of positions to exclude
def get_pos_to_exclude(bed_files):

    pos_to_exclude = {}
    if bed_files:

        for bed_path in bed_files:
            bed_path = Path(bed_path)

            with bed_path.open() as file:
                table = [position.strip().split('\t') for position in file] # Create a list of lists, each representing a line with the entries as items

                for i in range(1,len(table)):           # Iterate over every element of [table] -> line in BED
                    StartPosition = int(table[i][1])    
                    EndPosition = int(table[i][2])

                    ranges_of_coordinates = [i for i in range(StartPosition, EndPosition)]  # Create a list of every position of the interval of a BED line

                    for pos in ranges_of_coordinates:   # If the position of the range is not yet in the dict, add it
                        if pos not in pos_to_exclude:
                            pos_to_exclude[pos] = ''
                        else:
                            continue
            
    return dict(sorted(pos_to_exclude.items()))
   
        
# Create a dictionary where the keys are the genomic positions and values are the bases of the reference genome
def fasta2dict(args):

    fasta_dict = {}

    record = SeqIO.read(str(args.reference), "fasta")    # Create object with key info of the sequence
    number_of_refbases = len(record.seq)

    for i in range(number_of_refbases):     # Add the each refbase as value to the key (genomic position)
        fasta_dict[i+1]=record.seq[i]

    return(fasta_dict)

# Create a dictionary from the VCFs where POS is the key and REF, ALT, QUAL and AF(as list) are values
def vcf2dict(vcf_files):

    vcf_dict = {}   
    
    for vcf_path in vcf_files:
        vcf_path = Path(vcf_path)

        with open_vcf(vcf_path) as vcf_file:

            for row in vcf_file:

                if row.startswith("#"):
                    continue

                row = row.strip().split('\t')

                position = row[1]
                ref = row[3]
                alt = row[4]
                qual = float(row[5])

                # Compute the AF from AO and RO for one or multiple ALT alleles with AF = AO / AO + RO
                value_fields = row[-1].split(":") # Split column with the values at the :
                RO = int(value_fields[2])    # As int to do calculations with it
                AO_list = [int(ao) for ao in value_fields[4].split(",")]  # A list in case calculation of AF has do be done for multiple alleles

                AF = [ao / (sum(AO_list) + RO) for ao in AO_list]   # sum(AO_list) to accomodate possibility of many alleles

                vcf_dict[int(position)] = [ref,alt,qual,AF]
        
    return(vcf_dict)

# Parse the depth file to get a nested dictionary with with the depths of each genomic position for each BAM file 
# {{ BAM0: {pos1: depth1, pos2: depth2, ...}, BAM1: {pos1: depth1, pos2: depth2, ...}}
def get_sample_depth(depth_file_path, column_index):
    """
    column_index: the integer index of the sample (0, 1, 2...) 
    This corresponds to the order Galaxy passed the VCFs.
    """
    sample_depths = {}
    with open(depth_file_path, 'rt') as f:
        # Skip the header
        header = f.readline() 
        
        # We add 2 to column_index because:
        # Col 0 = CHROM, Col 1 = POS, Col 2 = First Sample (Index 0)
        actual_col = column_index + 2
        
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            # Pull the POS and the specific DEPTH column
            sample_depths[int(parts[1])] = int(parts[actual_col])
            
    return sample_depths

# Define the ambiguity bases to return ambiguity base if SNP between AF 10% and 90%
def ambiguity_code(ref, alt):

    ambiguity_base = ''

    if ref == 'C' and alt == 'T':
        ambiguity_base = 'Y'

    elif ref == 'T' and alt == 'C':
        ambiguity_base = 'Y'

    elif ref == 'A' and alt == 'G':
        ambiguity_base = 'R'

    elif ref == 'G' and alt == 'A':
        ambiguity_base = 'R'

    elif ref == 'A' and alt == 'T':
        ambiguity_base = 'W'

    elif ref == 'T' and alt == 'A':
        ambiguity_base = 'W'

    elif ref == 'G' and alt == 'C':
        ambiguity_base = 'S'

    elif ref == 'C' and alt == 'G':
        ambiguity_base = 'S'

    elif ref == 'T' and alt == 'G':
        ambiguity_base = 'K'

    elif ref == 'G' and alt == 'T':
        ambiguity_base = 'K'

    elif ref == 'C' and alt == 'A':
        ambiguity_base = 'M'

    elif ref == 'A' and alt == 'C':
        ambiguity_base = 'M'

    elif ref == alt:
        ambiguity_base = ref
    
    return(ambiguity_base)

# Main function to loop over VCFs in input directory and create the respective consensus genomes

def main(args):
    
    REF = fasta2dict(args)
    EXCLUDED_POS = get_pos_to_exclude(args.bed_files)
    print(EXCLUDED_POS)

    # Define the vcf_dir from the argument -v/vcf, sort the vcf files alphanumerically, extract all keys (columns of depth file) alphanumerically sorted
    
    vcf_files = args.vcf_files
    

    # Loop over each vcf file with the matching column from the depth file
    for sample_idx, vcf_file in enumerate(vcf_files):

        VCF = vcf2dict([vcf_file])
        VCF_DEPTH = get_sample_depth(args.depth, sample_idx)


        fasta_sequence = []

        i = 1

        # Go over each position in the reference sequence and modify the base if it is in the VCF
        while(i <= len(REF)):

            if i in EXCLUDED_POS:   # Check if the position falls in the BED file

                fasta_sequence += ['N']
                i += 1
                

            else:

                if i not in VCF:    # Position not in VCF -> must be ancestral or deletion

                    if i in VCF_DEPTH:

                        if 1 <= VCF_DEPTH[i] <=5:       # Covered by less than 5 reads -> N
                            fasta_sequence += ['N']
                            i += 1
                        
                        elif VCF_DEPTH[i] == 0:         # Covered by 0 reads -> deletion
                            fasta_sequence += ['-']
                            i += 1

                        elif VCF_DEPTH[i] > 5:          # Not in VCF but covered >5 reads -> ancestral base
                            fasta_sequence += [REF[i]]
                            i += 1
                        
                elif i in VCF:  # Position is in VCF -> Variant

                    # for i define:
                    reference_base = VCF[i][0]
                    alternative_base = VCF[i][1]
                    quality = VCF[i][2]
                    allele_freq = VCF[i][3]

                    if quality >= 20:   # Variants with a 99% confidence (phred-score >20)

                        if len(allele_freq) == 1:  # Only one ALT allele

                            allele_freq = allele_freq[0]

                            if len(alternative_base) == 1 and len(reference_base) == 1: # REF and ALT = 1 -> SNP

                                if allele_freq >= 0.90: # Take alt base from VCF
                                    fasta_sequence +=[alternative_base]
                                    i += 1
                                    
                                elif 0.10 <= allele_freq < 0.90:    # Take ambiguity base
                                    fasta_sequence += [ambiguity_code(reference_base,alternative_base)]
                                    i += 1

                                elif allele_freq < 0.10:    # Take the ancestral base
                                    fasta_sequence += [REF[i]]
                                    i += 1
                            

                            elif len(alternative_base) < len(reference_base):   # Small deletions

                                if len(alternative_base) == 1 and len(reference_base) > 1:
                                    small_deletion_length = len(reference_base)    # Get length of REF to skip correct number of bases

                                    if allele_freq >=0.90:  # Encode with "-"
                                        fasta_sequence += [alternative_base[0]]
                                        fasta_sequence += split('-'*(small_deletion_length-1))
                                        i += small_deletion_length

                                    elif 0.10 <= allele_freq < 0.90:    # Take ambiguity base
                                        fasta_sequence += split('N'*small_deletion_length)
                                        i += small_deletion_length

                                    elif allele_freq < 0.10:    # take reference base
                                        fasta_sequence += [REF[i]]
                                        i += 1
                                
                                if len(alternative_base) > 1 and len(reference_base) > 1:   # If both REF & ALT are more than one base but still REF>ALT
                                    fasta_sequence += [REF[i]]
                                    i += 1

                            elif len(alternative_base) > len(reference_base):   # Small insertion
                                fasta_sequence += [REF[i]]
                                i += 1
                            
                            elif (len(alternative_base) > 1 ) and (len(reference_base) > 1) and (len(alternative_base) == len(reference_base)): # MNP

                                length_mnp = len(alternative_base)

                                if allele_freq >= 0.90: # Take alt bases from VCF
                                    fasta_sequence += split(alternative_base)
                                    i += length_mnp
                                    
                                elif 0.10 <= allele_freq < 0.90:    # Take ambiguity bases
                                    
                                    for f, b in zip(reference_base, alternative_base):
                                        fasta_sequence += [ambiguity_code(f, b)]
                                    i += length_mnp

                                elif allele_freq < 0.10:    # Take the ancestral bases
                                    fasta_sequence += split(reference_base)
                                    i += len(reference_base)

                        elif len(allele_freq) > 1:  # More than one ALT allele, we consider the ALT allele with the highest allele frequency

                            highest_AF = max(allele_freq)   # Find highest AF

                            alternative_base = alternative_base.split(',')  # Get different alleles split by ','

                            index_of_nucleotide_with_highest_AF = allele_freq.index(highest_AF) # Index of the allele with highest AF

                            ALT_allele_with_highest_frequency = alternative_base[index_of_nucleotide_with_highest_AF]   # ALT allele for highest AF

                            print(i,VCF[i],highest_AF,split(ALT_allele_with_highest_frequency),len(ALT_allele_with_highest_frequency),reference_base,REF[i])

                            # Now again check for SNPs, deletions, insertions and MNPs
                            if len(ALT_allele_with_highest_frequency) == 1 and len(reference_base) == 1: # SNP

                                if highest_AF >= 0.90: # Take ALT base
                                    fasta_sequence += [ALT_allele_with_highest_frequency]
                                    i += 1

                                elif 0.10 <= highest_AF < 0.90: # take ambiguity_base
                                    fasta_sequence += [ambiguity_code(reference_base,ALT_allele_with_highest_frequency)]
                                    i += 1

                                elif highest_AF < 0.10: # Take reference base
                                    fasta_sequence += [REF[i]]
                                    i +=1

                            elif len(ALT_allele_with_highest_frequency) < len(reference_base):  # Small deletion 

                                if len(ALT_allele_with_highest_frequency) == 1 and len(reference_base) > 1: # Small deletion with ALT allele = 1bp
                                    small_deletion_length = len(reference_base) # To move correct amount of bases

                                    if highest_AF >= 0.90:  # Encode with a '-'
                                        fasta_sequence += [ALT_allele_with_highest_frequency[0]]
                                        fasta_sequence += split('-'*(small_deletion_length-1))
                                        i += small_deletion_length

                                    elif 0.10 <= highest_AF < 0.90: # Encode with N
                                        fasta_sequence += split('N'*small_deletion_length)
                                        i += small_deletion_length

                                    elif highest_AF < 0.10: # Take reference base
                                        fasta_sequence += [REF[i]]
                                        i += 1

                                if len(ALT_allele_with_highest_frequency) > 1 and len(reference_base) > 1:  # Although ALT < REF, ALT > 1
                                    fasta_sequence += [REF[i]]
                                    i += 1

                            elif len(ALT_allele_with_highest_frequency) > len(reference_base): # Small insertion, no matter AF -> take anc base
                                fasta_sequence += [REF[i]]
                                i += 1

                            elif (len(ALT_allele_with_highest_frequency) > 1) and (len(reference_base) > 1) and (len(ALT_allele_with_highest_frequency) == len(reference_base)): # MNP

                                length_mnp = len(ALT_allele_with_highest_frequency)

                                if highest_AF >= 0.90:  # Take alt base
                                    fasta_sequence += split(ALT_allele_with_highest_frequency)
                                    i += length_mnp

                                elif 0.10 <= highest_AF < 0.90: # Take ambiguity base
                                    for f, b in zip(reference_base, ALT_allele_with_highest_frequency):
                                        fasta_sequence += [ambiguity_code(f, b)]
                                    i += length_mnp

                                elif highest_AF < 0.10: # Take reference base
                                    fasta_sequence += split(reference_base)
                                    i += len(reference_base)
                                    

                    else:   # QUAL < 20
                        fasta_sequence += ['N']
                        i += 1
                    
                else:   # Pos i not in VCF
                    print("i not in VCF",i,VCF[i])
                    break


        print("Lenght of sequence:",len(flatten(fasta_sequence)))

        # Find and print empty strings in the list fasta_sequence
        for b in range(len(fasta_sequence)):
            if len(fasta_sequence[b]) == 0:
                print(b, fasta_sequence[b])

        # Check if the script is running in test mode (shorter genomes)
        if args.test_mode:

            print("Running in test mode - smaller genomes accepted")
            reference_length = get_reference_length(args)

            # Created consensus sequence too short
            if len(fasta_sequence) < reference_length:
                print("Fasta file not created")
                print("Fasta has less than {reference_length} bp:", len(fasta_sequence))

            # Consensus sequence has correct length
            elif len(fasta_sequence) == reference_length:
                
                sample_name = get_basename(vcf_file)
                output_filename = f"{sample_name}.consensus.fasta"
                output_path = os.path.join(args.output_dir, output_filename)
                with open(output_path, 'w') as output_file:
                    output_file.write(f">{sample_name}\n")
                    output_file.write("".join(fasta_sequence))
        
            # Consensus sequence too long
            else:
                print("Fasta file not created")
                print("Fasta more than {reference_length} bp:", len(fasta_sequence))

        else: 
            # The genomes have to have the length of the MTB genome = 4411532 bp
            # Created consensus sequence too short
            if len(fasta_sequence) < 4411532:
                print("Fasta file not created")
                print("Fasta has less than 4411532 bp:", len(fasta_sequence))

            # Consensus sequence has correct length
            elif len(fasta_sequence) == 4411532:

                sample_name = get_basename(vcf_file)
                output_filename = f"{sample_name}.consensus.fasta"
                output_path = os.path.join(args.output_dir, output_filename)
                with open(output_path, 'w') as output_file:
                    output_file.write(f">{sample_name}\n")
                    output_file.write("".join(fasta_sequence))

        
            # Consensus sequence too long
            else:
                print("Fasta file not created")
                print("Fasta more than 4411532 bp:", len(fasta_sequence))
        

        


        
        
                    

                            
















