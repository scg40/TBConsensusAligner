import os
import glob
import argparse
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
import gzip
import consensus_galaxy
import snp_aligner_galaxy

'''

This script is the main script for the "insert name here". 
It is built on top of consensus.py and snp_aligner.py.

The user can choose the following with the -s option.
-s all:         Input = VCFs, output = SNP alignment -> runs consensus.py and snp_aligner.py subsequently.
-s consensus:   Input = VCFs, output = consensus fasta files -> runs just consensus.py
-s alignment:   Input = consensus fasta files, output = SNP alignment -> runs just snp_aligner.py.
When -s all is chosen, -m option gives choice whether to just output the SNP alignment
with -m alignment_only or both the SNP alignment and the consensus fasta files with  -m everything

'''

# Define arguments used in the script
def get_args():

    parser = argparse.ArgumentParser(description='Main script for consensus.py and snp_aligner.py')

    # Use dest to set the name of the argument for further handling and clarity 
    parser.add_argument('-s', choices=['consensus', 'all'],dest='step', help='Run either consensus, or both consensus and snp_aligner', required=True)
    parser.add_argument('-m', choices=['alignment_only', 'everything'],dest='mode', help='If -s all then output either just SNP alignment or both the alignment and the consensus fasta files', required=False)
    parser.add_argument('-v', action='append', dest='vcf_files',help='path to the input directory with all the vcf files', required=True)
    parser.add_argument('-r', dest='reference',help='path to reference genome file', required=True)
    parser.add_argument('-d',dest='depth',help='Depth file per position, output of samtools depth', required= True)
    parser.add_argument('-c', dest='outgroup_vcf', help='Outgroup VCF required for variable alignment', required=False)
    parser.add_argument('-n', dest='outgroup_name', type=str, help='Clean name for outgroup header in Galaxy')
    parser.add_argument('-b', dest='bed_files', action='append', default=[], help='Optional BED files to mask certain genomic regions', required=False)
    parser.add_argument('-o', dest='output_dir', help='Output directory for consensus files and the variable alignment' ,required=False)
    parser.add_argument('-g', dest='undefined_states', help='Percentage of undefined states allowed per polymorphic position in the alignment', type=float, default=0.9, required=False)
    parser.add_argument('-t', dest='test_mode',help='allows to use smaller files with less genomic positions', action='store_true')


    # Run the parser and place data in the parser object for later use
    args = parser.parse_args()

    return args


# Check if all arguments for the consensus script are there
def check_consensus_args(args):
    
    required_consensus_args = ['vcf_files', 'reference', 'depth']

    # List comprehension to iterate over field in required_consensus_args and if the field is missing, adds it to the list
    missing_consensus_args = [field for field in required_consensus_args if getattr(args, field, None) is None]

    # If Missing_consensus_args has elements, it raises and error and outputs what arguments are missing, separated by a ','
    if missing_consensus_args:
        raise ValueError(f"Missing required arguments: {', '.join(missing_consensus_args)}")

# Check if all arguments for the snp_aligner script are there
def check_snp_aligner_args(args):

    required_snp_aligner_args = ['outgroup_vcf', 'undefined_states']

    # List comprehension to iterate over field in required_snp_aligner_args and if the field is missing, adds it to the list
    missing_snp_aligner_args = [field for field in required_snp_aligner_args if getattr(args, field, None) is None]

    # If missing_snp_aligner_args has elements, it raises and error and outputs what arguments are missing, separated by a ','
    if missing_snp_aligner_args:
        raise ValueError(f"Missing required arguments: {', '.join(missing_snp_aligner_args)}")
    
# Main logic of the program
def main():

    args = get_args()

    if args.output_dir is None:
        args.output_dir = "output"

    os.makedirs(args.output_dir, exist_ok=True)

    # Find input files and attribute them to arguments
    #args.reference, args.depth, args.outgroup_vcf, args.vcf_files, args.bed_files = find_input_files(args.input_dir)
    #print("REFERENCE:", args.reference)
    #print("DEPTH:", args.depth)
    #print("OUTGROUP:", args.outgroup_vcf)
    #print("VCFS:", args.vcf_files)
    #print("BEDS:", args.bed_files)


    # -s consensus: Input = VCFs, output = consensus fasta files -> runs just consensus.py
    if args.step == 'consensus':

        # Adress output paths
        # args.consensus_dir = os.path.join(args.output_dir, "consensus")
        # os.makedirs(args.consensus_dir, exist_ok=True)

        check_consensus_args(args)

        consensus_galaxy.main(args)

    # -s all: Input = VCFs, output = SNP alignment -> runs consensus.py and snp_aligner.py subsequently
    elif args.step == 'all':

                # Adress output paths
        # args.consensus_dir = os.path.join(args.output_dir, "consensus")
        # os.makedirs(args.consensus_dir, exist_ok=True)
        
        check_consensus_args(args)
        check_snp_aligner_args(args)

        consensus_galaxy.main(args)
        snp_aligner_galaxy.main(args)

                # Adress output paths
        # args.alignment_dir = os.path.join(args.output_dir, "alignment")
        # os.makedirs(args.alignment_dir, exist_ok=True)

        

        

        # After consensus is run and -m alignment_only is chosen, delete the fasta files again
        if args.mode == 'alignment_only':

            for file in glob.glob(os.path.join(args.output_dir, "*.fasta")):
                os.remove(file)

if __name__ == '__main__':
    main()