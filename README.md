# TBConsensusAligner
TBConsensusAligner is a tool designed to create a multi-sequence alignment via consensus 
genomes created from VCF files. 
It can currently be used in different modes as follows
- VCF to SNP alignment `-s all`
    - input = VCFs, output = SNP alignment
    - runs `consensus_galaxy.py` and `snp_aligner_galaxy.py`

- VCF to consensus FASTA files `-s consensus`
    - input = VCFs, output = consensus FASTA files
    - runs only `consensus_galaxy.py`

If the script produces a SNP alignment from VCFs, the user can choose whether to ouptut just the SNP alignment
or the used consensus FASTA files as well with the `-m` option. 

When running the `consensus_galaxy.py` the user has to provide the VCF files, the reference genome used in their creation and the multisample
depth file from bamtools depth. Optionally, the user can provide one or several BED files to mask certain regions of the genome.

## Usage
The script is run via the command line.

***usage***:

```
Usage:   TBConsensusAligner [options] 

Options: -s STR         ['consensus' or 'all'] Create either consensus files or consensus files and a variable alignment

         -m STR         ['everything' or 'alignment_only'] If -s all is chosen, output either consensus files and the alignment or just the alignment

         -v list[STR]   Path to the input VCF files

         -r STR         Path to the reference genome

         -d STR         Path to the depth file created by samtools depth

         -c STR         Path to the outgroup VCF for the alignment

         -b list[STR]   Path to the BED files to mask certain genomic regions

         -o STR         Path to the output directory

         -g FLOAT       Percentage of undefined states allowed per polymorphic position in the alignment default at 90 percent value=0.9

         -t BOOLEAN     Test Mode to run the script with smaller genomes

```

## Example usages command line

### From 3 VCFs to alignment, output consensus FASTAs and SNP alignment, masked by two bedfiles 

```
python yourfolder/galaxy_main.py \
-s all \
-m everything \
-v path/to/vcf1/vcf1.vcf -v path/to/vcf2/vcf2.vcf -v path/to/vcf2/vcf2.vcf \
-r path/to/reference/reference_genome.reference.fasta \
-d path/to/depthfile/depthfile.tabular \
-c path/to/outgroupvcf/outgroup.vcf \
-b path/to/bedfile1/bed1.bed -b path/to/bedfile2/bed2.bed \
-o output
-g 0.9

```

## Logic consensus.py 

This script produces consensus FASTA files from VCFs. The user provides the VCFs, the reference genome,
the multisample depth file obtained from `bamtools depth` and optionally BED files to mask certain genomic regions.

***algorithm***

We loop through each position in the reference genome and build the consensus sequence base per base with the following rules:

```
- SNPs and MNPs with a frequency >90% are encoded with the alternative base from the VCF.
- SNPs and MNPs with a frequency between 10% and 90% are encoded with the ambiguity base.
- SNPs and MNPs with a frequency <10% are encoded with the ancestral base from the reference genome.

- Large deletions (coverage 0) are encoded by dashes (-).
- Small deletions with frequency >90% are encoded by dashes (-).
- Small deletions with frequency between 10% and 90% are encoded with a 'N'.
- Small deletions with frequency <10% are encoded with the ancestral base from the reference genome.

- Small insertions (alternative bases longer than reference base) are encoded with the ancestral state from the reference genome.

- Sites to exclude specified in the bed files are encoded with a 'N'.
- Variants that have a phred score < 20 are encoded with a 'N'.
- Sites that are not in the VCF and are covered by less than 5 reads (via depth file) are encoded with a 'N'.

- If more than one ALT allele is present, we consider the one with the highest allele frequency.
```
Since the VCF for which the script is tailored to does not have the allele frequency `AF` calculated, we calculate it using 
`RO = REF allele occurance` and `AO = ALT allele occurance` with 
`AF = AO / ( sum(all AO's) + RO )`.

## Logic snp_aligner.py

This script produces a SNP alignment i.e. multi-sequence alignment of polymorphic positions in FASTA format from consensus FASTA files.
The consensus FASTAs are generated in the previous step by `consensus_galaxy.py` wrapped in `TBConsensusAligner`.

***algorithm***

The script creates a multi-sequence alignment from the consensus FASTAs in form of an 2D array where each row represents a sequence
and each column a genomic position. Each column is checked for the abundance of a polymorphism. If such a polymorphism is detected,
the column is appended to the alignment of polymorphic positions. For each polymorphic position, the corresponding nucleotide from 
outgroup VCF is retrieved to get the outgroup's sequence.

The algorithm to populate the multi-sequence alignment of polymorphic positions is asd follows.
The user can control the proportion of gaps or undefined states (`-`or`N`) for a polymorphic position
to be kept in the alignment (argument `-g`). By default, this is set to `0.9`, meaning that if a polymorphic position has 
more than 90% gaps or undefined states, it will not be in the final alignment.

## Directory structure used for development

* TBConsensusAligner
    * test_data
        - snp_alignment.fasta
        - test_G77777.consensus.fasta
        - test_G77777.vcf.gz
        - test_G88888_k1.consensus.fasta
        - test_G88888_k1.vcf.gz
        - test_G99999.k2.consensus.fasta
        - test_G99999.k2.vcf.gz
        - test_Galaxy_multiple_depths_header.tabular
        - test_reference_200bp.reference.fasta
        - test_regions_blindspots_modlin_farhat_and_PE_PPE_PGRS.bed
        - test.outgrou.all.pos.vcf.gz
    - consensus_galaxy.py
    - main_galaxy.py
    - README.md
    - snp_aligner_galaxy.py
    - TBConsensusAligner.xml



