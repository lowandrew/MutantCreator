import subprocess
import argparse
import textwrap
import random
import os
from Bio import SeqIO
from accessoryFunctions import accessoryFunctions


def mutate_site(sequence, site):
    if site > len(sequence) or site < 0:
        raise ValueError('Specified site is not within sequence length.')
    current_nucleotide = sequence[site]
    new_nucleotide = current_nucleotide
    while current_nucleotide == new_nucleotide:
        random_number = random.randint(0, 3)
        if random_number == 0:
            new_nucleotide = 'A'
        elif random_number == 1:
            new_nucleotide = 'C'
        elif random_number == 2:
            new_nucleotide = 'G'
        elif random_number == 3:
            new_nucleotide = 'T'
    new_sequence = sequence[:site] + new_nucleotide + sequence[site + 1:]
    return new_sequence


def get_sequence_from_fasta(fasta_file):
    total_sequence = ''
    contigs = SeqIO.parse(fasta_file, 'fasta')
    for contig in contigs:
        contig_sequence = str(contig.seq)
        contig_sequence = contig_sequence.upper()  # Non-uppercase characters could potentially cause problems later.
        total_sequence += contig_sequence
    return total_sequence


def mutate_genome(genome_sequence, num_sites_to_mutate):
    mutated_sites = list()
    while len(mutated_sites) < num_sites_to_mutate:
        site_to_mutate = random.randint(0, len(genome_sequence))
        if site_to_mutate not in mutated_sites:
            genome_sequence = mutate_site(genome_sequence, site_to_mutate)
            mutated_sites.append(site_to_mutate)
    return genome_sequence


def write_mutated_genome_to_fasta(genome_sequence, output_fasta):
    with open(output_fasta, 'w') as outfile:
        outfile.write('>sequence\n')
        seq = textwrap.fill(genome_sequence)
        outfile.write(seq + '\n')


def create_fastq_from_fasta(fasta_file, output_fastq, coverage_depth):
    cmd = 'art_illumina -ss MSv1 -i {mutated_fasta} -l 250 -na -p -f {depth} -m 350 -s 10' \
          ' -o {output_fastq}'.format(mutated_fasta=fasta_file, depth=coverage_depth, output_fastq=output_fastq)
    subprocess.call(cmd, shell=True)
    # Move the files created to a more normal naming scheme.
    cmd = 'mv {base_name}1.fq {base_name}_R1.fastq'.format(base_name=output_fastq)
    os.system(cmd)
    cmd = 'mv {base_name}2.fq {base_name}_R2.fastq'.format(base_name=output_fastq)
    os.system(cmd)
    # Gzip the files.
    cmd = 'gzip {base_name}_R1.fastq {base_name}_R2.fastq'.format(base_name=output_fastq)
    os.system(cmd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_fasta',
                        type=str,
                        required=True,
                        help='Base Fasta file you would like to have a mutant created for.')
    parser.add_argument('-o', '--output_fasta',
                        type=str,
                        required=True,
                        help='Output file for your mutated sequence.')
    parser.add_argument('-n', '--number_of_mutations',
                        type=int,
                        default=500,
                        help='Number of sites to mutate.')
    parser.add_argument('-fq', '--fastq_base',
                        default=None,
                        help='Base name for an output FASTQ file.')
    parser.add_argument('-d', '--coverage_depth',
                        default=60,
                        type=int,
                        help='Depth of coverage desired for FASTQ file.')
    args = parser.parse_args()
    accessoryFunctions.dependency_check('art_illumina')
    # Extract reference sequence from FASTA
    reference_sequence = get_sequence_from_fasta(args.input_fasta)
    # Mutate the input genome.
    new_sequence = mutate_genome(reference_sequence, args.number_of_mutations)
    # Write output genome to file.
    write_mutated_genome_to_fasta(new_sequence, args.output_fasta)
    # Simulate reads, if user wanted that.
    if args.fastq_base:
        create_fastq_from_fasta(args.output_fasta, args.fastq_base, args.coverage_depth)
