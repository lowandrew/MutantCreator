# Mutant Creator

The `MutantCreator.py` script provides functionality for creating an fasta that differs by
a set number of SNVs compared to an input fasta file. It can also create paired-end FASTQ files
that correspond to the mutant genome created if you want to do any read analysis on them.

## Installation

You'll need to have the biopython module installed, and have the `art_illumina` executable from ART
available on your $PATH if you want to simulate reads based on your mutant genome.

Clone this repository: `git clone https://github.com/lowandrew/MutantCreator.git`

You can then access the script from the cloned repository.

## Usage

Generate a fasta file called mutant.fasta with 1000 SNVs when compared to input.fasta: `python MutantCreator.py -i input.fasta -o mutated.fasta -n 1000`

Generate FASTQ files to a depth of 30X of a genome 200 SNVs away from input.fasta: `python MutantCreator.py -i input.fasta -o mutant.fasta -fq mutant -n 200 -d 30`

