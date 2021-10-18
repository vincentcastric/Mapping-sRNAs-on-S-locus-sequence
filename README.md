# Mapping-sRNAs-on-S-locus-sequence
This script maps sRNAs on a S-locus sequence and a reference genome, and returns separate alignment files for sRNAs that are only mapping to the S-locus ("unique") and sRNAs that also map elsewhere in the genome ("non-unique")

The following input files are required : 
./data/Slocus/S_locus_sequence.fasta (the raw sequence, eg. of a BAC clone or a genome assembly scaffold)
./data/Slocus/S_locus_sequence.bed (the masking file to mask everything in the fasta file that is outside the S-locus)
./data/RefGenome/AlyrataGenome.fasta (the raw sequence of the reference genome, e.g. of Arabidopsis lyrata)
./data/RefGenome/AlyrataGenome.bed (the masking file to mask the S-locus in the assembly)
./data/sRNAreads/data.fastq.gz (the compressed .gz archive containing the raw sRNA reads)
./sRNAreads/adapters.ultimate.fa (the Illumina adapter sequence)
