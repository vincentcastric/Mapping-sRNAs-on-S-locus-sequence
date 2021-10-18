#!usr/bin/env python

##############################################
#           WHAT THIS SCRIPT DOES        # 
##############################################

# This script takes as input a reference genome, a bed to maks the S-locus of the ref genome, a fasta file containing the focal S-allele and a .bed to mask s$
# It outputs an alignement file 


##############################################
#           MODULES                          # 
##############################################

import os
import os.path
import subprocess
import eep_modules 
import gzip
import sys
import re
from Bio import SeqIO
eep_modules.load("bedtools")
eep_modules.load("samtools")
eep_modules.load("bowtie")
eep_modules.load("trimmomatic")

##############################################
#           ARCHITECTURE ET LIENS            # 
##############################################

# A.lyrata genome
RefGenome = "./data/RefGenome/Alyrata_107_v9/Alyrata_107.fa"
# .bed file containing the coordinates of the S-locus in the reference genome
MaskingBed = "./data/RefGenome/slocuslyrata.bed"

#Slocus sequence 
SLocus_BAC = "./data/Slocus/Aha_21P24_Ah29.fas"
#SLocus_BAC = "./data/Slocus/BA_polished_contigs.fasta"
#SLocus_BAC = "./data/Slocus/Ah02_v2_Slocus.fa"
#SLocus_BAC = "./data/Slocus/Ah10_Slocus.fa"
#SLocus_BAC = "./data/Slocus/Ah20_Slocus.fa"

# .bed file containing the coordinates of the two S-locus flanking regions  
MaskingBedSlocus = "./data/Slocus/Aha_21P24_Ah29.bed"
#MaskingBedSlocus = "./data/Slocus/BA_polished_contigs.bed"
#MaskingBedSlocus = "./data/Slocus/Ah02_v2_Slocus.bed"
#MaskingBedSlocus  = "./data/Slocus/Ah10_Slocus.bed"
#MaskingBedSlocus = "./data/Slocus/Ah20_Slocus.bed"
SAllele = "Ah29"

# .fastqz file containing the raw sRNA-seq data
RawsRNA_reads = "./data/sRNAreads/data/210817_NB501473_A_L1-4_ALFV-4_R1.fastq.gz"
SLocusGenotype = "S01S29"

##############################################
#           COMMANDES                        # 
##############################################


#Mask the S-locus in the reference genome and index the masked genome
MaskedGenome = RefGenome.split(".fa")[0] + "_SLOCUSMASKED.fa"
if os.path.exists(MaskedGenome):
	print ("The masked file {0} already exists. Skipping this step.".format(MaskedGenome))
else:
	Command1 = "bedtools maskfasta -fi {0} -bed {1} -fo {2}".format(RefGenome, MaskingBed, MaskedGenome)
	os.system(Command1)
	Command2 = ("samtools faidx {0}".format(MaskedGenome))
	os.system(Command2)
	print ("Masked file {0} created and indexed.".format(MaskedGenome))

#Mask the  flanking regions in the sequence of the S-allele and index the masked sequence
Masked_SLocus = (SLocus_BAC.split(".fas")[0] + "_MASKED.fas")
if os.path.exists(Masked_SLocus):
	print ("The masked file {0} already exists. Skipping this step.".format(Masked_SLocus))
else: 
	Command3 = "bedtools maskfasta -fi {0} -bed {1} -fo {2}".format(SLocus_BAC, MaskingBedSlocus, Masked_SLocus)
	os.system(Command3)
	Command4 = ("samtools faidx {0}".format(Masked_SLocus))
	os.system(Command4)
	print ("Masked file {0} created and indexed.".format(Masked_SLocus))

#Trim the raw reads
TrimmedsRNA_reads = "./res/" + (RawsRNA_reads.split('.fastq.gz')[0]+"_TRIMMED.fastq.gz").split('data/sRNAreads/data/')[1]
if os.path.exists(TrimmedsRNA_reads):
	print ("The trimmed sRNA file {0} already exists. Skipping this step.".format(TrimmedsRNA_reads))
else: 
	print ("Trimming file {0} and creating {1}. Done.".format(RawsRNA_reads, TrimmedsRNA_reads))
	Command5 = ("trimmomatic SE -threads 24 -phred33 {0} {1} ILLUMINACLIP:./data/sRNAreads/adapters.ultimate.fa:2:30:10 AVGQUAL:25" .format(RawsRNA_reads, TrimmedsRNA_reads))
	os.system(Command5)

#Make a fasta file containing only non-redundant reads and a table containing their abundance
Trimmed_NR_sRNA_reads_txt = "./res/" + (RawsRNA_reads.split('.fastq.gz')[0]+"_TRIMMED_NR.txt").split('data/sRNAreads/data/')[1]
Trimmed_NR_sRNA_reads_fasta = "./res/" + (RawsRNA_reads.split('.fastq.gz')[0]+"_TRIMMED_NR.fasta").split('data/sRNAreads/data/')[1]
if os.path.exists(Trimmed_NR_sRNA_reads_txt):
        print ("The .txt file containing the set of non-redundant reads {0} already exists. Skipping this step.".format(Trimmed_NR_sRNA_reads_txt))
else:
	print("Busy creating the set of non-redundant reads...")
	f =  open('sequences.txt', 'w')
	fastq = SeqIO.parse(gzip.open(TrimmedsRNA_reads,'rt',encoding='utf-8'),"fastq")
	for rec in fastq:
		sequence = rec.seq
		print(sequence, file=f)
	f.close()
	print("Busy sorting the reads")
	Command6 = ("sort sequences.txt | uniq -c > {0}").format(Trimmed_NR_sRNA_reads_txt)
	os.system(Command6)
	Command7 = ("rm sequences.txt")
	os.system(Command7)

	print("Busy creating the FASTA file of non-redundant reads...")
	h = open(Trimmed_NR_sRNA_reads_fasta, 'w')
	lines = []
	with open (Trimmed_NR_sRNA_reads_txt, 'r') as g:
		lines = g.readlines()
	count = 0
	for line in lines:
		seq = line.split(" ")[-1].rstrip("\n")
		count = count+1
		occ = re.findall('\d+', line)[0]
		print(">sRNA{0}_occ{1}".format(count,occ), file = h)
		print(seq, file = h)
	h.close()
	print("... done !")


#Map the non-redundant reads on the reference genome, and map the reads on the S-locus
SlocusBowtieIndex = SLocus_BAC.split(".fa")[0]
AlyrataBowtieIndex = MaskedGenome.split(".fa")[0]
AlignmentOnSlocus =  "./res/" + (RawsRNA_reads.split("L1-4_")[1]).split("_R1")[0] + "_" + SAllele + ".sam"
AlignmentOnLyratagenome = "./res/" + (RawsRNA_reads.split("L1-4_")[1]).split("_R1")[0] + "_AlyGenome.sam".format()

if os.path.exists("{0}.1.ebwt".format(SlocusBowtieIndex)):
	print("An index for {0} already exists. Skipping this task...".format(SAllele))
else:
	print("Indexing {0} for Bowtie.".format(SAllele))
	command8 = ("bowtie-build {0} {1}".format(SLocus_BAC,SlocusBowtieIndex))
	os.system(command8)
if os.path.exists("{0}.1.ebwt".format(AlyrataBowtieIndex)):
	print("An index for the A. lyrata genome  already exists. Skipping this task...")
else:
	print("Indexing the A. lyrata genome for Bowtie.")
	command9 = ("bowtie-build {0} {1}".format(MaskedGenome, AlyrataBowtieIndex))
	os.system(command9)

if os.path.exists(AlignmentOnSlocus):
        print ("Mapping of sRNA reads on {0} has been done already. Skipping this task...".format(SAllele))
else:
	print("Launching Bowtie on {0} ...".format(SAllele))
	command10 = ("bowtie {0} -v 0 -a -M 10 --best --strata -y -S -p 4 -f {1} | samtools view -h -F 4 -S -  1>{2}").format(SlocusBowtieIndex, Trimmed_NR_sRNA_reads_fasta, AlignmentOnSlocus)
	os.system(command10)

if os.path.exists(AlignmentOnLyratagenome):
        print ("Mapping of sRNA reads on the A.lyrata genome has been done already. Skipping this task...") 
else:
	print("Launching Bowtie on the A. lyrata genome...")
	command11 = ("bowtie {0} -v 0 -a -M 10 --best --strata -y -S -p 4 -f {1} | samtools view -h -F 4 -S -  1>{2}").format(AlyrataBowtieIndex, Trimmed_NR_sRNA_reads_fasta, AlignmentOnLyratagenome)
	os.system(command11)


#Exclude the reads that are shared between these two libraries and flag those that are uniquely mapping to the S-locus  
UniqueAlignmentOnSlocus = AlignmentOnSlocus.split(".sam")[0] + "_{0}_UNIQ.sam".format(SLocusGenotype)
UniqueAlignmentOnSlocusBAM =  AlignmentOnSlocus.split(".sam")[0] + "_{0}_UNIQ.bam".format(SLocusGenotype)
UniqueSortAlignmentOnSlocusBAM =  AlignmentOnSlocus.split(".sam")[0] + "_{0}_UNIQSORT.bam".format(SLocusGenotype)
NotUniqueAlignmentOnSlocus = AlignmentOnSlocus.split(".sam")[0] + "_{0}_NOTUNIQ.sam".format(SLocusGenotype)
NotUniqueAlignmentOnSlocusBAM = AlignmentOnSlocus.split(".sam")[0] + "_{0}_NOTUNIQ.bam".format(SLocusGenotype)
NotUniqueSortAlignmentOnSlocusBAM = AlignmentOnSlocus.split(".sam")[0] + "_{0}_NOTUNIQSORT.bam".format(SLocusGenotype)

dictlyr = {}
with open (AlignmentOnSlocus, 'r') as g, open(AlignmentOnLyratagenome, 'r')  as lyr, open(UniqueAlignmentOnSlocus, 'w') as uniq, open(NotUniqueAlignmentOnSlocus, 'w') as notuniq:
	lineslyr = lyr.readlines()
	for line in lineslyr:
		sRNA_id = line.split("\t")[0].rstrip("\n")
		dictlyr[sRNA_id] = 1

	lines = g.readlines()
	for line in lines:
		if line[0] == "@":
			print(line.rstrip("\n"), file = notuniq)
			print(line.rstrip("\n"), file = uniq)
		else:
			sRNA_id = line.split("\t")[0].rstrip("\n")
			if sRNA_id in dictlyr.keys():
				print(sRNA_id)
				print("not unique")
				print(line.rstrip("\n"), file = notuniq)
			else:
				print(sRNA_id)
				print("unique")
				print(line.rstrip("\n"), file = uniq)

print("converting sam to bam")
command12 = ("samtools view -b {0} -o {1}").format(UniqueAlignmentOnSlocus, UniqueAlignmentOnSlocusBAM)
os.system(command12)
command13 = ("samtools view -b {0} -o {1}").format(NotUniqueAlignmentOnSlocus, NotUniqueAlignmentOnSlocusBAM)
os.system(command13)

print("sorting bam files")
command14 = ("samtools sort {0} -o {1}").format(UniqueAlignmentOnSlocusBAM, UniqueSortAlignmentOnSlocusBAM)
os.system(command14)
command15 = ("samtools sort {0} -o {1}").format(NotUniqueAlignmentOnSlocusBAM, NotUniqueSortAlignmentOnSlocusBAM)
os.system(command15)

print("indexing bam files")
command16 = ("samtools index {0}").format(UniqueSortAlignmentOnSlocusBAM)
os.system(command16)
command17 = ("samtools index {0}").format(NotUniqueSortAlignmentOnSlocusBAM)
os.system(command17)

print("...done !")

