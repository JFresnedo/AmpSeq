Changing directory to SNP calling
```bash
cd AmpSeq/2_SNP_calling/
```
Copying fastq.gz files to be processed through the pipeline
```bash
cp ../Original_Files/Pa* Files/ & cp ../Original_Files/Pr* Files/
```
Listing the just-copied files in the directory Files
```bash
ls Files/
```
Chaning directory to Indexed_Genome in order to pre-process reference genome infromation (generation of indices)
```bash
cd Indexed_Genome
```
Concatenating fasta files with sequences of chromosmes 9, 12, 14, 14 and 19 of grapevine
```bash
singularity exec ../../../AmpSeq.sig cat *.gz | gzip -d > grapevine_12Xv2_Chr9_12_13_14_19.fa
```
Indexing the fasta file with chromosme sequences using BWA
```bash
singularity exec ../../../AmpSeq.sig bwa index grapevine_12Xv2_Chr9_12_13_14_19.fa
```
Indexing the fasta file with chromosme sequences using samtools
```bash
singularity exec ../../../AmpSeq.sig samtools faidx grapevine_12Xv2_Chr9_12_13_14_19.fa
```
Indexing the fasta file with chromosme sequences using Picard
```bash
singularity exec ../../../AmpSeq.sig java -jar /usr/local/bin/picard/build/libs/picard.jar CreateSequenceDictionary R=grapevine_12Xv2_Chr9_12_13_14_19.fa O=grapevine_12Xv2_Chr9_12_13_14_19.dict
```
Returning to 2_SNP_calling directory
```bash
cd ..
```
Executing the SNp calling routine based on GATK
```bash
singularity exec ../../AmpSeq.sig perl run_gatk2.pl Files/ Output Indexed_Genome/grapevine_12Xv2_Chr9_12_13_14_19.fa
```
