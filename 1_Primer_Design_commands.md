Changing directory to Primer Desing
```bash
cd AmpSeq-Pre/1_Primer_Design/
```
Decompressing SAM file
```bash
singularity exec ../../AmpSeq.sig bgzip -dc Horizon_x_Illinois_547_1.sam.gz > Horizon_x_Illinois_547_1.sam
```
Producing a BAM file
```bash
singularity exec ../../AmpSeq.sig samtools view -bS -o Horizon_x_Illinois_547_1.bam Horizon_x_Illinois_547_1.sam
```
Sorting BAM file
```bash
singularity exec ../../AmpSeq.sig samtools sort Horizon_x_Illinois_547_1.bam -o Horizon_x_Illinois_547_1.sorted.bam
```
Indexing sorter BAM file
```bash
singularity exec ../../AmpSeq.sig samtools index Horizon_x_Illinois_547_1.sorted.bam
```
Decompressing VCF file
```bash
singularity exec ../../AmpSeq.sig bgzip -dc Horizon_x_Illinois_547_1.vcf.gz > Horizon_x_Illinois_547_1.vcf
```
Decompressing fasta file with the reference genome (in this case, chromosme 14 of grape genome)
```bash
singularity exec ../../AmpSeq.sig bgzip -dc grapevine_12Xv2_Chr14.fa.gz > grapevine_12Xv2_Chr14.fa
```
Indexing fasta file
```bash
singularity exec ../../AmpSeq.sig samtools faidx grapevine_12Xv2_Chr14.fa
```
Calculation of Linkage Disequilibrium (correlation) between markers in the region of interest with respect to flanking markers and anchor marker (in test.marker.txt)
```bash
singularity exec ../../AmpSeq.sig perl calculated_LD_distribution.pl test_marker.txt Horizon_x_Illinois_547_1.vcf
```
Parsing sequencing data in the region of interest to identify additional SNPs
```bash
singularity exec ../../AmpSeq.sig perl parse_bam.pl Horizon_x_Illinois_547_1.sorted.bam grapevine_12Xv2_Chr14.fa test_marker.txt_statistic.out 1e-15
```
Design or primers for amplicons of 45bp considering the SNPs markers identified from the BAM file
```bash
singularity exec ../../AmpSeq.sig perl primer3.pl grapevine_12Xv2_Chr14.fa selected_sites 22
```
