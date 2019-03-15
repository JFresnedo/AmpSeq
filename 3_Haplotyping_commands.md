Changing directory to Haplotyping directory
```bash
cd AmpSeq/3_Haplotyping/
```
Executing the routine for haplotyping based on Picard and ClustalOmega (step 8 is skipped, see the wiki for details)
```bash
singularity exec ../../AmpSeq.sig perl analyze_amplicon.pl -s sample_file.txt -k key_file.txt -o Output -m clustalo:clustal -i 8
```
