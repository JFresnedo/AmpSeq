Changing directory to Haplotyping directory
```bash
cd AmpSeq-Pre/3_Haplotyping/
```
Copying fast.gz files into the directory to proccess them
```bash
cp ../Original_Files/* .
```
Executing the routine for haplotyping based on Picard and ClustalOmega (step 8 is skipped)
```bash
singularity exec ../../AmpSeq.sig perl analyze_amplicon.pl -s sample_file.txt -k key_file.txt -o Output -m clustalo:clustal -i 8
```
