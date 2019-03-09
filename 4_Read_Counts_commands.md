Chaging directory to Read Counts directory
```bash
cd AmpSeq-Pre/4_Read_Counts/
```
Copying fastq.gz files files to be processed
```bash
cp ../Original_Files/* Files/
```
Changing directory to the Files directory
```bash
cd Files/
```
Executing the perl script that will produce a counting of the read of a given amplicon given their primer sequences
```bash
singularity exec ../../../AmpSeq.sig perl ../tag_presence.pl ../Primers.txt
```
