Chaging directory to Read Counts directory
```bash
cd AmpSeq/4_Read_Counts/
```
Executing the perl script that will produce a counting of the read of a given amplicon given their primer sequences
```bash
singularity exec ../../AmpSeq.sig perl tag_presence.pl Files/ Primers.txt
```
