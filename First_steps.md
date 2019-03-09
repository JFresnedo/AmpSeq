Log in your HPC service and change directory (cd) to the working directory

Get the information from the AmpSeq GitHub repository by executing:
```bash
git clone https://github.com/JFresnedo/AmpSeq-Pre.git
```
Get the Singularity container with all the software needed
```bash
singularity pull --name AmpSeq.sig shub://JFresnedo/AmpSeq:ampseq
```
Load Singularity in your HPC service

Test the container works by executing a simple command:
```bash
singularity exec vcftools
```
