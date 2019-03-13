Log in your HPC service and change directory (cd) to the working directory

Get the information from the AmpSeq GitHub repository by executing:
```bash
git clone https://github.com/JFresnedo/AmpSeq.git
```
Get the Singularity container with all the software needed
```bash
singularity pull --name AmpSeq.sig shub://JFresnedo/AmpSeq:ampseq
# You will see a warning which you can disregard, actually, a progress bar will show you the download progress of the Singularity container
```
Load Singularity in your HPC service

Test the container works by executing a simple command:
```bash
singularity exec vcftools
```
