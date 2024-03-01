# Pipeline Project: Differential Expression

**PipelineProject.log is not complete, had computer problems - issue with how I was making the database for the betaherpesvirinae. All scripts are there.**

This pipeline executes the comparison of human herpesvirus 5 (HCMV) transcriptomes between two patient donors (Donor 1 and Donor 3) at two time points (2- and 6- days post infection (dpi))

### Dependencies

- Python3
  - SRAtoolkit
  - logging
  - os
  - subprocess
  - Biopython: Entrez, SeqIO, Blast - NCBIWWW, NCBIXML
  - Pandas
- Blast+ - See NCBI for download
- Kallisto - Documentation: https://pachterlab.github.io/kallisto/about
- R
  - Sleuth - Documentation: https://pachterlab.github.io/sleuth/about 


### Executing program

**To run this script, please clone this github repository with:**
```
git clone https://github.com/Aligon1/comp483pipeline.git
```

**Before running:** <br>
- Enter your email for Entrez operations
- Enter your 'input_path' in wrapper
- Enter the NCBI accession number in wrapper for transcriptome indexing
 
**To run:**
```
python3 wrapper.py 
```

### Files in repo

- wrapper.py - the main pipeline file
- sample subset files:
  - SRR5660030_1_subset.fastq
  - SRR5660030_2_subset.fastq
  - SRR5660033_1_subset.fastq
  - SRR5660033_2_subset.fastq
  - SRR5660044_1_subset.fastq
  - SRR5660044_2_subset.fastq
  - SRR5660045_1_subset.fastq
  - SRR5660045_2_subset.fastq
- PipelineProject.log - main pipeline results
- local betaherpesvirinae 
