# Pipeline Project: Differential Expression

This pipeline executes the comparison of human herpesvirus 5 (HCMV) transcriptomes between two patient donors (Donor 1 and Donor 3) at two time points (2- and 6- days post infection (dpi))

### Dependencies

- Python3
  - SRAtoolkit
  - logging
  - os
  - subprocess
  - Biopython: Entrez, SeqIO, Blast - NCBIWWW, NCBIXML
  - Pandas
- Kallisto - Documentation: https://pachterlab.github.io/kallisto/about
- R
  - Sleuth - Documentation: https://pachterlab.github.io/sleuth/about 


### Executing program

**To run this script, please clone this github repository with:**
```
git clone https://github.com/Aligon1/comp483pipeline.git
```

**Before running:** <br>
- Make sure to enter your email for Entrez 
 
**To run:**
```
python3 wrapper.py 
```
