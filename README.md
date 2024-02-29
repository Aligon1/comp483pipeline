# Pipeline Project: Differential Expression

This pipeline executes the comparison of human herpesvirus 5 (HCMV) transcriptomes between two patient donors (Donor 1 and Donor 3) at two time points (2- and 6- days post infection (dpi))

### Dependencies

- Python
  - os
  - subprocess
  - Biopython: Entrez, SeqIO, Blast - NCBIWWW, NCBIXML
  - Pandas
- R
  - Sleuth 

### Executing program

**Before running:** <br>
- Update config.py for your specific configuration variables
- In the R script, make sure to set the correct path to the kallisto directory (line 94) 

**To run:** <br>
The script can be run by calling run_pipeline() at the end of the script or to run from commandline/terminal:

```
python pipeline_script.py
```
or <br>
```
python3 pipeline_script.py
```
