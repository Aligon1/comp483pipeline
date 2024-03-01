import os
import logging
import subprocess
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML

directory = 'PipelineProject_Alexa_Ligon'
if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)


def download_convert_fastq(SRR):
    '''Retrieves donor transcriptomes based on the SRR numbers and converts
    them to fastq files'''
    # Download the SRA file using wget
    download_command = ['wget', f'https://sra-download.st-va.ncbi.nlm.nih.gov/sos1/sra-pub-run-12/{SRR}.sra']
    subprocess.run(download_command)

    # Convert the SRA file to fastq using fastq-dump
    convert_command = ['fastq-dump', '-I', '--split-files', SRR]
    subprocess.run(convert_command)


def build_transcriptome_index_kallisto(accession):
    '''Builds the transcriptome index for Kallisto'''
    # Set the email address 
    Entrez.email = 'your email address'

    # Retrieve GenBank format file
    gb_handle = Entrez.efetch(db='nucleotide', id=accession, rettype='gb', retmode='text')
    cds_count = 0

    # Open log file
    with open('PipelineProject.log', 'a') as log:
        # Iterate over records and write CDS sequences to the output file
        with open(f'{accession}_CDS.fasta', 'w') as out_file:
            for record in SeqIO.parse(gb_handle, 'genbank'):
                for feature in record.features:
                    if feature.type == 'CDS':
                        cds_count += 1
                        protein_id = feature.qualifiers['protein_id'][0]
                        sequence = feature.location.extract(record).seq
                        out_file.write(f'>{protein_id}\n{sequence}\n')

        # Write to log file
        log.write(f'The HCMV genome {accession} has {cds_count} CDS.\n')

    # Build kallisto index
    subprocess.run(['kallisto', 'index', '-i', f'{accession}_index.idx', f'{accession}_CDS.fasta'])


def kallisto_quantify(SRR_list, input_path, accession):
    '''Quantify the TPM of each CDS in each transcriptome using Kallisto'''
    # write in header for results
    logging.info("sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm")

    for SRR in SRR_list:
        # run kallisto quantification
        subprocess.run(['kallisto', 'quant', '-i', f'{accession}_index.idx', '-o', f'./{SRR}', '-b', '30', f'{input_path}{SRR}_1.fastq', f'{input_path}{SRR}_2.fastq'])

        # read abundance.tsv file
        try:
            df = pd.read_csv(f'./{SRR}/abundance.tsv', sep='\t')

            # calculate min, median, mean, and max TPM
            min_tpm = df['tpm'].min()
            med_tpm = df['tpm'].median()
            mean_tpm = df['tpm'].mean()
            max_tpm = df['tpm'].max()

            condition1 = '2dpi'
            condition2 = '6dpi'

            condition = condition1 if int(SRR[3:]) % 2 == 0 else condition2
            # write to log file
            logging.info(f'{SRR}\t{condition}\t{min_tpm}\t{med_tpm}\t{mean_tpm}\t{max_tpm}')

        except FileNotFoundError:
            logging.error(f"Abundance file not found for {SRR}")


def sleuth_input(SRR_list):
    '''Generating input for Sleuth'''
    with open('covariate.txt', 'w') as covariate_file:
        condition1 = '2dpi'
        condition2 = '6dpi'
        covariate_file.write('sample\tcondition\tpath\n')  # initial line in file
        for SRR in SRR_list:
            path = f'./{SRR}'
            condition = condition1 if int(SRR[3:]) % 2 == 0 else condition2
            covariate_file.write(f'{SRR}\t{condition}\t{path}\n')


def run_sleuth(input_path):
    '''Runs sleuth in R and reads the sleuth output to add to the log file''' 
    subprocess.run(['Rscript', f'{input_path}sleuthscript.R'])

    # Read Sleuth output file
    try:
        with open('fdr05_results.txt', 'r') as output:
            for line in output:
                logging.info(line.strip())
    except FileNotFoundError:
        logging.error('Sleuth output file not found')


def extract_most_diff_exp(filename):
    '''From sleuth output, extract the target id of the protein with lowest qval'''
    # Read the file into a DataFrame
    df = pd.read_csv(filename, sep='\t')

    # Find the row with the lowest qval
    min_qval_row = df.loc[df['qval'].idxmin()]

    # Extract the target_id from the row
    target_id = min_qval_row['target_id']

    return target_id


def retrieve_protein_fasta(target_protein_id, output_file):
    ''' Function to retrieve protein fasta file from target protein id'''
    Entrez.email = 'enter your email'
    handle = Entrez.efetch(db='protein', id=target_protein_id, rettype='fasta', retmode='text')
    with open(output_file, 'w') as out_handle:
        out_handle.write(handle.read())
    handle.close()


def create_betaherpesvirinae_database(output_directory):
    '''make local database of Betaherpesvirinae nucleotide sequences'''
    Entrez.email = 'enter your email' 
    # Search for nucleotide sequences in the Betaherpesvirinae subfamily
    handle = Entrez.esearch(db='nucleotide', term='Betaherpesvirinae', retmax=100000, entrez_query='nr')  # Increase retmax if needed
    record = Entrez.read(handle)
    ids = record['IdList']

    try:
        # Create the output directory if it doesn't exist
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        # Fetch the sequences and write them to FASTA files in the output directory
        for id in ids:
            handle = Entrez.efetch(db='nucleotide', id=id, rettype='fasta', retmode='text')
            # Read the fetched sequence
            sequence_record = SeqIO.read(handle, 'fasta')
            if sequence_record:
                # Write sequence in FASTA format
                filename = os.path.join(output_directory, f'{sequence_record.id}.fasta')
                with open(filename, 'w') as outfile:
                    outfile.write(f'>{sequence_record.id}\n{sequence_record.seq}\n')
    except Exception as e:
        logging.error(f'Error fetching sequences: {e}')
    finally:
        handle.close()  # Close the handle after all sequences have been fetched and written


###

def run_blast(input_fasta, database_dir, output_xml):
    ''' Function to run Blast+'''
    # Concatenate all FASTA files in the database directory into a single file
    database_file = os.path.join(database_dir, 'betaherpesvirinae.fasta')
    with open(database_file, 'w') as db_handle:
        for filename in os.listdir(database_dir):
            if filename.endswith('.fasta'):
                record = SeqIO.parse(os.path.join(database_dir, filename), 'fasta')
                SeqIO.write(record, db_handle, 'fasta')

    # Run Blast+
    result_handle = NCBIWWW.qblast('blastp', db=database_file, sequence=input_fasta.read(), format_type='XML', hitlist_size=10)
    with open(output_xml, 'w') as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()

    # Remove the temporary database file
    #os.remove(database_file)


def parse_blast_xml(output_xml):
    '''Function to parse Blast+ XML output and write to log file'''
    result_handle = open(output_xml, 'r')
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            hsp = alignment.hsps[0]  # Keep only the best alignment (HSP)
            logging.info(f'{alignment.accession}\t{hsp.identities / hsp.align_length * 100:.2f}\t{hsp.align_length}\t{hsp.query_start}\t{hsp.query_end}\t{hsp.sbjct_start}\t{hsp.sbjct_end}\t{hsp.bits}\t{hsp.expect}\t{alignment.title}')
    result_handle.close()


def wrapper():
    '''Using SRR id number(s) as arguments for command line usage'''

    #SRR_list = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

    # SRR data files download and split (commented out after first use to save time)
    # if repo is cloned, test files can be found there    
    #download_convert_fastq(SRR_list)

    accession_num = 'NC_006273.2'  # accession number for index
    build_transcriptome_index_kallisto(accession_num)

    input_path = '/change/to/your/path/'  # change this to your path
    kallisto_quantify(SRR_list, input_path, accession_num)
    sleuth_input(SRR_list)
    run_sleuth(input_path)    

    # extracting the target id from sleuth table
    filename = 'fdr05_results.txt'
    target_id = extract_most_diff_exp(filename)

    # creating local non repeating nucleotide database for Betaherpesvirinae sub family
    output_dir = 'betaherpesvirinae'
    #create_betaherpesvirinae_database(output_dir)  # use file if possible

    # Define paths and filenames
    protein_fasta_file = f'{target_id}.fasta'
    blast_output_xml = 'blast_output.xml'
    retrieve_protein_fasta(target_id, protein_fasta_file)

    # Run Blast+ using the local database
    run_blast(open(protein_fasta_file), 'betaherpesvirinae', blast_output_xml)

    # Parse Blast+ XML output and write to log file
    parse_blast_xml(blast_output_xml)

    # Clean up
    os.remove(protein_fasta_file)
    os.remove(blast_output_xml)

if __name__ == '__main__':
    logFormatter = '%(message)s'
    logging.basicConfig(filename='PipelineProject.log', format=logFormatter, level=logging.INFO)
    wrapper()
