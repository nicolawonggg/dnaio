import dnaio


#Read chunks of complete FASTA or FASTQ records from a file

def chunk_single(fastq, buffer_size):
    '''
    PARAMETERS:
    - fastq (RawIOBase) – File with FASTA or FASTQ reads; must have been opened in binary mode
    - buffer_size (int) – Largest allowed chunk size
    '''
    return dnaio.read_chunks(fastq,buffer_size)

