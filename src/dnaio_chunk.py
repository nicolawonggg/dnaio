# from dnaio import read_chunks, read_paired_chunks
import random

import dnaio
from xopen import xopen
from multiprocessing import Pool
from functools import partial


#Read chunks of complete FASTA or FASTQ records from a file
def chunk_single(fastq, buffer_size):
    '''
    RETURNS:
        generator
    PARAMETERS:
        fastq (RawIOBase) – File with FASTA or FASTQ reads; must have been opened in binary mode
        buffer_size (int) – Largest allowed chunk size
    '''
    with dnaio.open(fastq, mode="r") as f1:
        i = 0
        for chunk in dnaio.read_chunks(f=f1, buffer_size=buffer_size):
            i += 1
            list_fn = fastq.split(".")
            filename = list_fn[0].join(f"_{i}").join(list_fn[1:])
            print("check fn: {filename}")

            with open(filename, mode="w") as writer:
                writer.write(chunk)
                print(f"file saved: {filename}")


def parse_chunk(chunk_byte):
    """
    This function parses all reads (in bytes) per chunk into DNAIO SequenceRecord objects
    Arguments:  chunk_byte(RawIOBase) from memoryview objects returned by  dnaio.read_paired_chunks
    Yields:     SequenceRecord
    """
    def reset():
        name = '';
        sequence = '';
        smth = '';
        qualities = ''
        return name, sequence, smth, qualities
    chunk_seq_records = []
    line_count = 0
    name, sequence, smth, qualities = reset()
    for r in chunk_byte.split(b'\n'):
        line_count += 1
        if line_count % 4 == 1:
            name = r.decode("utf-8")
        elif line_count % 4 == 2:
            sequence = r.decode("utf-8")
        elif line_count % 4 == 3:
            smth = r.decode("utf-8")
        elif line_count % 4 == 0:
            qualities = r.decode("utf-8")
            yield dnaio.SequenceRecord(name, sequence, qualities)
            # print(f"name: {name}, seq: {sequence}, smth: {smth}, qualities: {qualities}")
            name, sequence, smth, qualities = reset()


def chunk_pair(fastq1, fastq2, buffer_size):
    """
    This function chunks a pair of gzipped fastq files
    and parses reads in each chunk into SequenceRecord objects

    :param fastq1: read1 compressed fastq
    :param fastq2: read2 compressed fastq
    :param buffer_size:  4*1024**2
    :return: list of tuple (SequenceRecord, SequenceRecord)
    """
    resolve_input = []
    with xopen(fastq1, "rb") as r1:
        with xopen(fastq2, "rb") as r2:
            chunk = 0
            total_reads = 0
            for c1, c2 in dnaio.read_paired_chunks(r1,r2, buffer_size):
                chunk += 1
                r1_in_chunk = parse_chunk(bytes(c1[:]))
                r2_in_chunk = parse_chunk(bytes(c2[:]))
                for read1, read2 in zip(r1_in_chunk, r2_in_chunk):
                    resolve_input.append((read1, read2))
                    total_reads += 1
            print(f"total number of chunks: {chunk}, reads: {total_reads}")
    return resolve_input




