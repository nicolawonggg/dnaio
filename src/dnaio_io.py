import dnaio


# Reading FQSTQ with dnaio

# dnaio.open returns a FastqReader object
# - interator of a SequenceRecord obj
def read_single(fastq):
    with dnaio.open(fastq) as reader:
        bp = 0
        for record in reader:
            bp += len(record)
    print(f"file read: {fastq}")
    print(f"The input file contains {bp/1E6:.1f} Mbp")


## by providing to reads, dnaio.open returns a TwoFilePairedEndReader
# - iterator of a tuple of 2 SequenceRecord objs
def read_pair(fastq1, fastq2):
    with dnaio.open(fastq1, fastq2) as reader:
        bp = 0
        for r1, r2 in reader:
            bp += len(r1) + len(r2)
        print(f"files read: {fastq1}, {fastq2}")
        print(f"The paired-end input contains {bp/1E6:.1f} Mbp")


# Writing FQSTQ with dnaio
# one file in, two files out
def split(fastq_in):
    with dnaio.open(fastq_in) as reader, \
            dnaio.open("../output/CEG9330132-19-01_S12_L001_R1_001_1.fastq.gz",
                       "../output/CEG9330132-19-01_S12_L001_R1_001_2.fastq.gz", mode="w") as writer:
        for record in reader:
            r1_1 = record[:int(len(record)/2)]
            r1_2 = record[-int(len(record)/2):]
            writer.write(r1_1, r1_2)

    ## checking the output FASTQs
    read_pair("../output/CEG9330132-19-01_S12_L001_R1_001_1.fastq.gz",
              "../output/CEG9330132-19-01_S12_L001_R1_001_2.fastq.gz")

# Writing FQSTQ with dnaio
# two files in, one file out
def combine(fastq_in1, fastq_in2):
    filename = "../output/combined.fastq.gz"
    with dnaio.open(fastq_in1, fastq_in2) as reader, \
            dnaio.open(filename, mode="w") as writer:
        sequence = ''
        qualities = ''
        for r1, r2 in reader:
            ##Some fake resolve logic
            #use first half of r1
            sequence.join(r1[:int(len(r1)/2)].sequence)
            qualities.join(r1[:int(len(r1)/2)].qualities)
            #use second half of r2
            sequence.join(r2[-int(len(r2)/2):].sequence)
            qualities.join(r2[:int(len(r2)/2)].qualities)
            # print(f"processed sequence full length:  {len(sequence)}")
            #create dummy resolved single-end read
        out = dnaio.SequenceRecord(name=r1.name.join("_comb"), sequence=sequence, qualities=qualities)
        writer.write(out)
        print(f"file saved: {filename}")

    ## checking the output FASTQs
    read_single("../output/combined.fastq.gz")
