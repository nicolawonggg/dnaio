import dnaio

#
## Reading FQSTQ with dnaio
#
# # dnaio.open returns a FastqReader object
# # - interator of a SequenceRecord obj
# with dnaio.open("../data/CEG9330132-19-01_S12_L001_R1_001.fastq.gz") as reader:
#     bp = 0
#     for record in reader:
#         bp += len(record)
# print(f"The input file contains {bp/1E6:.1f} Mbp")
#

# ## by providing to reads, dnaio.open returns a TwoFilePairedEndReader
# # - iterator of a tuple of 2 SequenceRecord objs
# with dnaio.open("../data/CEG9330132-19-01_S12_L001_R1_001.fastq.gz",
#                 "../data/CEG9330132-19-01_S12_L001_R2_001.fastq.gz") as reader:
#     bp = 0
#     for r1, r2 in reader:
#         bp += len(r1) + len(r2)
#     print(f"The paired-end input contains {bp/1E6:.1f} Mbp")
#
#
# # Writing FQSTQ with dnaio
# # one file in, two files out
# with dnaio.open("../data/CEG9330132-19-01_S12_L001_R1_001.fastq.gz") as reader, \
#         dnaio.open("../output/CEG9330132-19-01_S12_L001_R1_001_1.fastq.gz",
#                    "../output/CEG9330132-19-01_S12_L001_R1_001_2.fastq.gz", mode="w") as writer:
#     for record in reader:
#         r1_1 = record[:int(len(record)/2)]
#         r1_2 = record[-int(len(record)/2):]
#         writer.write(r1_1, r1_2)
#
# ## checking the output FASTQs
# with dnaio.open("../output/CEG9330132-19-01_S12_L001_R1_001_1.fastq.gz",
#                 "../output/CEG9330132-19-01_S12_L001_R1_001_2.fastq.gz") as reader:
#     bp1 = 0
#     bp2 = 0
#     for r1_1, r1_2 in reader:
#
#         bp1 += len(r1_1)
#         bp2 += len(r1_2)
#     print(f"{r1_1.name} contains {bp1/1E6:.1f} Mbp in total")
#     print(f"{r1_2.name} contains {bp2/1E6:.1f} Mbp in total")


# Writing FQSTQ with dnaio
# two files in, one file out
with dnaio.open("../data/CEG9330132-19-01_S12_L001_R1_001.fastq.gz",
                "../data/CEG9330132-19-01_S12_L001_R2_001.fastq.gz") as reader, \
        dnaio.open("../output/combined.fastq.gz", mode="w") as writer:
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
        print(f"processed sequence full length:  {len(sequence)}")
        #create dummy resolved single-end read
        out = dnaio.SequenceRecord(name=r1.name.join("_comb"), sequence=sequence, qualities=qualities)
        writer.write(out)

## checking the output FASTQs
with dnaio.open("../output/combined.fastq.gz") as reader:
    bp = 0
    for record in reader:
        bp += len(record)
    print(f"{record.name} contains {bp/1E6:.1f} Mbp in total")
