import dnaio_io


read1 = "../data/CEG9330132-19-01_S12_L001_R1_001.fastq.gz"
read2 = "../data/CEG9330132-19-01_S12_L001_R2_001.fastq.gz"

# Testing read-write capabilities
dnaio_io.read_single(read1)
dnaio_io.read_pair(read1, read2)

dnaio_io.split(read1)
#TODO: understand SequenceRecord object instantiation
dnaio_io.combine(read1,read2)

## Testing chunked reading capabilities


