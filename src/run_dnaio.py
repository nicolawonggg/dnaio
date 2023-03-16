import random
import dnaio_chunk
import dnaio_io
from multiprocessing import Pool


read1 = "../data/CEG9330132-19-01_S12_L001_R1_001.fastq.gz"
read2 = "../data/CEG9330132-19-01_S12_L001_R2_001.fastq.gz"
cpu_count = 4
buffer_size =  4*1024**2 #maximum length for a chunk


### mimicking discard / resolve logic
def _resolve(read1, read2):
    dis_r1 = [];
    dis_r2 = [];
    res = []
    if random.random() <= 0.3:
        dis_r1.append(read1)
        dis_r2.append(read2)
    else:
        res.append(read1)
    return dis_r1, dis_r2, res

if __name__ == '__main__':

## Testing read-write capabilities
    # dnaio_io.read_single(read1)
    # dnaio_io.read_pair(read1, read2)
    #
    # dnaio_io.split(read1)
    # dnaio_io.combine(read1,read2)

## Testing chunked reading capabilities
    # read1_chunks = dnaio_chunk.chunk_single(read1,128)
    # dnaio_io.write_fastqs(read1_chunks, "../output/read1_chunked.fastq.gz")
    # dnaio_chunk.chunk_single(read1, 4)

## Testing chunk read pairs, then output discard pairs and resolved reads
    resolve_input = dnaio_chunk.chunk_pair(read1, read2,buffer_size)
    discard_r1 = []
    discard_r2 = []
    resolved = []
    # create a pool with the number of processing matching the available cpus
    with Pool(cpu_count) as pool:
        for dis_r1, dis_r2, res in pool.starmap(
            _resolve, resolve_input
        ):
            discard_r1.extend(dis_r1)
            discard_r2.extend(dis_r2)
            resolved.extend(res)

    print(f"discarded r1 :{len(discard_r1)},  r2 :{len(discard_r2)}")
    print(f"resolved :{len(resolved)}")

    if dnaio_io.check_mates(discard_r1, discard_r2):
        print("passed check: discard files pair align with each other")

    dnaio_io.write_fastqs(discard_r1, "../output/test/discarded_read1.fastq.gz")
    dnaio_io.write_fastqs(discard_r2, "../output/test/discarded_read2.fastq.gz")
    dnaio_io.write_fastqs(resolved, "../output/test/resolved.fastq.gz")

    # dnaio_io.read_single("../output/test/resolved.fastq.gz")