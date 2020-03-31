import gzip
import os
import Bio.bgzf
import sys
import time

ukb = os.environ['UKB']
print("Opening file")
with Bio.bgzf.BgzfReader(f"{ukb}/str_imputed/runs/first_pass/vcfs/chr21.vcf.gz", 'b') as chr21:
    print("File opened")
    nlines = 0
    for line in chr21:
        if line[0] != int.from_bytes(b'#', "little"):
            break
    start = time.time()
    while True:
        byte = chr21.read(1)
        if byte == b'\n':
            nlines += 1
            time_per_line = (time.time() - start)/nlines
            print("End of line {}. Time per line {:0.2f}".format(nlines,
                                                                time_per_line))
            byte = chr21.read(1)
