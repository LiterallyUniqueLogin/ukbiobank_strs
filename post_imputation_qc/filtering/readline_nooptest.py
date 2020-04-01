import gzip
import os
import Bio.bgzf
import sys
import time

ukb = os.environ['UKB']
print("Opening file")
with Bio.bgzf.BgzfReader(f"{ukb}/str_imputed/runs/first_pass/vcfs/chr21.vcf.gz") as chr21:
    print("File opened")
    nlines = 0
    for line in chr21:
        if line[0] != int.from_bytes(b'#', "little"):
            break
    start = time.time()
    for line in chr21:
        nlines += 1
        if nlines % 10 == 0:
            time_per_line = (time.time() - start)/nlines
            print("End of line {}. Time per line {:0.4f}".format(nlines,
                                                                 time_per_line))
