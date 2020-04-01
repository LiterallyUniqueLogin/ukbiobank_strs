import gzip
import os
import sys
import gzip
import time

ukb = os.environ['UKB']
print("Opening file")
with gzip.open(f"{ukb}/str_imputed/runs/first_pass/vcfs/chr21.vcf.gz") as chr21:
    print("File opened")
    nlines = 0
    for line in chr21:
        if line[0] != int.from_bytes(b'#', "little"):
            break
    start = time.time()
    for line in chr21:
        for char in line:
            if char == b'\n':
                print("Hello!")
        nlines += 1
        if nlines % 10 == 0:
            time_per_line = (time.time() - start)/nlines
            print("End of line {}. Time per line {:0.2f}".format(nlines,
                                                                time_per_line))
