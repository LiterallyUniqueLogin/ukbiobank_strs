import os
import time
import vcf

ukb = os.environ['UKB']
print("Opening file")
chr21 = vcf.Reader(filename=ukb+"/str_imputed/runs/first_pass/vcfs/chr21.vcf.gz")
print("File opened")
nlines = 0
start = time.time()
for line in chr21:
    nlines += 1
    time_per_line = (time.time() - start)/nlines
    print("End of line {}. Time per line {:0.2f}".format(nlines,
                                                        time_per_line))
