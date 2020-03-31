import os
import time
import cyvcf2

ukb = os.environ['UKB']
print("Opening file")
vcf = cyvcf2.VCF(f"{ukb}/str_imputed/runs/first_pass/vcfs/chr21.vcf.gz")
print("File opened")

start = time.time()
nlines = 0
for variant in vcf("21:15449726-15454050"):
    print("POS:", variant.POS, "shape:", variant.format("AP1").shape)
    nlines += 1
    if nlines % 10 == 0:
        time_per_line = (time.time() - start) / nlines
        print("End of line {}. Time per line {:0.2f}".format(nlines,
                                                            time_per_line))
