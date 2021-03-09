import os

import cyvcf2
import numpy as np

ukb = os.environ['UKB']

vcf = cyvcf2.VCF(
    f'{ukb}/str_imputed/runs/first_pass/vcfs/annotated_strs/chr1.vcf.gz'
)
vcf_samples = np.array(
    list(sample.split("_")[0] for sample in vcf.samples)
)

bgen_samples = []
with open(f'{ukb}/microarray/ukb46122_hap_chr1_v2_s487314.sample') as samplefile:
    for num, line in enumerate(samplefile):
        if num <= 1:
            # skip first two lines
            continue
        bgen_samples.append(line.split()[0])
bgen_samples = np.array(bgen_samples)

assert np.all(vcf_samples == bgen_samples)

