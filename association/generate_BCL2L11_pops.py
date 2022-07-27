import numpy as np
import cyvcf2 ;

vcf = cyvcf2.VCF('str_imputed/runs/first_pass/vcfs/annotated_strs/chr2.vcf.gz')
var = next(vcf('2:111878544'))

len_alleles = [len(var.REF)] + [len(alt) for alt in var.ALT]
n_samples = len(vcf.samples)
dosage_gts = {
    _len: np.zeros((n_samples, 2)) for _len in np.unique(len_alleles)
}
for p in (1, 2):
    # todo genotype dosages
    ap = var.format(f'AP{p}')
    dosage_gts[len_alleles[0]][:, (p-1)] += \
            np.maximum(0, 1 - np.sum(ap, axis=1))
    for i in range(ap.shape[1]):
        dosage_gts[len_alleles[i+1]][:, (p-1)] += ap[:, i]
dosage_gts_arr = np.stack([dosage_gts[len_] for len_ in dosage_gts])
best_guess_coords = np.argmax(dosage_gts_arr, axis=0)
best_guess_lens = np.unique(len_alleles)[best_guess_coords]
sample_subsets = {int(np.floor(len_/3)): np.any(best_guess_lens == len_, axis=1) for len_ in np.unique(len_alleles)}

for len_, subset in sample_subsets.items():
    with open(f'sample_qc/subpops/BCL2L11_len_{len_}.txt', 'w') as f:
        f.write('ID\n')
        samples = np.array(vcf.samples)[subset]
        samples = [sample[:7] for sample in samples]
        f.write('\n'.join(np.array(samples)))
