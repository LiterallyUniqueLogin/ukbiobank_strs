import bgen_reader

def get_alleles(variants, pos):
    same_pos_vars = variants[variants['pos'] == pos]
    allele_sets = []
    for _, var in same_pos_vars.iterrows():
        allele_sets.append(set(var['allele_ids'].split(',')))
    return allele_sets

for chrom in range(1, 23):
    print("Working on chrom ", chrom, flush=True)
    imputed_variants = \
            bgen_reader.read_bgen(f"ukb_imp_chr{chrom}_v3.bgen")['variants'].compute()
    array_variants = \
            bgen_reader.read_bgen(f"../microarray/ukb_hap_chr{chrom}_v2.bgen")['variants'].compute()
    assert (imputed_variants.nalleles == 2).all()
    assert (array_variants.nalleles == 2).all()

    # assert distinct rsids
    imp_rsids = set(imputed_variants['rsid'].values)
    array_rsids = set(array_variants['rsid'].values)
    for rsid in array_rsids:
        if rsid in imp_rsids:
            print("Shared rsid ", rsid, flush=True)

    # assert distinct variants (naively)
    imp_poses = set(imputed_variants['pos'].values)
    array_poses = set(array_variants['pos'].values)
    for pos in array_poses:
        if pos in imp_poses:
            imp_allele_sets = get_alleles(imputed_variants, pos)
            for allele_set in get_alleles(array_variants, pos):
                if allele_set in imp_allele_sets:
                    print(f"Shared variant {pos} {allele_set}", flush=True)
