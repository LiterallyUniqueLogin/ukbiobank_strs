import bgen_reader

for chrom in {10}:#range(1, 23):
    bgen = bgen_reader.read_bgen(f"ukb_imp_chr{chrom}_v3.bgen")
    n_bgen_vars = bgen['variants'].shape[0].compute()
    with open(f"ukb_mfi_chr{chrom}_v3.txt") as mfi:
        assert len(mfi.readlines()) == n_bgen_vars
        print(f"Done with chr{chrom}")
