import os  


def produce_list_from_sample_qc(run_name, list_name, keep, column):
    """
    In/exclude participants marked by col <column> in the sample qc file.

    Parameters
    ----------
    run_name
        the name of the run to output the file to
    list_name
        the name of the output file (w/o extension)
    keep
        whether the participants identified by this column should be kept
        or excluded
    column
        the column of the sample qc file in the misc/EGA directory to look at

    Returns
    -------
    None
        Writes out a list of samples to
        $UKB/sample_qc/runs/run_name/keep/list_name.sample
        (or with remove instead of keep depending)

    """
    if keep:
        keep_or_remove = "keep"
    else:
        keep_or_remove = "remove"
    ukb = os.environ["UKB"]

    sample_floc = f"{ukb}/microarray/ukb46122_hap_chr1_v2_s487314.sample"
    with open(sample_floc) as sample_file:
        sample_header = sample_file.readline()

    sqc_floc = f"{ukb}/misc_data/EGA/ukb_sqc_v2.txt"
    fam_floc = f"{ukb}/microarray/ukb46122_cal_chr1_v2_s488282.fam"
    out_floc = (f"{ukb}/sample_qc/runs/{run_name}/"
                f"{keep_or_remove}/{list_name}.sample")
    with open(sqc_floc) as sqc_file, \
            open(fam_floc) as fam_file, \
            open(out_floc, 'w') as out_file:
        out_file.write(sample_header)

        # Assumes already validated that csv and fam are the same length
        for line in sqc_file:
            fam_tokens = fam_file.readline().split()
            if bool(int(line.split()[column])):
                out_file.write("{} {} {} {}\n".format(
                    fam_tokens[0],
                    fam_tokens[1],
                    fam_tokens[2],
                    fam_tokens[4]))

