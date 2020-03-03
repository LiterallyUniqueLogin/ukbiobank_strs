"""
List samples with a sex mismatch.

Those are defined to be participants whose inferred biological sex
differs from their reported sex.

See UKB nature paper supplementary information section S 3.6
https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0579-z/MediaObjects/41586_2018_579_MOESM1_ESM.pdf

This list heavily overlaps, but is not the same as, the sex aneuploidy list.
"""

import list_from_sample_qc as sqc


def sex_mismatch(row):
    return row[9] != row[10]


sqc.produce_list_from_sample_qc('sex_mismatch', False, func=sex_mismatch)
