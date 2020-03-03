"""
List putative sex aneuploidy samples.

See UKB nature paper supplementary information section S 3.6
https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0579-z/MediaObjects/41586_2018_579_MOESM1_ESM.pdf

This list heavily overlaps, but is not the same as, the sex mismatch list.
"""

import list_from_sample_qc as sqc


sqc.produce_list_from_sample_qc('sex_aneuploidy', False, column=19)
