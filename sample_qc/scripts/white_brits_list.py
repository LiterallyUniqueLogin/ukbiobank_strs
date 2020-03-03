"""
List white samples with British ancestry.

See UKB nature paper supplementary information section S 3.4
https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0579-z/MediaObjects/41586_2018_579_MOESM1_ESM.pdf
"""

import list_from_sample_qc as sqc


sqc.produce_list_from_sample_qc('white_brits', True, column=23)
