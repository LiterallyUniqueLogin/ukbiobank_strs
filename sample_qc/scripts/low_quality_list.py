"""
List low quality samples to exclude.

I predict most analyses will exclude these individuals
low quality and should be excluded form most analyses
See UKB nature paper supplementary information
https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0579-z/MediaObjects/41586_2018_579_MOESM1_ESM.pdf
This includes (supplementary section number)
* calls are overly missing (S 3.5)
* outlier in number of heterozygous loci (S 3.5)
* seems excessively related (S 3.7.1)
"""

import list_from_sample_qc as sqc


sqc.produce_list_from_sample_qc('low_quality', False, column=21)
