import os
import sys

ukb = os.environ['UKB']

sys.path.insert(0, f'{ukb}/../trtools/repo')

import trtools.utils.tr_harmonizer
