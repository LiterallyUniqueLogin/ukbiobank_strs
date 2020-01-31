import os
import Bio.bgzf
import gzip

import random
import string

ranstr = ""
for _ in range(int(5e5)):
	ranstr += random.choice(string.ascii_letters)

ranstr = ranstr.encode()

tmpdir = os.environ['TMPDIR']


with open(f"{tmpdir}/ran.txt", 'wb') as txt:
	txt.write(ranstr)

with gzip.open(f"{tmpdir}/ran.txt.gz", 'w') as gz:
	gz.write(ranstr)

with bio.bgzf.open(f"{tmpdir}/ran.txt.bgz.gz", 'w') as bgz:
	bgz.write(ranstr)

"""
to compare, cd into tmpdir and run
gunzip -c ran.txt.gz > rangz.txt
cmp ran.txt rangz.txt
bgzip -cd ran.txt.bgz.gz > ranbgz.txt
cmp ranbgz.txt ran.txt
"""
pass
