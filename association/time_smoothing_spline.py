import time

import numpy as np
import numpy.random
import csaps

rng = numpy.random.default_rng(13)

n_iters = 5
print(f"Calculations for n_iters = {n_iters}")
n = 1000
while n < 5e5:
    total = 0
    for _ in range(n_iters):
        ys = rng.standard_normal(n)
        xs = np.sort(rng.standard_normal(n))
        #xs = np.sort(rng.choice(100, size=n))
        start = time.time()
        yi = csaps.csaps(xs, ys, range(100), smooth=0.8)
        total += time.time() - start
    print(f"Avg time for n={n}: {total/n_iters}")

    n = n*4

