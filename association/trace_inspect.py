import tracemalloc
'''
pre_loop = tracemalloc.Snapshot.load(
    f'{prefix}/21:48108278-48108278_loop_begin'
)
loop_start = tracemalloc.Snapshot.load(
    f'{prefix}/21:48108278-48108278_loop_48108278_start'
)
'''

print('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
prefix = 'runs/test_snps/profiling'
print("Loading first snapshot ... ", flush=True)
preOLS = tracemalloc.Snapshot.load(
    f'{prefix}/21:48108278-48108278_loop_48108278_preOLS'
)
print("Loading second snapshot ...", flush=True)
postOLS = tracemalloc.Snapshot.load(
    f'{prefix}/21:48108278-48108278_loop_48108278_postOLS'
)
print("Computing comparison ...", flush=True)
for stat in postOLS.compare_to(preOLS, 'lineno'):
    if abs(stat.size_diff) < 1024**2:
        continue
    print(stat)

for tb in postOLS.compare_to(preOLS, 'traceback'):
    if abs(tb.size_diff) < 1024**2:
        continue
    print('\n---------\n')
    print(f'{tb.size_diff // 1024**2}MiB diff, {tb.count_diff} change in num objects')
    for line in tb.traceback.format():
        print(line)
