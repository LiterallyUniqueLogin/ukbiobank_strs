import pandas as pd

import python_array_utils as utils

def get_relation(str_start, str_end, feature_start, feature_end, feature_direction):
    assert str_start <= str_end and feature_start <= feature_end
    assert feature_direction in {'-', '+'}
    if feature_start <= str_start and str_end <= feature_end:
        return 'inside'
    if str_start <= feature_start and feature_end <= str_end:
        raise ValueError(f'STR spans feature {str_start} {str_end} {feature_start} {feature_end}')
    if str_start <= feature_start:
        if feature_direction == '+':
            return 'upstream'
        else:
            return 'downstream'
    else:
        if feature_direction == '+':
            return 'downstream'
        else:
            return 'upstream'

def distance(start1, end1, start2, end2):
    if start1 <= start2 <= end1 or start2 <= start1 <= end2:
        # overlapping
        return '0'
    return str(min(abs(start1 - end2), abs(end1 - start2)))

def get_merged_annotations(signals, annotation_files, distance=False, bp_overlap=False, how='inner'):
    dfs = []
    names=[
        'chrom',
        'STR_pos',
        'STR_ID',
        'annotation_feature_type',
        'annotation_pos',
        'annotation_end_pos',
        'annotation_strand',
        'annotation_phase',
        'annotation_info'
    ]
    if distance:
        names.append('annotation_distance')
    if bp_overlap:
        names.append('bp_overlap')
    for fname in annotation_files:
        dfs.append(pd.read_csv(
            fname,
            delimiter='\t',
            names=names,
            dtype=utils.get_dtypes(fname, colnames = names),
            index_col=False
        ))
    full_df = pd.concat(dfs)
    merge = signals.merge(
        full_df,
        left_on=['chrom', 'pos'],
        right_on=['chrom', 'STR_pos'],
        how=how
    )
    return merge

def get_gff_kvp(kvp_string, key):
    if kvp_string[:len(key)] == key:
        splits = kvp_string.split(f'{key}=', maxsplit=1)
    else:
        splits = kvp_string.split(f';{key}=', maxsplit=1)

    if len(splits) == 1:
        return 'missing from gencode'
    elif len(splits) > 2:
        return 'too many splits'
    else:
        return str(splits[1].split(';', maxsplit=1)[0])

