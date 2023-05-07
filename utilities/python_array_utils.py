import io
import sys

import numpy as np
import polars as pl

def merge_arrays(a, b):
    '''
    Return a left outer join b.
    Join is performed on the first column.

    Assume first column of each array is id, but not necessarily same order

    Parameters
    ----------
    a, b: np.ndarray
        2D arrays
    '''

    assert len(a.shape) == 2 and len(b.shape) == 2
    intersection_size = len(set(a[:, 0]).intersection(b[:,0]))
    if intersection_size <= 1000:
        print(f"Working with a sample intersection of only {intersection_size} samples. Is this intentional, or possibly a coding bug?", file=sys.stderr)
    assert len(set(a[:, 0])) == a.shape[0]
    assert len(set(b[:, 0])) == b.shape[0]

    b = b[np.isin(b[:, 0], a[:, 0])] # drop all elements of b that aren't in a
    matches = np.isin(a[:, 0], b[:, 0]) # a indicies that are in b

    a_sort = np.argsort(a[matches, 0])
    b_match_sorted = np.searchsorted(a[matches, 0], b[:, 0], sorter=a_sort)

    new_data = np.full((a.shape[0], b.shape[1] - 1), np.nan) # to hold the data from b
    new_data[matches, :] = b[np.argsort(b_match_sorted), 1:][np.argsort(a_sort), :]

    return np.concatenate((
        a,
        new_data
    ), axis=1)

# from https://stackoverflow.com/questions/52579601/convert-dataframe-to-a-rec-array-and-objects-to-strings
def df_to_recarray(df):
    '''
    Parameters
    ----------
    df: pandas dataframe
    '''
    names = df.columns
    arrays = [df[col].values for col in names]

    formats = [ array.dtype if array.dtype != 'O'
                else f'{array.astype(str).dtype}' for array in arrays ]

    rec_array = np.rec.fromarrays(
        arrays,
        dtype={'names': names, 'formats': formats}
    )

    return rec_array

def get_dtypes(fname, type_dict = {}, colnames = None):
    import pandas as pd
    if colnames is None:
        df = pd.read_csv(
            fname,
            header=0,
            delimiter='\t',
            encoding='UTF-8',
            nrows=1
        )
    else:
        df = pd.read_csv(
            fname,
            delimiter='\t',
            encoding='UTF-8',
            nrows=1,
            names=colnames
        )
    dtypes = dict(df.dtypes)
    for key, val in type_dict.items():
        dtypes[key] = val
    return dtypes

