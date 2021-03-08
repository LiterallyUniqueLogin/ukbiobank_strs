
import numpy as np

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
    assert len(set(a[:, 0]).intersection(b[:,0])) > 1000
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

