from typing import Optional

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np

def PlotCDF(
        data: np.ndarray,
        fname: str,
        dist_name: str,
        measurement: str,
        data_len: Optional[int] = None,
        reverse: bool = False):
    spacing = 5e-3
    fig = plt.figure()
    ax = fig.add_subplot(111)
    measurable_data = data[~np.isnan(data) & ~np.isinf(data)]
    ax.set_xlim(np.min(measurable_data) - spacing, np.max(measurable_data) + spacing)
    assert len(data.shape) == 1
    _BetterCDF(data, ax, data_len=data_len, reverse=reverse)
    ax.set_xlabel(measurement, size=15)
    if not reverse:
        ax.set_ylabel(f"% of {dist_name} with at least this much \n {measurement}", size=15)
    else:
        ax.set_ylabel(f"% of {dist_name} with no more than this \n {measurement}", size=15)
    plt.tight_layout()
    fig.savefig(fname)
    plt.close()

# original idea https://stackoverflow.com/a/39729964/2966505
def _BetterCDF(data: np.ndarray,
               ax: matplotlib.axes.Axes,
               data_len: Optional[int] = None,
               reverse: bool = False):
    # assumes that axes are already set to (min, max)
    data = np.sort(data)
    x_axis_min, x_axis_max = ax.get_xlim()
    if not reverse:
        xfirst = x_axis_min
        xlast = x_axis_max
    else:
        data = data[::-1]
        xfirst = x_axis_max
        xlast = x_axis_min

    if data_len is None:
        data_len = len(data)
    n_points = len(data)
    has_quality_1_point = data[-1] == 1
    if has_quality_1_point:
        # don't print a drop off if the last data point(s)
        # have quality 1
        n_ones = sum(data == data[-1])
        data = np.hstack((
            [xfirst],
            data[0:(len(data) - n_ones)],
            [xlast]
        ))
        ys = np.hstack((
            [1],
            np.arange(n_points - 1, n_ones - 1, -1) / data_len,
            [n_ones / data_len]
        ))
    else:
        data = np.hstack((
            [xfirst],
            data,
            [xlast]
        ))
        ys = np.hstack((
            [1],
            np.arange(n_points - 1, -1, -1) / data_len,
            [0]
        ))
    ax.step(data, ys, where='post')
