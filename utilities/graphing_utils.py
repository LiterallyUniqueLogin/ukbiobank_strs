import math
from typing import Optional

import bokeh.models.plots
import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np

# matplotlib

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

# bokeh

def resize_font(obj, fontattr, ratio):
    font = obj.__getattribute__(fontattr)
    font_size, font_unit = font[:-2], font[-2:]
    obj.__setattr__(fontattr, str(math.floor(int(font_size)*ratio)) + font_unit)

def resize_int(obj, intattr, ratio):
    val = obj.__getattribute__(intattr)
    obj.__setattr__(intattr, math.floor(val*ratio))

def resize(fig, ratio, legend=True):
    resize_int(fig, 'width', ratio)
    resize_int(fig, 'height', ratio)
    resize_font(fig.title, 'text_font_size', ratio)
    for axis in fig.axis:
        resize_font(axis, 'axis_label_text_font_size', ratio)
        resize_font(axis, 'major_label_text_font_size', ratio)
        resize_int(axis, 'axis_label_standoff', ratio)
        resize_int(axis, 'axis_line_width', ratio)
        resize_int(axis, 'major_tick_line_width', ratio)
        resize_int(axis, 'major_tick_in', ratio)
        resize_int(axis, 'major_tick_out', ratio)
        resize_int(axis, 'major_label_standoff', ratio)
        resize_int(axis, 'minor_tick_line_width', ratio)
        resize_int(axis, 'minor_tick_in', ratio)
        resize_int(axis, 'minor_tick_out', ratio)
        if isinstance(axis, bokeh.models.CategoricalAxis):
            resize_int(axis, 'separator_line_width', ratio)
            resize_int(axis, 'separator_line_dash_offset', ratio)
            resize_int(axis, 'group_text_line_height', ratio)
            resize_font(axis, 'group_text_font_size', ratio)
            resize_int(axis, 'subgroup_text_line_height', ratio)
            resize_font(axis, 'subgroup_text_font_size', ratio)
    if legend:
        resize_font(fig.legend, 'label_text_font_size', ratio)
        resize_int(fig.legend, 'title_standoff', ratio)
        resize_int(fig.legend, 'border_line_width', ratio)
        resize_int(fig.legend, 'label_standoff', ratio)
        resize_int(fig.legend, 'label_width', ratio)
        resize_int(fig.legend, 'label_height', ratio)
        resize_int(fig.legend, 'glyph_width', ratio)
        resize_int(fig.legend, 'glyph_height', ratio)
        resize_int(fig.legend, 'padding', ratio)
        resize_int(fig.legend, 'spacing', ratio)
        resize_int(fig.legend, 'margin', ratio)
    for place in bokeh.core.enums.Place:
        for obj in getattr(fig, place):
            if isinstance(obj, bokeh.models.Title):
                resize_font(obj, 'text_font_size', ratio)
            if isinstance(obj, bokeh.models.ColorBar):
                resize_font(obj, 'title_text_font_size', ratio)
                resize_int(obj, 'title_text_line_height', ratio)
                resize_int(obj, 'title_standoff', ratio)
                resize_font(obj, 'major_label_text_font_size', ratio)
                resize_int(obj, 'major_label_text_line_height', ratio)
                resize_int(obj, 'label_standoff', ratio)
                resize_int(obj, 'major_tick_line_width', ratio)
                resize_int(obj, 'major_tick_in', ratio)
                resize_int(obj, 'major_tick_out', ratio)
                resize_int(obj, 'minor_tick_line_width', ratio)
                resize_int(obj, 'minor_tick_in', ratio)
                resize_int(obj, 'minor_tick_out', ratio)
                resize_int(obj, 'bar_line_width', ratio)
                resize_int(obj, 'border_line_width', ratio)
                resize_int(obj, 'margin', ratio)
                resize_int(obj, 'padding', ratio)
                if isinstance(obj.height, int):
                    resize_int(obj, 'height', ratio)
                if isinstance(obj.width, int):
                    resize_int(obj, 'width', ratio)

