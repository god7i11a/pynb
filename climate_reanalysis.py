#!/usr/bin/env python
# coding: utf-8

# # Resources
# # https://climatereanalyzer.org/clim/sst_daily/
# # https://www.pyngl.ucar.edu/
# # https://www.ncl.ucar.edu/External/
# # https://climatedataguide.ucar.edu/
#
# # why $5\sigma$ ???:    # https://home.cern/resources/faqs/five-sigma
#
# # TODOs:
# ## explore other ways of quantifying variations besides daily distribution around the 20 year mean.

import datetime
import ipywidgets as widgets
from numpy import array, linspace, vstack, arange
from pandas import DataFrame, concat, read_json, read_excel
from matplotlib import pyplot as plt, colors

# "text.latex.preamble": r"\usepackage{amsmath} \usepackage{amssymb} \usepackage[utf8]{inputenc}"
plt.rcParams.update({"text.usetex": True})

CID = None

datasetsD = {'T2': {'URL': 'https://climatereanalyzer.org/clim/t2_daily/json',
                    'end_offset': -4,
                    'mean_base': -2,
                    'fileN': 'era5_world_t2_day.json',
                    'data': None,
                    'name': 'T2',
                    'pos': {'loc': 'upper left', 'ncols': 3,
                            'sigma_caption_voffset': .125,
                            'caption_height': 11.1, 'ylim': [11, 18]},
                    'pos_jump': {'loc': 'upper left', 'ncols': 4,
                                 'sigma_caption_voffset': .125,
                                 'caption_height': 11.15, 'ylim': [11, 18]},
                    'jump_years': {1976: 'D', 1977: 'x', 1978: 's', 1979: '^', 1980: 'D', 1981: 'x', 1982: 's',
                                   1983: 's', 2022: '*', 2023: 'D', 2024: 'p', 2025: 'v'},
                    'jump_title': 'T2 jumps of 1976-1983 and 2022-2025',
                    'hl_years': {2022: 's', 2023: 'x', 2024: 'D', 2025: 'v'}
                    },
             'SST': {'URL': 'https://climatereanalyzer.org/clim/sst_daily/json_2clim',
                     'end_offset': -3,
                     'mean_base': -2,
                     'fileN': 'oisst2.1_world2_sst_day.json',
                     'data': None,
                     'name': 'SST',
                     'pos': {'loc': 'upper right', 'ncols': 4,
                             'sigma_caption_voffset': .05,
                             'caption_height': 19.54, 'ylim': [19.5, 21.5]},
                     'pos_jump': {'loc': 'upper right', 'ncols': 3,
                                  'sigma_caption_voffset': .05,
                                  'caption_height': 19.54, 'ylim': [19.5, 21.5]},
                     'jump_years': {1996: 'D', 1997: 'x', 1998: 's', 1999: '^',
                                    2022: '*', 2023: 'p', 2024: 'v', 2025: '<'},
                     'jump_title': 'SST jumps of 1996-1999 and 2022-2025',
                     'hl_years': {2009: 'o', 2010: 's', 2011: '^', 2012: 'v', 2013: '*', 2014: 'D', 2015: '<', 2016: '>',
                                  2023: 'x', 2024: 'D', 2025: 'v'}
                     },
             'SST_old': {'URL': 'https://climatereanalyzer.org/clim/sst_daily/json',
                         'end_offset': -4,
                         'fileN': 'oisst2.1_world2_sst_day.json',
                         'caption_height': 19.6,
                         'data': None
                         }
             }

markersL = ('1', '2', '3', '4', '+', 'x', '.', 'o', 'v', '<', '>',
            '8', 's', 'p', '*', 'h', 'H', 'd', 'D', 'P', 'X',
            4, 5, 6, 7, 8, 9, 10, 11)


class Dataset:

    def __init__(self, **kwD):
        self.init(kwD)

    def init(self, kwD):
        for kw, val in kwD.items():
            setattr(self, kw, val)


T2 = Dataset(**datasetsD['T2'])
SST = Dataset(**datasetsD['SST'])
SST_old = Dataset(**datasetsD['SST_old'])

# # get the latest data
#
# # global average   https://climatereanalyzer.org/clim/t2_daily/json/era5_world_t2_day.json
# # gridded data https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/


def update_data(ds):
    path = f"{ds.URL}/{ds.fileN}"
    df = read_json(path_or_buf=path)
    df.set_index("name", inplace=True)
    df = DataFrame(df["data"].to_list(), columns=list(range(1, 367)), index=df.index)
    df.to_excel(f'{ds.name}_raw_data.xlsx')


def get_data(ds):
    df = read_excel(f'{ds.name}_raw_data.xlsx')
    df.set_index('name', inplace=True, drop=True)
    ds.last_year = df.index.to_list()[ds.end_offset]
    ds.first_year = df.index.to_list()[0]
    ds.day_of_year = df.loc[ds.last_year, :].dropna().index.to_list()[-1]
    ds.last_date = (datetime.datetime(int(ds.last_year), 1, 1) + datetime.timedelta(ds.day_of_year - 1)).strftime('%m/%d/%Y')
    ds.yearL = range(int(ds.first_year), int(ds.last_year)+1)
    ds.end_data_period = df.index.to_list()[ds.end_offset+1:]
    ds.mean_baseL = df.index.to_list()[ds.mean_base].split('-')

    print(f'{ds.day_of_year=}')
    print(f'latest data point is {ds.last_date}')
    print(f'last year = {ds.last_year}')
    print(f'first year = {ds.first_year}')

    # # compute mean and standard deviations and display raw data

    cdf = df.T.loc[:, ds.first_year:]
    cdf['mean'] = cdf.loc[:, ds.mean_baseL[0]:ds.mean_baseL[1]].mean(axis=1, skipna=True)
    cdf['std'] = cdf.loc[:, ds.mean_baseL[0]:ds.mean_baseL[1]].std(axis=1, skipna=True)
    cdf['finalmean'] = cdf.loc[:, '2011':ds.last_year].mean(axis=1, skipna=True)
    cdf['finalstd'] = cdf.loc[:, '2011':ds.last_year].std(axis=1, skipna=True)
    cdf['mean+5std'] = cdf['mean']+5*cdf['std']
    cdf['mean+5finalstd'] = cdf['mean']+5*cdf['finalstd']
    cdf['mean+2std'] = cdf['mean']+2*cdf['std']
    cdf['mean-2std'] = cdf['mean']-2*cdf['std']
    ds.data = cdf


def display_latest_year_data(ds):

    cdf = ds.data

    _, ax_last = plt.subplots(1, 1, figsize=(10, 6), clear=True, num=1)
    ax_last.plot(array(cdf.index.to_list()), cdf[f'{ds.yearL[-1]}'], lw=1)
    plt.xlim([0, ds.day_of_year])
    plt.xlabel('day_of_year')
    plt.ylabel(r'$T\ \  (^{\circ} C)$', fontsize='xx-large')
    plt.title(f'latest year of data = {ds.last_year}')
    plt.savefig(f'{ds.name}_last_year.png')


# # Colormaps
# # https://colorcet.holoviz.org/
# # https://cmasher.readthedocs.io/user/usage.html#accessing-colormaps
def _display_color_chart(cmaps, cmap_types):
    # ## color schemes
    gradient = linspace(0, 1, 256)
    gradient = vstack((gradient, gradient))

    # Create figure and adjust figure height to number of colormaps

    ncols = len(cmap_types)
    nrows = max([len(arr) for arr in cmaps.values()])
    cbarfig, axsL = plt.subplots(nrows=nrows, ncols=ncols, sharey=True, sharex=True, figsize=(16, 16), num=20, clear=True)
    cbarfig.subplots_adjust(hspace=0.35, wspace=0.05)

    for axrow in axsL:
        for ax in axrow:
            # Turn off *all* ticks & spines, not just the ones with colormaps.
            ax.set_axis_off()

    for n, key in zip(range(ncols), cmaps.keys()):
        axs = [[j, axsL[j][n]] for j in range(nrows)]
        for (j, ax), cmap_name in zip(axs, cmaps[key]):
            if j == 0:
                ax.text(.5, 1.05, cmap_types[n], ha='center', fontsize=12, transform=ax.transAxes)
            ax.imshow(gradient, aspect='auto', cmap=cmap_name)
            ax.text(.5, -0.3, cmap_name, ha='center', fontsize=12, transform=ax.transAxes)
            # turn on just the box outline
            ax.set_axis_on()
            ax.set_xticks([])
            ax.set_yticks([])


def construct_colormap(show_chart=False):

    cmaps = {'Cyclic': ['twilight', 'twilight_shifted', 'hsv'],
             'Perceptually Uniform Seq': ['viridis', 'plasma', 'inferno', 'magma', 'cividis'],
             'Seq': ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                     'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                     'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'],
             'Seq2': ['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
                      'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                      'hot', 'afmhot', 'gist_heat', 'copper'],
             'Divirgent': ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                           'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'],
             'Qualitative': ['Pastel1', 'Pastel2', 'Paired', 'Accent',
                             'Dark2', 'Set1', 'Set2', 'Set3',
                             'tab10', 'tab20', 'tab20b', 'tab20c'],
             'Misc': ['flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
                      'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
                      'gist_rainbow', 'rainbow', 'jet', 'turbo', 'nipy_spectral',
                      'gist_ncar']}

    cmap_types = list(cmaps.keys())
    cmap_names = [cname for ctype in cmap_types for cname in cmaps[ctype]]

    if show_chart:
        _display_color_chart(cmaps, cmap_types)

    w = widgets.Dropdown(
        options=cmap_names,
        value='jet',
        description='ColorMap Choice:',
    )

    return w


# ## graph
# isDrawn1 = 0
# isDrawn0 = 0


def display_main_plot(ds, w):
    global PICK, CID
    SHOWSIG = True
    ACCENTS = True
    PICK = False
    if (not PICK) and CID:
        print('disconnecting')
        CID.disconnect()

    def on_pick(event):
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata
        ax.text(x, y, event.artist.get_label())

    cmapStr = w.value

    cdf = ds.data

    if SHOWSIG:    # 15,11
        fig = plt.figure(figsize=(15, 11), clear=True, num=2, layout='constrained')
        (ax, ax1) = fig.subplots(2, 1, sharex=True, height_ratios=[.8, .2])
    else:
        fig, ax = plt.subplots(1, 1, figsize=(16, 16), num=2, clear=True)

    if PICK:
        CID = fig.canvas.mpl_connect('pick_event', on_pick)

    ax.plot(cdf.index, cdf['mean'], 'b', lw=4, label=f'{ds.mean_baseL[0]}-{ds.mean_baseL[1]} mean')
    ax.plot(cdf.index, cdf['mean+2std'], 'k', lw=4, label=rf'{ds.mean_baseL[0]}-{ds.mean_baseL[1]} mean $\pm 2\sigma$')
    ax.plot(cdf.index, cdf['mean-2std'], 'k', lw=4)
    ax.plot(cdf.index, cdf['mean+5std'], 'r', lw=4, label=rf'{ds.mean_baseL[0]}-{ds.mean_baseL[1]} mean +$5\sigma$')

    num_colors = len(ds.yearL)
    cmap = plt.get_cmap(cmapStr, num_colors)

    _skip = 10
    ms = 5

    labels = None
    for i, _year in enumerate(ds.yearL):
        if PICK:
            labels = _year

        ax.plot(cdf.index, cdf.loc[:, f'{_year}'], c=cmap(i), lw=1, label=labels, picker=PICK, pickradius=5)
        if ACCENTS:
            if _year in ds.hl_years.keys():
                day_num = cdf.index.to_list()[0::_skip]
                data = cdf[f'{_year}'].iloc[0::_skip]
                ax.plot(day_num, data, ds.hl_years[_year], c=cmap(i), ms=ms, label=f'{_year}')

    norm = colors.Normalize(vmin=int(ds.first_year), vmax=ds.yearL[-1])
    ds.sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # check that the color chosen in is the middle of the color segment
    plt.colorbar(ds.sm, ax=ax, ticks=linspace(int(ds.first_year), ds.yearL[-1], num_colors),
                 boundaries=arange(int(ds.first_year)-.5, ds.yearL[-1]+1, 1))

    posD = ds.pos
    ax.set_ylabel(r'$T\ \  (^{\circ} C)$', fontsize='xx-large')
    ax.set_title(f'{ds.first_year}-{ds.last_year} daily {ds.name}', fontsize='xx-large')
    ax.text(275, posD['caption_height'],  f'\\it{{data current to {ds.last_date}}}', fontsize='xx-large')
    ax.set_xlim(0, 365)
    ax.set_ylim(*posD['ylim'])
    if not PICK:
        _legend = ax.legend(ncols=posD['ncols'], loc=posD['loc'])

    textStr = rf'$\sigma = \sigma({ds.mean_baseL[0]}-{ds.mean_baseL[1]})$'

    if SHOWSIG:
        ax1.plot(cdf.index, cdf['std'], 'b.', ms=4)
        ax1.set_ylabel(r'$\sigma\ \  (^{\circ} C)$', fontsize='xx-large')
        ax1.set_xlabel(r'day \#', fontsize='xx-large')
        ax1.set_title(r'daily $\sigma$', fontsize='xx-large')
        # plt.tight_layout()
    else:
        ax.set_xlabel(r'day \#', fontsize='xx-large')

    LegendHandler(fig, ax, _legend, posD, (textStr, ))
    fig.savefig(f'{ds.name}-color.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{ds.name}-color.png', dpi=300, bbox_inches='tight')


# # SST/T2 jumps
def display_jumps(ds, w):
    SHOWSIG = True

    cdf = ds.data

    if SHOWSIG:
        fig = plt.figure(figsize=(15, 11), clear=True, num=3, layout='constrained')
        (ax, ax1) = fig.subplots(2, 1, sharex=True, height_ratios=[.8, .2])
    else:
        fig, ax = plt.subplots(1, 1, figsize=(16, 16), num=3, clear=True)

    ax.plot(cdf.index, cdf['mean'], 'b', lw=4, label=f'{ds.mean_baseL[0]}-{ds.mean_baseL[1]} mean')
    ax.plot(cdf.index, cdf['mean+2std'], 'k', lw=4, label=rf'{ds.mean_baseL[0]}-{ds.mean_baseL[1]} mean $\pm 2\sigma$')
    ax.plot(cdf.index, cdf['mean-2std'], 'k', lw=4)
    ax.plot(cdf.index, cdf['mean+5std'], 'r', lw=4, label=rf'{ds.mean_baseL[0]}-{ds.mean_baseL[1]} mean +$5\sigma$')
    ax.plot(cdf.index, cdf['mean+5finalstd'], 'g', lw=4, label=f'{ds.mean_baseL[0]}-{ds.mean_baseL[1]}'r' mean +$5\sigma_{final}$')

    num_colors = len(ds.yearL)
    cmap = plt.get_cmap(w.value, num_colors)

    _skip = 10
    ms = 5

    for i, _year in enumerate(ds.yearL):
        if _year in ds.jump_years.keys():
            ax.plot(cdf.index, cdf.loc[:, f'{_year}'], c=cmap(i), lw=1)
            day_num = cdf.index.to_list()[0::_skip]
            data = cdf[f'{_year}'].iloc[0::_skip]
            ax.plot(day_num, data, ds.jump_years[_year], c=cmap(i), ms=ms, label=f'{_year}')

    posD = ds.pos_jump
    ax.set_ylabel(r'$T\ \  (^{\circ} C)$', fontsize='xx-large')
    ax.set_title(ds.jump_title, fontsize='xx-large')
    ax.text(275, posD['caption_height'],  f'\\it{{data current to {ds.last_date}}}', fontsize='xx-large')
    ax.set_xlim(0, 365)
    ax.set_ylim(*posD['ylim'])
    _legend = ax.legend(ncols=posD['ncols'], loc=posD['loc'])

    if SHOWSIG:
        ax1.plot(cdf.index, cdf['std'], 'r', ms=4)
        ax1.plot(cdf.index, cdf['finalstd'], 'g', ms=4)
        ax1.set_ylabel(r'$\sigma\ \  (^{\circ} C)$', fontsize='xx-large')
        ax1.set_xlabel(r'day \#', fontsize='xx-large')
        ax1.set_title(r'daily $\sigma$', fontsize='xx-large')
    else:
        ax.set_xlabel(r'day \#', fontsize='xx-large')

    textStr0 = rf'$\sigma = \sigma({ds.mean_baseL[0]}-{ds.mean_baseL[1]})$'
    textStr1 = r'$\sigma_{final} = \sigma(2011-2024)$'

    LegendHandler(fig, ax, _legend, posD, (textStr0, textStr1))

    fig.savefig(f'{ds.name}-color-jumps.pdf', bbox_inches='tight')
    fig.savefig(f'{ds.name}-color-jumps.png', bbox_inches='tight')


class LegendHandler:
    drawn = 0

    def __init__(self, fig, ax, legend, posD, textL):
        self.legend = legend
        self.ax = ax
        self.fig = fig
        self.posD = posD
        self.drawn = 0
        self.textL = textL
        fig.canvas.mpl_connect('draw_event', self.on_draw)
        fig.canvas.draw()

    def on_draw(self, event):
        self.drawn = self.drawn+1
        if self.drawn == 1:
            bbox = self.legend.get_window_extent().transformed(self.ax.transData.inverted())
            if 0:
                print(bbox, f'{self.drawn=}')
            voff = self.posD['sigma_caption_voffset']
            for i, textStr in enumerate(self.textL):
                self.ax.text(bbox.x0, bbox.y0-(i+1)*voff, textStr, fontsize='medium')


def display_as_one_curve(ds, w, rolling_only=False):
    if not rolling_only:
        fig_all_T, ax_all_T = plt.subplots(1, 1, figsize=(16, 14), clear=True, num=4)
    else:
        fig_all_T, ax_all_T = plt.subplots(1, 1, figsize=(16, 14), clear=True, num=5)
    num_colors = len(ds.yearL)
    cmap = plt.get_cmap(w.value, num_colors)
    cdf = ds.data

    if not rolling_only:
        for i, _year in enumerate(ds.yearL):
            if i == 0:
                label = f'daily {ds.name} data'
            else:
                label = None
            ax_all_T.plot(array(cdf.index.to_list())/365.+i+int(ds.first_year), cdf[f'{_year}'], '-', c=cmap(i), ms=1, lw=1, label=label)

    onedf = DataFrame(cdf[f'{ds.yearL[0]}'].to_list())
    for i, _year in enumerate(ds.yearL[1:]):
        df1 = DataFrame(cdf[f'{_year}'].to_list())
        onedf = concat([onedf, df1]).reset_index(drop=True)

    center_window = False
    for num_year_roll, _color in zip((1, 2, 3, 4), ('k', 'g', 'r', 'b')):
        window_length = num_year_roll*365
        min_periods = window_length-num_year_roll
        lab = f'rolling {num_year_roll}-year average'
        ax_all_T.plot(array(onedf.index.to_list())/365.+int(ds.first_year),
                      onedf[0].rolling(window=window_length,
                                       center=center_window,
                                       min_periods=min_periods).mean(),
                      color=_color, ms=1, label=lab)
    title = 'centered' if center_window else 'right_aligned'
    title = title + ' rolling window'
    plt.title(title)
    ax_all_T.set_ylabel(r'$T\ \  (^{\circ} C)$', fontsize='xx-large')
    plt.legend(numpoints=5)

    if not rolling_only:
        plt.colorbar(ds.sm, ax=ax_all_T, ticks=linspace(int(ds.first_year), ds.yearL[-1], num_colors),
                     boundaries=arange(int(ds.first_year)-.5, ds.yearL[-1]+1, 1))
    fN = f'{ds.name}-one-curve-rolling-only.png' if rolling_only else f'{ds.name}-one-curve.png'
    plt.savefig(fN)


# # Just checking my stuff against Climate Reanalyzer graph

def display_shifting_averages(ds):
    markerL = ('g>', 'k<', 'c^')
    cdf = ds.data

    fig, ax = plt.subplots(1, 1, figsize=(16, 8), num=6, clear=True)

    cdf['3'] = cdf.loc[:, '2001':'2020'].mean(axis=1, skipna=True)
    cdf['4'] = cdf.loc[:, '2011':'2025'].mean(axis=1, skipna=True)
    cdf['5'] = cdf.loc[:, '2021':'2025'].mean(axis=1, skipna=True)

    ms = 4

    sl = slice(0, len(cdf.index), 10)
    date = cdf.index[sl]

    num = len(ds.end_data_period)
    for end_data, marker in zip(ds.end_data_period, markerL[0:num]):
        ax.plot(date, cdf[end_data][sl], marker, ms=ms, label=f'{end_data} given daily mean')

    ax.plot(date, cdf['3'][sl], 'y<', ms=ms, label='2001-2020 daily mean')
    ax.plot(date, cdf['4'][sl], 'b<', ms=ms, label='2011-2025 daily mean')
    ax.plot(date, cdf['5'][sl], 'c<', ms=ms, label='2021-2025 daily mean')

    lw = 1
    sl = slice(0, len(cdf.index), 1)
    date = cdf.index[sl]
    ax.plot(date, cdf['mean-2std'][sl], 'k-.', lw=lw, ms=ms, label=rf'{ds.mean_baseL[0]}-{ds.mean_baseL[1]} mean -$2\sigma$')
    ax.plot(date, cdf['mean'][sl], 'b-.', lw=lw, ms=ms, label=rf'{ds.mean_baseL[0]}-{ds.mean_baseL[1]} mean')
    ax.plot(date, cdf['mean+2std'][sl], 'k-.', lw=lw, ms=ms, label=rf'{ds.mean_baseL[0]}-{ds.mean_baseL[1]} mean +$2\sigma$')
    ax.plot(date, cdf['mean+5std'][sl], 'r-.', lw=lw, ms=ms, label=rf'{ds.mean_baseL[0]}-{ds.mean_baseL[1]} mean +$5\sigma$')
    ax.plot(date, cdf['mean+5finalstd'][sl], 'g', lw=lw, ms=ms, label=f'{ds.mean_baseL[0]}-{ds.mean_baseL[1]} 'r'mean +$5\sigma_{final}$')
    ax.legend(numpoints=4, ncols=2)
    ax.set_xlim(0, 365)
    ax.set_title(ds.name)
    plt.savefig(f'{ds.name}-change-in-averages')


def make_all_plots(ds):
    get_data(ds)
    w = construct_colormap(show_chart=True)
    display_latest_year_data(ds)
    display_main_plot(ds, w)
    display_jumps(ds, w)
    display_as_one_curve(ds, w)
    display_as_one_curve(ds, w, rolling_only=True)
    display_shifting_averages(ds)


if __name__ == '__main__':
    _ds = T2
    if 1:
        update_data(_ds)
        make_all_plots(_ds)

    plt.show()
