import ast
import copy
import os
import os.path
import time
from typing import Dict, Tuple, Optional, Set

import bokeh.colors
import bokeh.colors.named
import bokeh.embed
import bokeh.events
import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.models.callbacks
import bokeh.models.tools
import bokeh.plotting
import bokeh.resources
import colorcet
import numpy as np
import pandas as pd
import polars as pl

import python_array_utils as utils

from graphing_utils import export

ukb = os.environ['UKB']

chr_lens = np.genfromtxt(
    f'{ukb}/misc_data/genome/chr_lens.txt',
    skip_header=1,
    usecols=(1),
    dtype=int
)

max_p_val = 300 # in -log10

def replace_last(text, sub, rep):
    return rep.join(text.rsplit(sub, 1))

def fix_cols(header):
    return replace_last(replace_last(header, '0.05_significance_CI', 'foo'), '5e-8_significance_CI', 'bar').split('\t')

def get_conditioned_strs(condition):
    splits = condition.split('_')
    return [int(STR) for STR in splits[(splits.index('STR')+1):(splits.index('ISNP')-1)]]

def get_conditioned_isnps(condition):
    splits = condition.split('_')
    return splits[(splits.index('ISNP')+1):(splits.index('ASNP')-1)]

def create_source_dict(
        data: np.recarray,
        start_chrom: Optional[int],
        cols_to_skip: Set[str] = set(),
        cols_to_include: Set[str] = set(),
        chrs_to_var_signals = None,
        start_p_val_cap = None
    ) -> Tuple[bokeh.models.ColumnDataSource, Dict[int, bokeh.models.ColumnDataSource]]:
    if not start_chrom:
        cds_range = range(1, 23)
        start_chrom = 1
    else:
        cds_range = [start_chrom]

    sources = {}
    assert len(cols_to_skip) == 0 or len(cols_to_include) == 0
    for field in 'chr', 'pos', 'p_val':
        if field not in data.dtype.names:
            print(field, flush=True)
            assert False
    for chrom in cds_range:
        chrom_dict = {}
        idx = data['chr'] == chrom
        for name in data.dtype.names:
            if len(cols_to_skip) > 0:
                if name in cols_to_skip:
                    continue
            elif len(cols_to_include) > 0:
                if name not in cols_to_include:
                    continue
            chrom_dict[name] = data[name][idx]
        sources[chrom] = bokeh.models.ColumnDataSource(chrom_dict)
        if start_p_val_cap is not None:
            sources[chrom].data['display_p_val'] = np.minimum(
                sources[chrom].data['p_val'], start_p_val_cap
            )
        else:
            sources[chrom].data['display_p_val'] = sources[chrom].data['p_val'].copy()

        if chrs_to_var_signals:
            colnames_to_merge_on = list(chrs_to_var_signals[chrom].columns)[:-1]
            merge_cols = pd.DataFrame(data[colnames_to_merge_on][idx])

            merged_data = merge_cols.merge(
                chrs_to_var_signals[chrom].astype({'pos': int}),
                on=colnames_to_merge_on,
                how='left'
            )
            assert merged_data.shape[1] == len(chrs_to_var_signals[chrom].columns)
            assert merged_data.shape[0] == merge_cols.shape[0]
            sources[chrom].data['FINEMAP_pcausal'] = merged_data['pcausal']

    copy_source = bokeh.models.ColumnDataSource(copy.deepcopy(sources[start_chrom].data))

    return copy_source, sources

my_results_rename = {
    0: 'chr',
    3: 'filtered',
    4: 'p_val',
    5: 'coeff_phenotype'
}
my_str_results_rename = {
    '0.05_significance_CI': 'CI5e_2SingleDosagePhenotype',
    '5e-8_significance_CI': 'CI5e_8SingleDosagePhenotype',
}

def load_my_str_results(phenotype, binary, unconditional_results_fname, conditional_results_fname = None):
    print(f"Loading my STR results for {phenotype} ... ", end='', flush=True)
    start_time = time.time()
    with open(unconditional_results_fname) as tsv:
        header = tsv.readline().strip()
    unconditional_results = pl.scan_csv(
        unconditional_results_fname,
        sep='\t',
        skip_rows=1,
        has_header=False,
        with_column_names = lambda _: fix_cols(header),
        dtypes={'alleles': str, 'locus_filtered': str}
    ).filter(pl.col(f'p_{phenotype}') <  5e-5).collect().to_pandas()

    if not conditional_results_fname:
        results = unconditional_results
    else:
        results = pd.read_csv(
            conditional_results_fname,
            header=0,
            delimiter='\t',
            encoding='UTF-8',
            dtype=utils.get_dtypes(conditional_results_fname, {'locus_filtered': str})
        )

        unconditional_results[f'p_{phenotype}'] = np.maximum(unconditional_results[f'p_{phenotype}'], 1 / 10**max_p_val)
        unconditional_results[f'p_{phenotype}'] = -np.log10(unconditional_results[f'p_{phenotype}'])
        unconditional_results.rename(
            columns = {f'p_{phenotype}': 'unconditional_p'}, inplace=True
        )
        unconditional_results = unconditional_results[['chrom', 'pos', 'unconditional_p']]

        results = results.merge(
            unconditional_results,
            on=['chrom', 'pos'],
            how= 'inner'
        ) # subsets to only those which passed the p-val threshold in the unconditional run

    if binary == 'logistic':
        results.rename(columns={'firth?': 'firth'}, inplace=True)

    rename_dict = {}
    for idx, name in my_results_rename.items():
        rename_dict[results.columns[idx]] = name
    rename_dict.update(my_str_results_rename)
    for colname in ('total_per_allele_dosages', 'total_hardcall_alleles',
                'subset_total_per_allele_dosages', 'subset_total_hardcall_alleles',
                'subset_allele_dosage_r2'):
        # convert allele lens from strings to floats, in addition round allele lens and values, but not NaN values
        new_col = np.array(list(map(
            lambda dict_str: {round(float(allele_len), 2): (round(val, 2) if val != 'NaN' else val) for allele_len, val in ast.literal_eval(dict_str).items()},
            results[colname]
        )))
        # convert allele_lens to ints if they are close enough
        new_col = np.array(list(map(
            lambda d: str({(int(key) if key == int(key) else key) : val for key, val in d.items()}),
            new_col
        )))
        results[colname] = new_col
    results.rename(columns=rename_dict, inplace=True)
    results = utils.df_to_recarray(results)
    results['p_val'] = np.maximum(results['p_val'], 1 / 10**max_p_val)
    results['p_val'] = -np.log10(results['p_val'])
    if conditional_results_fname:
        for STR in get_conditioned_strs(conditional_results_fname):
            results = results[results['pos'] != STR]
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return results

def load_my_snp_results(phenotype, binary, chrom, start, end):
    print(f"Loading my {phenotype} snp results chrom={chrom} start={start} end={end} ... ", end='', flush=True)
    start_time = time.time()

    snp_results_fname = f'{ukb}/association/plots/input/{phenotype}/my_imputed_snp'
    if binary:
        snp_results_fname += '_' + binary
    snp_results_fname += f'_chr{chrom}'
    if start:
        snp_results_fname += f'_{start}_{end}'
    snp_results_fname += '_results.tab'

    type_overrides = {
        'locus_filtered': str,
        f'p_{phenotype}': float,
        '0.05_significance_CI': str,
        '5e-8_significance_CI': str
    }
    if not binary:
        type_overrides.update({f'mean_{phenotype}_per_single_dosage': str})
    else:
        type_overrides.update({f'fraction_{phenotype}_per_single_dosage': str})

    if binary == 'logistic':
        type_overrides.update({'unused_col': str, 'firth?': str})

    my_snp_results = utils.df_to_recarray(pd.read_csv(
        snp_results_fname,
        header=0,
        delimiter='\t',
        encoding='UTF-8',
        dtype=utils.get_dtypes(snp_results_fname, type_overrides)
    ))

    names = list(my_snp_results.dtype.names)
    for idx, name in my_results_rename.items():
        names[idx] = name
    my_snp_results.dtype.names = names

    my_snp_results['p_val'] = -np.log10(my_snp_results['p_val'])

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return my_snp_results

def load_plink_results(phenotype, binary, unconditional_results_fname, conditional_results_fname=None):
    # TODO remove conditional snps
    # Load plink SNP results
    print(f"Loading plink SNP results for {phenotype} ... ", end='', flush=True)
    start_time = time.time()

    if binary:
        binary_colnames = {
            'A1_CASE_CT': 'alt_case_count',
            'A1_CTRL_CT': 'alt_control_count',
            'FIRTH?': 'firth?'
        }
    else:
        binary_colnames = {}
    start_time = time.time()
    unconditional_results = pl.scan_csv(
        unconditional_results_fname,
        sep='\t',
        null_values='NA'
    ).filter(
        pl.col('P') <  5e-5
    ).rename({
        '#CHROM': 'chr',
        'POS': 'pos',
        'ID': 'id',
        'REF': 'ref',
        'ALT': 'alt',
        'P': 'p_val',
        'ERRCODE': 'error',
        # these last three only occur in logistic regression
        **binary_colnames
    }).select([
        pl.col(col) for col in ['chr', 'pos', 'id', 'ref', 'alt', 'p_val', 'error', *binary_colnames.values()]
    ]).collect().to_pandas()

    if not conditional_results_fname:
        results = unconditional_results
    else:
        results = pl.scan_csv(
            conditional_results_fname,
            sep='\t',
            null_values='NA'
        ).rename({
            '#CHROM': 'chr',
            'POS': 'pos',
            'ID': 'id',
            'REF': 'ref',
            'ALT': 'alt',
            'P': 'p_val',
            'ERRCODE': 'error',
            # these last three only occur in logistic regression
            **binary_colnames
        }).select([
            pl.col(col) for col in ['chr', 'pos', 'id', 'ref', 'alt', 'p_val', 'error', *binary_colnames.values()]
        ]).collect().to_pandas()

        unconditional_results['p_val'] = np.maximum(unconditional_results['p_val'], 1 / 10**max_p_val)
        unconditional_results['p_val'] = -np.log10(unconditional_results['p_val'])
        unconditional_results.rename(
            columns={'p_val': 'unconditional_p'},
            inplace=True
        )
        unconditional_results = unconditional_results[['chr', 'pos', 'unconditional_p']]

        results = results.merge(
            unconditional_results,
            on=['chr', 'pos'],
            how = 'inner'
        ) # subsets to only those which passed the p-val threshold in the unconditional run

    if binary == 'logistic':
        results.rename(columns={'firth?': 'firth'}, inplace=True)

    results = utils.df_to_recarray(results)

    results = results[results['error'] != 'CONST_OMITTED_ALLELE']
    if binary == 'logistic':
        # in theory could keep unfinished error codes and just note them,
        # but easier to ignore
        results = results[
            (results['error'] != 'FIRTH_CONVERGE_FAIL') &
            (results['error'] != 'UNFINISHED')
        ]
    results['p_val'] = np.maximum(results['p_val'], 1 / 10**max_p_val)
    results['p_val'] = -np.log10(results['p_val'])

    # we've already filtered all the spots that had errors in the unconditional run
    # having a VIF_TOO_HIGH or CORR_TOO_HIGH only in the conditional run just means that
    # SNP is extremely correlated with the conditioning variants, which means
    # its p-value should be very small, so this isn't an issue.
    if not conditional_results_fname:
        if not np.all(results['error'] == '.'):
            print(np.unique(results['error']))
            assert False
    else:
        assert np.all(
            (results['error'] == '.') |
            (results['error'] == 'VIF_TOO_HIGH') |
            (results['error'] == 'CORR_TOO_HIGH')
        )
    # rename for readability
    results['error'][results['error'] == '.'] = 'none'
    results['p_val'][results['error'] == 'VIF_TOO_HIGH'] = 0

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return results

def load_gwas_catalog(phenotype):
    if phenotype not in {'height', 'total_bilirubin'}:
        return None, None
    # Load the NHGRI-EBI GWAS catalog
    print("Loading GWAS catalog results ... ", end='', flush=True)
    start_time = time.time()
    catalog = np.loadtxt(
        f'{ukb}/misc_data/snp_summary_stats/catalog/catalog_hg19.tsv',
        usecols=[11, 12, 27],
        delimiter='\t',
        skiprows=1,
        dtype=object
    ) # chrom, pos, pvalue
    # omit results that are not mapped to a recognizable chromosome
    catalog_names = np.loadtxt(
        f'{ukb}/misc_data/snp_summary_stats/catalog/catalog_hg19.tsv',
        usecols=[7, 21, 34],
        delimiter='\t',
        skiprows=1,
        dtype=object
    ).astype('U') # cols 7, 34 describe the 'height', 21 is the snp rsid
    catalog_names = np.char.lower(catalog_names)
    filter_weird_chroms = np.isin(catalog[:, 0],
                                  list(str(chrom) for chrom in range(1, 23)))
    catalog_names = catalog_names[filter_weird_chroms, :]
    catalog = catalog[filter_weird_chroms, :]
    catalog = catalog.astype(float)

    known_assocs = {}
    known_assoc_ids = {}
    height_rows = np.logical_and(
        catalog_names[:, 2] == 'body height',
        catalog_names[:, 0] != "pericardial adipose tissue adjusted for height and weight"
    )
    known_assocs['height'] = catalog[height_rows, :]
    known_assoc_ids['height'] = catalog_names[height_rows, 1]
    known_assocs['total_bilirubin'] = \
            catalog[catalog_names[:, 2] == 'bilirubin measurement', :]
    known_assoc_ids['total_bilirubin'] = \
            catalog_names[catalog_names[:, 2] == 'bilirubin measurement', 1]

    for assoc_phen in known_assocs:
        zero_idx = known_assocs[assoc_phen][:, 2] == 0
        known_assocs[assoc_phen][~zero_idx, 2] = -np.log10(known_assocs[assoc_phen][~zero_idx, 2])
        known_assocs[assoc_phen][zero_idx, 2] = 50.12345 #choose an arbitrary number

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)

    return known_assocs[phenotype], known_assoc_ids[phenotype]

def load_finemap_signals(finemap_signals):
    # finemap signals must be sorted
    print("Loading FINEMAP causality results ... ", end='', flush=True)
    start_time = time.time()
    chrs_to_snp_signals = {}
    chrs_to_str_signals = {}
    chrs_to_regions = {}
    for i in range(1, 23):
        chrs_to_snp_signals[i] = []
        chrs_to_str_signals[i] = []
        chrs_to_regions[i] = []
    for signal in finemap_signals:
        region = signal.split('/')[-2]
        chrom, start, end = (int(val) for val in region.split('_'))
        chrs_to_regions[chrom].append(f'{start}_{end}')
        with open(signal) as per_var_output:
            next(per_var_output)
            for line in per_var_output:
                _id, pos, pcausal = np.array(line.split())[[1, 3, 10]]
                pos = int(pos)
                split = _id.split('_')
                assert int(split[1]) == pos
                if pcausal == 'NA':
                    pcausal = np.nan
                else:
                    pcausal = float(pcausal)
                if _id[:4] == 'STR_':
                    chrs_to_str_signals[chrom].append((pos, pcausal))
                elif _id[:4] == 'SNP_':
                    chrs_to_snp_signals[chrom].append(
                        (pos, split[2], split[3], pcausal)
                    )
                else:
                    raise ValueError(f'Found uninterpretable id {_id}')

    chrs_to_snp_signals = {
        key: (
            pd.DataFrame(np.stack(val), columns=('pos', 'ref', 'alt', 'pcausal'))
            if len(val) > 0
            else pd.DataFrame(      [], columns=('pos', 'ref', 'alt', 'pcausal'))
        )
        for key, val in chrs_to_snp_signals.items()
    }
    chrs_to_str_signals = {
        key: (
            pd.DataFrame(np.stack(val), columns=('pos', 'pcausal'))
            if len(val) > 0
            else pd.DataFrame(      [], columns=('pos', 'pcausal'))
        )
        for key, val in chrs_to_str_signals.items()
    }

    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return (chrs_to_snp_signals, chrs_to_str_signals, chrs_to_regions)

def get_my_str_run_date(run_readme_file):
    with open(run_readme_file) as README:
        date_line = next(README)
        return date_line.split(' ')[2]

def get_plink_snp_run_date(run_log_file):
    with open(run_log_file) as log_file:
        for line in log_file:
            line = line.strip()
            if line.startswith('Start time: '):
                return line[12:]
    assert False

def display_my_str_run_date(manhattan_plot, my_str_run_date, fontsize=18):
    manhattan_plot.add_layout(bokeh.models.Title(
        text=f"My STR code run date: {my_str_run_date}",
        align='right',
        text_font_size=f'{fontsize}px'
    ), 'below')

def display_plink_snp_run_date(manhattan_plot, plink_snp_run_date, fontsize=18):
    manhattan_plot.add_layout(bokeh.models.Title(
        text=f"Plink SNP run date: {plink_snp_run_date}",
        align='right',
        text_font_size=f'{fontsize}px'
    ), 'below')

def load_and_merge_peaks_into_dfs(peaks_fname, my_str_results,  plink_snp_results):
    print("Adding peak data ... ", end='', flush=True)
    start_time = time.time()

    peaks = pd.read_csv(
        peaks_fname,
        delimiter='\t',
        header=0,
        dtype=utils.get_dtypes(peaks_fname, {'ref_(snp_only)': object, 'alt_(snp_only)': object})
    )
    peaks.rename(columns={'chrom': 'chr', 'ref_(snp_only)': 'ref', 'alt_(snp_only)': 'alt'}, inplace=True)
    peaks['marker'] = 1

    str_peaks = peaks[peaks['variant_type'] == 'STR']
    my_str_results = pd.DataFrame.from_records(my_str_results)
    my_str_results = my_str_results.merge(
        str_peaks[['chr', 'pos', 'marker']],
        how='left',
        on=['chr', 'pos']
    )
    my_str_results['is_peak'] = ~np.isnan(my_str_results['marker'])
    my_str_results = utils.df_to_recarray(my_str_results)

    snp_peaks = peaks[peaks['variant_type'] == 'SNP']
    plink_snp_results = pd.DataFrame.from_records(plink_snp_results)
    plink_snp_results = plink_snp_results.merge(
        snp_peaks[['chr', 'pos', 'ref', 'alt', 'marker']],
        how='left',
        on=['chr', 'pos', 'ref', 'alt']
    )
    plink_snp_results['is_peak'] = ~np.isnan(plink_snp_results['marker'])
    plink_snp_results = utils.df_to_recarray(plink_snp_results)
    print(f"done ({time.time() - start_time:.2e}s)", flush=True)
    return my_str_results, plink_snp_results

def full_genome_x_axis(plot):
    pre_chr_sums = np.cumsum([0, *chr_lens[:-1]])
    mid_points = [int(num) for num in pre_chr_sums + (chr_lens//2)]
    plot.xaxis.ticker = mid_points
    plot.xaxis.major_label_overrides = {
        mid_points[chrom - 1]: str(chrom) for chrom in range(1, 23)
    }

def full_genome_pandas_df(df):
    df['plot_pos'] = df['pos']
    for chrom in range(2, 23):
        df['plot_pos'][df['chr'] >= chrom] += chr_lens[chrom - 2]

def full_genome_polars_df(df):
    df = df.with_column(
        pl.col('pos').alias('plot_pos')
    )
    for chrom in range(2, 23):
        df = df.with_column(pl
            .when(pl.col('chr') >= chrom)
            .then(pl.col('plot_pos') + int(chr_lens[chrom - 2]))
            .otherwise(pl.col('plot_pos'))
            .alias('plot_pos')
        )
    return df

def add_x_navigate(bokeh_plot):
    wheel_zoom = bokeh.models.tools.WheelZoomTool(dimensions="width")
    bokeh_plot.add_tools(wheel_zoom)
    bokeh_plot.toolbar.active_scroll = wheel_zoom

    pan = bokeh.models.tools.PanTool(dimensions="width")
    bokeh_plot.add_tools(pan)
    bokeh_plot.toolbar.active_drag = pan

    bokeh_plot.add_tools(bokeh.models.tools.ZoomInTool(dimensions='width'))
    bokeh_plot.add_tools(bokeh.models.tools.ZoomOutTool(dimensions='width'))
