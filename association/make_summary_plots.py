import argparse
import os

from matplotlib.colors import CSS4_COLORS as colors
import matplotlib.pyplot as plt
import matplotlib.image
import numpy as np
import numpy.ma

import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
import scipy.ndimage
import matplotlib.patches

import load_and_filter_genotypes

ukb = os.environ['UKB']
GENOME_WIDE_TRANS_SIG = -np.log10(5e-8)
HEIGHT_CUTOFF = 20

# assume data.shape is (loci, 3)
# where columns are chrom, pos, pval
# will return sorted from least to most significant
# with values in the transformed axis (-np.log10(pval))
def prep_data(data):
    data = data[~np.isnan(data[:, 2]), :]
    pvals = data[:, 2]
    pvals[pvals == 0] = min(pvals[pvals != 0])
    pvals[:] = -np.log10(pvals)
    return data[np.argsort(data[:, 2])]


# gets an array of x-axis vals of the corresponding length
def qq_x_axis(length):
    uniform = (np.arange(length) + 0.5)/length
    trans_uniform = -np.log10(uniform)[::-1]
    return trans_uniform


def plot_qq(ax, prepped_pvals, label, color, height_cutoff):
    if height_cutoff:
        prepped_pvals = prepped_pvals.copy()
        prepped_pvals[prepped_pvals > HEIGHT_CUTOFF] = HEIGHT_CUTOFF
    xs = qq_x_axis(len(prepped_pvals))
    ax.scatter(
        xs,
        prepped_pvals,
        color=colors[color],
        marker='.',
        label=label
    )


def make_qq_plots(phenotypes, results, snp_summary_stats, run_name,
                  height_cutoff):
    # plot QQ plots
    nphenotypes = len(phenotypes)
    fig, axs = plt.subplots(nphenotypes, figsize=(10, 5*nphenotypes))
    fig.suptitle('QQ plots')
    fig.tight_layout(pad=5.0)
    plt_count = 0
    for phenotype in phenotypes:
        nresults = results[phenotype].shape[0]
        if len(phenotypes) > 1:
            ax = axs[plt_count]
        else:
            ax = axs
        ax.set_title(f'{phenotype}')
        ax.set_ylabel('-log_10(pval)')
        ax.set_xlabel(f'-log_10(Uniform dist) (# STR loci = {nresults})')

        if phenotype in snp_summary_stats:
            max_len = max(snp_summary_stats[phenotype].shape[0], nresults)
            plot_qq(ax, snp_summary_stats[phenotype][:, 2],
                    'SNPs', 'royalblue',
                    height_cutoff)
        else:
            max_len = nresults

        plot_qq(ax, results[phenotype][:, 2], 'STRs',
                'mediumvioletred', height_cutoff)
        plot_qq(ax, qq_x_axis(max_len),
                'y=x, expectation under the null', 'gray',
                height_cutoff)

        crosses_gws = np.argmax(results[phenotype][:, 2]
                                > GENOME_WIDE_TRANS_SIG)
        xy = (qq_x_axis(nresults)[crosses_gws],
              results[phenotype][crosses_gws, 2])
        n_sig_loci = np.sum(results[phenotype][:, 2]
                            > GENOME_WIDE_TRANS_SIG)
        ax.annotate(
            f"{n_sig_loci} loci with p <= 5e-8",
            xy=xy,
            xytext=(xy[0]-1.5, xy[1]+5),
            arrowprops=dict(facecolor='black', shrink=0.05),
        )

        ax.legend()

        plt_count += 1

    cutoff_txt = ''
    if not height_cutoff:
        cutoff_txt = '_no_cutoff'
    plt.savefig(f'{ukb}/association/runs/{run_name}/qq_plot{cutoff_txt}.png')


def plot_manhattan(ax, data, label, color_pair, chr_lens, height_cutoff):
    if height_cutoff:
        data = data.copy()
        data[data[:, 2] > HEIGHT_CUTOFF, 2] = HEIGHT_CUTOFF
    first = True
    for chrom in range(1, 23):
        if not first:
            label = None
        chrom_data = data[data[:, 0] == chrom, :]
        xs = chrom_data[:, 1]
        xs += np.sum(chr_lens[0:(chrom-1)])
        ax.scatter(
            xs,
            chrom_data[:, 2],
            color=colors[color_pair[chrom % 2]],
            label=label,
            marker='.'#,
            #alpha=0.5
        )
        first = False


def make_genome_manhattan_plots(phenotypes, results, snp_summary_stats,
                                snp_ss_description, known_assocs, run_name,
                                height_cutoff):
    # plot Manhattan plots
    nphenotypes = len(phenotypes)
    fig, axs = plt.subplots(nphenotypes, figsize=(300, 5*nphenotypes))
    fig.suptitle('Manhattan plots')
    fig.tight_layout(pad=5.0)
    chr_lens = np.genfromtxt(
        f'{ukb}/misc_data/genome/chr_lens.txt',
        skip_header=1,
        usecols=(1),
        dtype=int
    )
    max_len = np.sum(chr_lens)

    plt_count = 0
    for phenotype in phenotypes:
        if len(phenotypes) > 1:
            ax = axs[plt_count]
        else:
            ax = axs
        ax.set_title(f'{phenotype}')
        ax.set_ylabel('-log_10(pval)')
        ax.set_xlabel('Position (1e7 bp)')

        ax.plot([0, max_len], [GENOME_WIDE_TRANS_SIG] * 2, label='p = 5e-8',
                 color=colors['red'])

        if phenotype in snp_summary_stats:
            desc = f'SNPs ({snp_ss_description[phenotype]})'
            plot_manhattan(ax, snp_summary_stats[phenotype], desc,
                           ('royalblue', 'cornflowerblue'), chr_lens,
                           height_cutoff)

        plot_manhattan(ax, results[phenotype],
                       'STRs (Ancestry: European, n=500,000)',
                       ('mediumvioletred', 'deeppink'), chr_lens,
                       height_cutoff)

        if phenotype in known_assocs:
            plot_manhattan(ax, known_assocs[phenotype], 'GWAS Catalog Hits',
                           ('goldenrod', 'gold'), chr_lens, height_cutoff)


        ax.legend()

        ticks = []
        tick_labels = []
        for chrom, length in enumerate(chr_lens):
            chrom += 1
            new_ticks = np.arange(0, length, int(5e7))
            new_ticks += np.sum(chr_lens[0:(chrom-1)])
            ticks.extend(new_ticks)
            tick_labels.append(f'chr{chrom}:0')
            tick_labels.extend(str(int(tick)) for tick in np.arange(5, length/1e7, 5))
        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels)

        plt_count += 1

    cutoff_txt = ''
    if not height_cutoff:
        cutoff_txt = '_no_cutoff'
    plt.savefig(f'{ukb}/association/runs/{run_name}/manhattan_plot_all{cutoff_txt}.png')


def make_region_plot(results, snp_summary_stats,
                     snp_ss_description, known_assocs, run_name,
                     imputation_run_name, region):
    chrom, locus, phenotype = region.split(":")
    chrom = int(chrom)
    if '-' in locus:
        start, end = (int(x) for x in locus.split('-'))
    else:
        pos = int(locus)
        start, end = pos - int(5e5), pos + int(5e5)

    os.makedirs(f'{ukb}/association/runs/{run_name}/regional_plots', exist_ok=True)
    fig_height = 10
    fig_width = 8
    width_sides = 0.1
    height_sides = width_sides * fig_width/fig_height
    absolute_height_spacing = 0.15 * 8
    height_spacing = absolute_height_spacing/fig_height
    absolute_width_spacing = 0.8
    width_spacing = absolute_width_spacing/fig_width
    cbar_width_abs = 0.5
    cbar_width = cbar_width_abs/fig_width
    fig = plt.figure(figsize = (fig_width, fig_height))
    manhattan_ax = fig.add_axes([
        width_sides,
        fig_width/fig_height*(1 - 2*width_sides - width_spacing - cbar_width) + height_spacing,
        1 - 2*width_sides - width_spacing - cbar_width,
        1 - (fig_width/fig_height*(1 - 2*width_sides - width_spacing - cbar_width) + height_spacing) - fig_width/fig_height*width_sides,
    ])
    #from 
    # https://matplotlib.org/gallery/axisartist/demo_floating_axes.html
    '''
    tr = Affine2D().scale(1, 1).rotate_deg(45)
    grid_helper = floating_axes.GridHelperCurveLinear(
        tr, extremes=(-0.5, 3.5, 0, 4),
        grid_locator1=MaxNLocator(nbins=4),
        grid_locator2=MaxNLocator(nbins=4))
    '''
    heatmap_ax = fig.add_axes([
        width_sides,
        height_sides,
        1 - 2*width_sides - width_spacing - cbar_width,
        fig_width/fig_height*(1 - 2*width_sides - width_spacing - cbar_width)
    ])
    cbar_ax = fig.add_axes([
        1 - width_sides - cbar_width,
        height_sides,
        cbar_width,
        fig_width/fig_height*(1 - 2*width_sides - width_spacing - cbar_width)
    ])
    
    fig.suptitle(f'Associations with {phenotype}', fontsize='large')
    # fig.tight_layout(pad=5.0)
    #'''
    make_region_manhattan_plot(
        manhattan_ax,
        results,
        snp_summary_stats,
        snp_ss_description,
        known_assocs,
        chrom,
        start,
        end,
        phenotype
    )
    #'''
    make_ld_heatmap(
        heatmap_ax,
        cbar_ax,
        imputation_run_name,
        chrom,
        start,
        end
    )

    plt.savefig(f'{ukb}/association/runs/{run_name}/regional_plots/{phenotype}_{chrom}_{start}_{end}.png')

def make_region_manhattan_plot(ax, results, snp_summary_stats,
                         snp_ss_description, known_assocs,
                         chrom, start, end, phenotype):
    chr_lens = np.genfromtxt(
        f'{ukb}/misc_data/genome/chr_lens.txt',
        skip_header=1,
        usecols=(1),
        dtype=int
    )

    def subset(arr, chrom, start, end):
        arr = arr[arr[:, 0] == chrom, :]
        arr = arr[arr[:, 1] >= start, :]
        return arr[arr[:, 1] <= end, :]

    ax.set_ylabel('-log_10(pval)')
    ax.set_xlabel('Position')

    offset = np.sum(chr_lens[0:(chrom-1)])
    ax.plot([start+offset, end+offset], [GENOME_WIDE_TRANS_SIG] * 2, label='p = 5e-8',
             color=colors['red'])

    if phenotype in snp_summary_stats:
        desc = f'SNPs ({snp_ss_description[phenotype]})'
        plot_manhattan(ax,
                       subset(snp_summary_stats[phenotype], chrom, start, end),
                       desc,
                       ('royalblue', 'cornflowerblue'),
                       chr_lens, False)

    plot_manhattan(ax,
                   subset(results[phenotype], chrom, start, end),
                   'STRs (Ancestry: European, n=500,000)',
                   ('mediumvioletred', 'deeppink'),
                   chr_lens, False)

    if phenotype in known_assocs:
        plot_manhattan(ax,
                       subset(known_assocs[phenotype], chrom, start, end),
                       'GWAS Catalog Hits',
                       ('goldenrod', 'gold'),
                       chr_lens, False)
    ax.legend()


def make_ld_heatmap(ax, cbar_ax, imputation_run_name, chrom, start, end):
    '''
    gts_per_locus_list = []
    poses = []
    count = 0
    for len_gts, _, pos, _, _ in \
            load_and_filter_genotypes.filtered_strs(
                imputation_run_name, f'{chrom}:{start}-{end}'):
        count += 1
        print(f"loaded str {count}", end="\r")
        gts_per_locus_list.append(len_gts)
        poses.append(pos)
        if count == 20:
            end = pos
            break
    gts_per_locus = np.stack(gts_per_locus_list, axis=-1)
    # shape is now participant x haplotype x locus
    # all haplotypes are equally useful here, regardless of which person
    # they come from, so combine the first two dimensions into one
    gts_per_locus = gts_per_locus.reshape(gts_per_locus.shape[0]*2,
                                          gts_per_locus.shape[2])
    gts_per_locus = numpy.ma.masked_invalid(gts_per_locus)
    # shape is hapltoype x locus
    corrcoefs = np.square(numpy.ma.corrcoef(gts_per_locus, rowvar=False))
    '''
    corrcoefs = np.array([
        [1.0, 0.058071465626636495, 0.31039026622741855,         0.11832356183390194, 0.2666239166247261, 0.3686336628398046,         0.036703468456753756, 0.12210677207762605, 0.07739021666281841,         0.3366440537450745, 0.04087319746675118, 0.2857763238357282,         0.1003312933352248, 0.050749336024172864, 0.11504098423537285,         0.10278600878593348, 0.0713853075335143, 0.0005882198401218372,         0.6596429911239123, 0.05612157977115239],
        [0.058071465626636495, 1.0, 0.17439041646473993,         0.144967015184881, 8.3677484131342e-05, 0.21559221216050953,         0.014825911232766132, 0.01857279872205382, 0.015255547303802117,         0.15968081145945595, 0.40815861087486754, 0.08914653789493271,         0.01911698601627397, 0.41693754653588977, 0.2672161253816999,         0.2067087037671558, 0.09272423301145619, 0.2659691819828519,         0.08194731624019555, 0.09222605754295815],
        [0.31039026622741855, 0.17439041646473993, 1.0,         0.38551567117737606, 0.7343243054661031, 0.8217461894541107,         0.15188934959254602, 0.4159141241942624, 0.2114553830145246,         0.8472791026094442, 0.04368056683910006, 0.0389213827552414,         0.369038804592888, 0.04449362332323188, 0.5081728406541677,         0.25876668357374333, 0.24239313465856385, 0.06495314281347224,         0.17532878446908992, 0.18753979540658586],
        [0.11832356183390194, 0.144967015184881, 0.38551567117737606,
         1.0, 0.7462206885167709, 0.24998359744811396,         0.22346954411248376, 0.6490160939335295, 0.14840730286125933,
         0.33516029122261903, 0.14777302415261417, 0.005907316802405638,
         0.564288706836558, 0.1417074234459837, 0.05189700886071217,         0.003655399326966783, 0.6225080542395888, 0.5752432475578055,
         0.02736581235072704, 0.5360220629255285],
        [0.2666239166247261, 8.3677484131342e-05, 0.7343243054661031,         0.7462206885167709, 1.0, 0.6700760945525388, 0.2642028779591854,
         0.7088775573885092, 0.2958783093536479, 0.7604641377025146,         6.983542768531956e-06, 0.00664560136458498, 0.5172200524665731,
         4.952655836873531e-05, 0.3247156995587446, 0.15739013337053775,         0.5382658075071728, 0.24996700110532735, 0.11861618397955083,
         0.4522720115026945],
        [0.3686336628398046, 0.21559221216050953, 0.8217461894541107,
         0.24998359744811396, 0.6700760945525388, 1.0,         0.10003124540013952, 0.3106301704836652, 0.21613028289487288,         0.9730567168361598, 0.18746774301413202, 0.0844664940169808,         0.28673785043178296, 0.19081473143416405, 0.7327025492269393,
         0.44308797119845167, 0.1487738115586195, 0.006543633620788875,         0.22763309745808125, 0.1112278947033518],
        [0.036703468456753756, 0.014825911232766132, 0.15188934959254602,    0.22346954411248376, 0.2642028779591854, 0.10003124540013952,         1.0, 0.26114186878884665, 0.09010575205392872,         0.11679798737566167, 0.013507489926201259, 0.11362917065992822,         0.12866754663076685, 0.014979075630485663, 0.04119134587106534,         0.00861892635131197, 0.1860675845121037, 0.05847253948401881,         0.08746186464625146, 0.17144604385012333],
        [0.12210677207762605, 0.01857279872205382, 0.4159141241942624,         0.6490160939335295, 0.7088775573885092, 0.3106301704836652,         0.26114186878884665, 1.0, 0.21403457661139338,         0.390111625473481, 0.16301451338044487, 0.0005304338043141785,         0.6629820336910666, 0.16583874933899864, 0.06454587244282098,         0.00010672131389127171, 0.5209823122141619, 0.44558330645628375,         0.04384179374835876, 0.4590256692595173],
        [0.07739021666281841, 0.015255547303802117, 0.2114553830145246,         0.14840730286125933, 0.2958783093536479, 0.21613028289487288,         0.09010575205392872, 0.21403457661139338, 1.0,         0.24757751947416348, 0.01363575507646984, 0.021542674459145515,         0.13024984905167544, 0.013569926260719381, 0.13109331215663364,         0.06384989157702367, 0.2138741444855752, 0.03931601647156595,         0.0568629174124288, 0.18176105430855347],
        [0.3366440537450745, 0.15968081145945595, 0.8472791026094442,         0.33516029122261903, 0.7604641377025146, 0.9730567168361598,         0.11679798737566167, 0.390111625473481, 0.24757751947416348,
         1.0, 0.1404875742798725, 0.05796312871685066,         0.34246178788637455, 0.14249308655987353, 0.6928852050292633,
         0.40927283216312144, 0.20694396926606445, 0.026127290331885752,
         0.20082206001798228, 0.1617429018277851],
        [0.04087319746675118, 0.40815861087486754, 0.04368056683910006,
         0.14777302415261417, 6.983542768531956e-06, 0.18746774301413202,
         0.013507489926201259, 0.16301451338044487, 0.01363575507646984,         0.1404875742798725, 1.0, 0.07814618589925182,
         0.16684308204409265, 0.9621898331200696, 0.34997387678094666,         0.5348639079954022, 0.09324185388311972, 0.4826224328807197,
         0.06967122699030351, 0.08393838158065155],
        [0.2857763238357282, 0.08914653789493271, 0.0389213827552414,
         0.005907316802405638, 0.00664560136458498, 0.0844664940169808,         0.11362917065992822, 0.0005304338043141785,
         0.021542674459145515, 0.05796312871685066, 0.07814618589925182,         1.0, 0.12949350351575606, 0.09031179276294585,         0.19091450250352562, 0.20319448214863198, 0.10190639880147831,         1.8531189180076446e-05, 0.7972804653394614, 0.0762162963091056],
        [0.1003312933352248, 0.01911698601627397, 0.369038804592888,                                                                                                                0.564288706836558, 0.5172200524665731, 0.28673785043178296,         0.12866754663076685, 0.6629820336910666, 0.13024984905167544,         0.34246178788637455, 0.16684308204409265, 0.12949350351575606,         1.0, 0.15730196543167624, 0.16945520238280368,         0.02160853389929099, 0.6111280544647345, 0.5782002029451881,         0.13121234407684249, 0.4992745922904188],
        [0.050749336024172864, 0.41693754653588977, 0.04449362332323188,         0.1417074234459837, 4.952655836873531e-05, 0.19081473143416405,         0.014979075630485663, 0.16583874933899864, 0.013569926260719381,         0.14249308655987353, 0.9621898331200696, 0.09031179276294585,         0.15730196543167624, 1.0, 0.3428637796344565,         0.5364709375683436, 0.08864814163740822, 0.4790811763638503,         0.07770298431330631, 0.08573821691306112],
        [0.11504098423537285, 0.2672161253816999, 0.5081728406541677,         0.05189700886071217, 0.3247156995587446, 0.7327025492269393,         0.04119134587106534, 0.06454587244282098, 0.13109331215663364,         0.6928852050292633, 0.34997387678094666, 0.19091450250352562,         0.16945520238280368, 0.3428637796344565, 1.0,         0.6350240818195402, 0.09372032543384298, 0.00019656741120168734,         0.17613645231728325, 0.06314593881929084],
        [0.10278600878593348, 0.2067087037671558, 0.25876668357374333,         0.003655399326966783, 0.15739013337053775, 0.44308797119845167,         0.00861892635131197, 0.00010672131389127171,
         0.06384989157702367, 0.40927283216312144, 0.5348639079954022,         0.20319448214863198, 0.02160853389929099, 0.5364709375683436,
         0.6350240818195402, 1.0, 0.056766270469723296,
         0.042838772162018465, 0.18882227157951142, 0.04049640326843064],
        [0.0713853075335143, 0.09272423301145619, 0.24239313465856385,
         0.6225080542395888, 0.5382658075071728, 0.1487738115586195,
         0.1860675845121037, 0.5209823122141619, 0.2138741444855752,         0.20694396926606445, 0.09324185388311972, 0.10190639880147831,
         0.6111280544647345, 0.08864814163740822, 0.09372032543384298,         0.056766270469723296, 1.0, 0.6140391801532973,
         0.1123206000388427, 0.8821447986410944],
        [0.0005882198401218372, 0.2659691819828519, 0.06495314281347224,
         0.5752432475578055, 0.24996700110532735, 0.006543633620788875,         0.05847253948401881, 0.44558330645628375, 0.03931601647156595,
         0.026127290331885752, 0.4826224328807197,         1.8531189180076446e-05, 0.5782002029451881, 0.4790811763638503,         0.00019656741120168734, 0.042838772162018465,         0.6140391801532973, 1.0, 0.0011227638445394537,
         0.5444783725254136],
        [0.6596429911239123, 0.08194731624019555, 0.17532878446908992,                                                                                                              0.02736581235072704, 0.11861618397955083, 0.22763309745808125,         0.08746186464625146, 0.04384179374835876, 0.0568629174124288,         0.20082206001798228, 0.06967122699030351, 0.7972804653394614,         0.13121234407684249, 0.07770298431330631, 0.17613645231728325,         0.18882227157951142, 0.1123206000388427, 0.0011227638445394537,         1.0, 0.1015663120524458],
        [0.05612157977115239, 0.09222605754295815, 0.18753979540658586,         0.5360220629255285, 0.4522720115026945, 0.1112278947033518,         0.17144604385012333, 0.4590256692595173, 0.18176105430855347,         0.1617429018277851, 0.08393838158065155, 0.0762162963091056,         0.4992745922904188, 0.08573821691306112, 0.06314593881929084,         0.04049640326843064, 0.8821447986410944, 0.5444783725254136,         0.1015663120524458, 1.0]])
    poses = [233669089, 233674690, 233685026, 233686648, 233705236, 233706170, 233707844, 233708271, 233712201, 233720636, 233728289, 233730244, 233732861, 233740056, 233744126, 233747763, 233752933, 
233760237, 233767702, 233768147]
    end = poses[-1]
    poses = np.array(poses)
    mesh1D = np.linspace(start, end, 500, endpoint=True)
    meshx, meshy = np.meshgrid(mesh1D, mesh1D)
    meshx_exp = np.expand_dims(meshx, axis=2)
    meshy_exp = np.expand_dims(meshy, axis=2)
    closestx = np.argmin(np.abs(meshx_exp - poses.reshape(1, 1, -1)), axis=2)
    closesty = np.argmin(np.abs(meshy_exp - poses.reshape(1, 1, -1)), axis=2)
    z = corrcoefs[closestx, closesty]
    z = z.reshape(meshx.shape)
    z = scipy.ndimage.rotate(z, 45)
    # do I need to reflect this across the x-axis?
    im = ax.imshow(z, extent=(start,end,start,end))
    im.set_clip_path(matplotlib.patches.Rectangle(
        (start, start), end-start, (end-start)/2,
        transform=ax.transData
    ))
    #tr = Affine2D().rotate_deg_around((end+start)/2, (end+start)/2, 45)
    #pcm = ax.imshow(z, extent=(start,end,start,end), transform = tr)
    #ax.set_xlim(-1e9, 1e9)
    #ax.set_ylim(-1e9, 1e9)
    ax.set_xlim(start, end)
    ax.set_ylim(start, (start+end)/2)
    ax.set_aspect('equal', adjustable='box')
    cbar = ax.figure.colorbar(im, cax=cbar_ax)
    cbar.ax.set_ylabel("Correlation between loci", rotation=-90, va="bottom")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_name")
    parser.add_argument("imputation_run_name")
    parser.add_argument(
        "--regions",
        help=("comma separated list of regions to make summary plots for. "
              "Each region is of the form chrom:start-end:phenotype or "
              "chrom:pos:phenotype in which case a 1Mb region around the "
              "position is shown ")
    )
    args = parser.parse_args()
    if args.regions:
        regions = args.regions.split(',')
    else:
        regions = None

    phenotypes = ['height', 'total_bilirubin']

    results = {}
    #'''
    for phenotype in phenotypes:
        results[phenotype] = np.loadtxt(
            f'{ukb}/association/runs/{args.run_name}/results/{phenotype}.txt',
            usecols=[0, 1, 4],
            skiprows=1
        )
    #TODO fix nan alleles!
    for phenotype in phenotypes:
        results[phenotype] = prep_data(results[phenotype])
    #'''
    snp_summary_stats = {}
    snp_ss_description = {}
    #'''
    snp_summary_stats['height'] = np.loadtxt(
        (f"{ukb}/misc_data/snp_summary_stats/height/"
         "Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt"),
        usecols=(0, 1, 8),
        skiprows=1
    )
    snp_ss_description['height'] = 'Ancestry: European, n=700,000'
    snp_summary_stats['total_bilirubin'] = np.loadtxt(
        (f"{ukb}/misc_data/snp_summary_stats/bilirubin/"
         "phenocode-TBil_GWAS_in_BBJ_autosome.tsv"),
        usecols=(0, 1, 6),
        skiprows=1,
        delimiter='\t'
    )
    snp_ss_description['total_bilirubin'] = 'Ancestry: Japanese, n=110,000'
    for phenotype in snp_summary_stats:
        snp_summary_stats[phenotype] = prep_data(snp_summary_stats[phenotype])
    
    # Load the NHGRI-EBI GWAS catalog
    catalog = np.loadtxt(
        f'{ukb}/misc_data/snp_summary_stats/catalog/catalog.tsv',
        usecols=[11, 12, 27],
        delimiter='\t',
        skiprows=1,
        dtype=object
    )
    # omit results that are not mapped to a recognizable chromosome
    catalog_names = np.loadtxt(
        f'{ukb}/misc_data/snp_summary_stats/catalog/catalog.tsv',
        usecols=[7],
        delimiter='\t',
        skiprows=1,
        dtype=object
    )
    filter_weird_chroms = np.isin(catalog[:, 0],
                                  list(str(chrom) for chrom in range(1, 23)))
    catalog_names = catalog_names[filter_weird_chroms]
    catalog = catalog[filter_weird_chroms, :]
    catalog = catalog.astype(float)
    #'''

    known_assocs = {}
    #'''
    known_assocs['height'] = catalog[np.equal(catalog_names, 'Height'), :]
    bil_catalog_trait_names = [
        'Total bilirubin levels',
        'Bilirubin levels'
    ]
    known_assocs['total_bilirubin'] = \
            catalog[np.isin(catalog_names, bil_catalog_trait_names), :]
    for phenotype in known_assocs:
        known_assocs[phenotype] = prep_data(known_assocs[phenotype])
    #'''

    if not regions:
        for cutoff in {True, False}:
            make_qq_plots(phenotypes, results, snp_summary_stats,
                          args.run_name, cutoff)
            make_genome_manhattan_plots(
                phenotypes,
                results,
                snp_summary_stats,
                snp_ss_description,
                known_assocs,
                args.run_name,
                cutoff
            )
    else:
        for region in regions:
            make_region_plot(
                results,
                snp_summary_stats,
                snp_ss_description,
                known_assocs,
                args.run_name,
                args.imputation_run_name,
                region
            )

if __name__ == "__main__":
    main()

