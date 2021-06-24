===================
Association Testing
===================

Preparing Phenotypes/Covariates
===============================

#. Loading shared covariates

   * Shared covariates are standardized (subtract mean, then divide by standard deviation)
     Note: this is done before sample subsetting - but this is good enough
     to achieve numeric stability, which was the intent.

   .. details:: sanity checks

       Confirmed that the sex encoding is male == 1 and female == 2 by comparing
       ``$UKB/microarray/ukb46122_cal_chr1_v2_s488282.fam`` to 
       ``$UKB/main_dataset/extracted_data/sex.tab`` and noting that
       `this <https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=9>`_ is the
       encoding for sex in that file.

#. Loading phenotypes and phenotype-specific covariates
  
  * For phenotypes which were recorded on multiple visits, only the first
    recording is used
  * An indicator covaraite is addded for each visit beyond the first where 
    recordings were taken from for any participant
  * Except: visits with number of first recordings less than 0.1% of the
    total participants, or less than 50, are excluded. The participants whose
    first recordings are in that visit are dropped from the sample set
    for this phenotype.
  * Ages are loaded corresponding to the visit where each phenotype was recorded
  * Subsequently, the full sample set is filtered to only contain white brits
    which passed the UKB's published quality filters, without sex aneuploidy
    or sex mismatch, and have not withdrawn themselves
  * Subsequent to that, the full sample set is filtered to contain only unrelated
    samples
  * Lastly, for the remaining samples, the phenotype is rank inverse normalized, 
    and all covariates are standardized to mean zero and variance 1.

  .. details: standardization

    Standardization is to help with numerical stability, and is required
    by plink for that reason in some cases. For effect sizes to be comparable 
    between plink SNP and my STR association results, we standardize the variables
    before any testing instead of one-off during testing.
    The only theoretical drawback to standardization is the need to unstandardize
    effect sizes for them to be interpretable. However, we are rank inverse
    normalizing the phenotype anyway, so effect sizes are not interpretable regardless.
    Aside from the shift in scale, the results should be identical.
    See `here <https://groups.google.com/g/plink2-users/c/midmoPgUntA>`_
    for plink's author agreeing that there is no issue with standardization.

Loading And Filtering Genotypes
===============================

.. details:: Sanity checks

   Confirmed that bgen sample order is same as STR VCFs. (TODO double
   check against this file: ``./ukbgene imp -m -c1 -a../main_dataset/raw_data/k29170.key``)

   File ``$UKB/side_analyses/confirm_sample_lists_same.py``. Contents:

   .. literalinclude:: ../side_analyses/confirm_sample_lists_same.py
       :language: python

STRs
----

.. details:: Sanity checks

   2021/02/17 - manually confirmed that the length allel dosage r2 is correct for (chr1,
   pos 1048570, STR_384) for the first 8 samples

   .. code:: python
       
       hard16 = np.array([0,0,0,1,1,1,1,1,0,1,0,0,0,0,0,0])
       hard15 = 1 - hard16
       prob16 = np.array([.04,0,.01,1,1,1,1,.97,.01,.99,0,0,0.2,0.06,0,0])
       prob15 = 1 - prob15
       prob15[7] = 0
       np.corrcoef(hard16,prob16)[0,1]**2
    
       > 0.989749155123994

       np.corrcoef(hard15,prob15)[0,1]**2

       > 0.9900357942862258

.. details:: Allelic Dosage R2

   Based on Beagle's R2 score
   Appendix 1 here
   https://www.cell.com/ajhg/fulltext/S0002-9297(09)00012-3#app1
   Browning, Brian L., and Sharon R. Browning. "A unified approach to genotype imputation and haplotype-phase inference for large data sets of trios and unrelated individuals." The American Journal of Human Genetics 84.2 (2009): 210-223.

   PerasonCorr(length allelic dosage, hardcall) within subset samples

   Beagle's article talks about PearsonCorr(True genotypes, Hardcalls)
   and then looks at TrueGenotypes|Dosages (assuming dosages are accurate). I'm not sure
   if that's exactly equivalent to this metric. Need to check.


Microarray SNPs
---------------

Imputed SNPs
------------

.. details:: sizing

    Total variants: 93095623
    Number variants per chrom:
    1 7402791
    2 8129063
    3 6696680
    4 6555871
    5 6070641
    6 5751712
    7 5405524
    8 5282223
    9 4066774
    10 4562904
    11 4628348
    12 4431052
    13 3270217
    14 3037521
    15 2767971
    16 3089612
    17 2660711
    18 2599579
    19 2087017
    20 2082571
    21 1261158
    22 1255683

.. details:: Sanity checks

   2021/02/11 - manually confirmed dosage loading in load_imputed_snps is correct. Still need
   to check dosages=False, info_thresh and call_thresh

.. details:: Thoughts on INFO threshold
   
    UKB paper suggests 0.3 :
    https://www.ukbiobank.ac.uk/wp-content/uploads/2014/04/imputation_documentation_May2015.pdf
    Neale lab suggests 0.8:
    http://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas

.. details:: full plink output

   4.4GB for chr21 for one phenotype. ``93095623/1261158 * 4.4GB = 324.8GB`` for an entire phenotype.
   Workable for a few phenotypes, not many, need to work on scaling down. Either filter the files
   and throw away the originals or don't request as much information from plink.

.. details:: Using plink to compute my own metrics on the sample subset

   :code:`--freqs` seems to work. columns :code:`altfreq,alteq,altnumeq` just differ by labeling (either
   <freq> or <allele>=<freq> or <allele_num>=<freq>).

   Comparison of MACH2 imputation metric from :code:`--freq` to IMPUTE v2's INFO metric: 
   This paper:
   Marchini, Jonathan, and Bryan Howie. "Genotype imputation for genome-wide association studies." Nature Reviews Genetics 11.7 (2010): 499-511.
   https://www.nature.com/articles/nrg2796
   (specifically, Figure 1 and supplementary information s3)
   linked to by this documentation https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#info_metric_details
   The INFO metric that plink calculates put forth by Mach2 is apparently very similar to IMPUTEv2
   so should be similar to the info metric given out by UKB. Not super similar to the Beagle info
   metric tho.

   Plink HWE only uses hardcalls (see bullet here: `proper support for dosages`
   https://www.cog-genomics.org/plink/2.0/ ) This is what I do when calculating HWE for 
   strs. However, plink hardcalls implicitly filter many genotypes where dosage isn't close
   enough to a clear hardcall, while I don't do this for STRs. Change?


Association and stats
=====================

.. details:: sanity checks

   2021/02/25 - confirmed that the single dosage means being calculated post linear regression are correct

   .. code:: python

       test_samples = np.isin(data[:, 0], [2497795, 2143467, 1288463, 2632032, 3667876, 3154457, 5713647, 2548437,                                
                 1218644, 3505384,])

       for _len in dosage_gts:
           dosage_gts[_len] = dosage_gts[_len][:10, :]

       mean_stats = statsmodels.stats.weightstats.DescrStatsW(    
           data[test_samples, col_names.index(f'{dep_var}_residual')],
           weights = dosages
       )


.. details:: cost

   Before an str run: 3001

Post association QC
===================

.. details:: First round chr21 plot

   Ways of interacting:

   * Scroll wheel: zoom in/out
   * Click and drag: scroll across the chromosome left/right
   * Use the p-value cap slider above the plot: change the height
     cap of the plot
   * Mouse over: see details about individual loci (if overlapping
     tooltips appear, zoom in enough that the points
     are separated underneath your mouse)
   * Click on a legend: hide that data source
   * Click on the save icon in the toolbar on the right: snap a picture of the plot
   * Click on the +/- zoom icons in the toolbar on the right: also zoom in/out

   Plot details:

   * All STRs are plotted, SNPs are only plotted if they have p values <= 10^-3
   * The 'my code' and 'Plink' SNP runs were both done by me on the same loci and
     same data - they should be close to equivalent
   * In the tooltips for the my code runs, total_hardcall_alleles refers to the
     allele distribution across the entire UKB population, whileas
     subset_total_hardcall_alleles refers to the distribution across only the samples
     on which the regression was run (e.g. qc'ed, unrelated, white, and with 
     a measurement of the phenotype of interest)
   * per_allele_dosage means the sum of the dosages of that allele across all samples
   * allele_dosage_r2 is the squared pearson correlation coefficient between the dosages
     and the hardcalls - numbers closer to 1 indicate that we were more confident in
     calling haplotypes with this length

   Thoughts:

   * The interactivity of the plots seems very useful for exploring the data, and
     now that I know what I'm doing making this again wouldn't be too much work. What
     do you think? Is it useful to you?
   * There's a bug with at least some of the STR tooltip data: locus 3009776 in the
     height plot (the dosage for allele len 14 is zero but there are ~800k hardcalls
     with that allele)
   * Assuming the plots are roughly accurate despite that bug, there don't seem to be
     many obviously spurious STR loci - they follow the trend of the SNP signals and
     published GWAS signals closely. Are there any spots that stand out to you as
     worrisome?
   * There are a few STR loci in the height plot which are significant and have many
     alleles - next step will be to plot them and see if any of them are 
     highly convincing.
   * Other next steps: Extend these plots to the entire genome. Write code to 
     automatically detect STR association signal regions. Read up on Bolt LMM
     and use it for both SNP and STR associations. Do this on more quantitative traits.

   .. raw:: html
      :file: ../association/plots/me_manhattan_height.html

   .. raw:: html
      :file: ../association/plots/me_manhattan_total_bilirubin.html

