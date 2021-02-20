===================
Association Testing
===================

Preparing Phenotypes
====================

.. details:: Files
   

#. Regressing out covariates

   For the phenotype height and total bilirubin I have included three categories of covariates:
   sex, age (in years) at time of measurement, the first 40 principal components of ancestry.
   I've attempted to regress out those covariates nonlinearly with four different ML models.
   I've concluded that regressing out covariates linearly is acceptable for now
   and potentially optimal. The nonlinear ML models take some time to optimize and in my
   attempts where I've rushed that process they are all doing slightly worse than the
   linear model. It's possible that when properly optimized they will do better than the
   linear model, but my intution is that even if that is so the gains will be minor.
   I still plan to go back to this at some point in the future and test that, as a training
   exercise for myself in applying different ML models and to learn if there are gains to be had
   with this type of data. I'll probably wait till I have access to Expanse where running
   large jobs will be faster and easier.
  
   .. details:: Thoughts on why more complex models aren't easily outperforming linear models

       (Including discussions with Nolan, a friend)

       The benefit of nonlinear models is "proportional" to the amount of nonlinearity
       in the true relationship between the covariates and the phenotype as compared to the
       the amount of noise that obscures the relationship (linear or othwerise).
       Biology isn't linear, there are definitely nonlinear relations to find here.
       But if their magnitude is overwhelmed by the magnitude of the noise, then even with a large
       dataset we're liable to overfit and have to carefully tune any model if we want to make
       any progress.

       Linear models can also be very good at modelling the relationship of binary predictors
       with a phenotype (they are perfect if the binary predictors are additive with respect to
       one another). This means the linear model is already likely accounting well for sex.
       It's unclear what the shape of the principal component data is: those are nominally
       fourty continuous axes, but if all they do is demarcate separate clusters of the data
       then its possible that the data is close to binary underneath. Again, a linear model
       is good for that sort of data. The only data that a nonlinear model seems clearly better
       for is age.

       Lots of data leads to slower training leads to less well trained models on fixed time
       budgets leads to worse performing nonlinear models.

   .. details:: Attempts so far

       I've worked with cubic smoothing splines implemented by this [CSAPS]_ package with an
       unknown author. I've also worked with three models from scikit learn:

       * kernel ridge regression (very similar to support vector machine regression) with the
         rbf kernel
       * Random forests
       * AdaBoosted decision stumps.

       I've been measuring success as minimizing RMSE: sqrt(mean((predictions - actual_values)**2))
       For reference, the standard deviation of untransformed height is 9.24cm, the validation RMSE
       for the linear model (309k (90%) training samples, 34k (10%) validation samples, averaged over 5 runs)
       is 6.29cm. The std of log(total bilirubin) is 0.391 and the validation RMSE for the linear
       model (294k (90%) training samples, 33k (10%) validation samples, averaged over 5 runs) is
       0.37 log(umol/L).

       .. details:: linear model

           .. literalinclude:: linear_rmse_README

       .. details:: cubic smoothing splines with CSAPS [CSAPS]_ (no longer using)

           This algorithm had two limitations. First, when handling multiple input covariates
           it needed them to be measured on a grid. This isn't the case with real data.
           So instead of giving it multiple covariates, I repeatedly used it to regress
           out a single covariate at a time, with the intent of reducing the residuals
           slightly each step. The second was that it needed each data point to be
           unique. Again, real data with finite percision doesn't do this, so instead
           I added a small amount of noise to all the data. Because there were so
           many data points, I needed to add 5e-4 * [-0.5, 0.5) to the data to get
           each coordinate unique.

           When running csaps, as I added features the RMSE of the residuals
           would increase. (Those features were PCs, presumably PCs with no
           significant association. Still, should remain near constant, not
           increase). Not sure why, couldn't fix this. So abandoned.
           Plausible explanations: 

           * There's a bug in this package, its not from a known source.
           * The jitter is adding too much noise relative to the signal
           * The smoothing parameters I chose from weren't fine enough
             ([0, 1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6,
             3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1],
             1 minus those values, and 0.5)

           .. details:: Sanity checks

               - Confirmed that csaps is deterministic and fast
                 ``$UKB/association/time_smoothing_spline.py``

                 .. details:: code

                     .. literalinclude:: ../association/time_smoothing_spline.py
                         :language: python

       .. details:: kernel ridge regression

           Not sure why, but the implementation of this memory and time both
           scale quadratically in the number of parameters being fit. So max
           number of training samples that will fit in memory is ~64k (122gb).
           (Time of this is 406 sec).
           Tried with 1.6k training, 400 validation, 5 folds, best RMSE for
           height was 6.49. Need to try with larger sample number. On TSCC
           for 5-fold validation and 40**2 metaparameter grid search that should
           take ~$100. (Param space [10**(i/8) for i in range(-80, -40)])

           Could swap out the rbf kernel for a linear kernel to make sure 
           this properly reproduces the linear model in that case.

       .. details:: random forests

           Using the same 90%/10% train/validation split as with the linear model,
           200 trees with min_samples_leaf = 10 gave height RMSE of 6.317 . This is
           very slow, would want to run with many trees parallelized for each fold.
           200 trees performed better than 50 (6.330) or 100 (6.321) indiciating there
           is room for at least some more improvement.

           Caveat: even if RMSE drops below linear, due to the discontinuities of this
           model some of the residuals may be much worse estimates

       .. details:: AdaBoosted decision stumps

           Same 90%/10% split as linear model. RMSE increases as number of stumps
           increase (50: 6.424, 100: 6.445, 200: 6.527). Overfitting? Maybe would
           need to lower learning rate to make this model applicable.


   .. details:: Sanity checks

       2021/02/08 - checked that for height and bilirubin in the get_residuals_linear
       method that the covariates are being properly loaded by comparing to the
       input files.

       2021/02/11 - checked that ranking is working correctly. Checked that inverse
       normalization corresponds to correct samples' ranks. Checked that inverse
       normalization are correct calculations: compared to normal distribution
       quantile function here: https://planetcalc.com/4986/

       .. code:: bash

           # pull out ranks first, residuals second
           paste <(cut -f57 covars_and_phenotypes.tab  | tail -n+2 | grep -v nan ) \
               <(cut -f55 covars_and_phenotypes.tab | tail -n +2 | grep -v nan) \
               | sort | head -n 10

           # matches sort with just residuals

           cut -f55 covars_and_phenotypes.tab | tail -n +2 | grep -v nan | sort -n | head -n 10

        .. code:: bash

            # pull out inverse normalization first, ranks second
            # show that smallest inverse normalization has rank 0
            paste <(cut -f59 covars_and_phenotypes.tab  | tail -n+2 | grep -v nan ) \
                <(cut -f57 covars_and_phenotypes.tab | tail -n +2 | grep -v nan) \
                | sort -n | head -n 10


Loading And Filtering Genotypes
===============================

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


