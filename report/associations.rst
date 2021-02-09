===================
Association Testing
===================

Preparing Phenotypes
====================

#. Regressing out covariates

   .. details:: Sanity checks

       2021/02/08 - checked that for height and bilirubin in the get_residuals_linear
       method that the covariates are being properly loaded by comparing to the
       input files.

   .. details:: cspas (no longer using)

       When running csaps, as I added features the RMSE of the residuals
       would increase. (Those features were PCs, presumably PCs with no
       significant association. Still, should remain near constant, not
       increase). Not sure why, couldn't fix this. So abandoned.
   
       .. details:: text

           Using the CSAPS [CSAPS]_ implementation of the approximation
           of data with cubic smoothing splines algorithm [CSAPS_algorithm]_

       .. details:: Sanity checks

           + Confirmed that csaps is deterministic and fast
             ``$UKB/association/time_smoothing_spline.py``

             .. details:: code

                 .. literalinclude:: ../association/time_smoothing_spline.py
                     :language: python


