===================
Association Testing
===================

Preparing Phenotypes
====================

#. Regressing out covariates

   .. details:: text

       Using the CSAPS [CSAPS]_ implementation of the approximation
       of data with cubic smoothing splines algorithm [CSAPS_algorithm]_

   .. details:: Sanity checks

       * Confirmed that csaps is deterministic and fast
         ``$UKB/association/time_smoothing_spline.py``

         .. details:: code

             .. literalinclude:: ../association/time_smoothing_spline.py
                 :language: python
