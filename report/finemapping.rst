===========
Finemapping
===========

Set up 
======

#. Generate 'clumps'/'signal regions'

   .. details:: Text

       Takes my STR and plink SNP results, centers a 250kb interval on each
       variant that passes the GWAS threshold  (125kb in either direction),
       and then merges all overlapping intervals. Similar to plink 1.9's
       clumping algorithm but simpler.

       Filters the extended MHC complex from choromosome 6 (25mb - 33.5mb) as it
       is more properly handled by specialized finemapping tools.

   .. details:: Implementation details

       * Didn't use plink for clumping because plink2 hasn't implemented the clumping
         feature yet and plink1.9 can't read from bgen files.

       * Workflow step :code:`clump_results`
   
       * Script at :code:`$UKB/finemapping/clump.py`
 
         .. details:: Contents

             .. literalinclude:: ../finemapping/clump.py
                 :linenos:
                 :language: python

       * Outputs to :code:`$UKB/finemapping/signal_clumps/`

#. Filtered SNPs that accounted for STR alleles

   .. details:: Text

       If two identical variants are included in a finemapping run, then the finemapper
       has no choice but to assign half the combined probability to each of them which dilutes
       its conclusions. For a multiallelic STR, having a SNP which represents one of STR's alternate
       alleles similarly dilutes the finemapping, albeit not as much as a completely duplicate variant.
       So before finemapping we've filtered out some of the UKB's imputed SNPs that overlap STRs we are
       testing so as to clearly demonstrate where causality is or is not contained by STRs.
       We do so conservatively; we only filter out SNPs that:

       * are near only a single STR and not a region which has two STRs back to back
       * are a clear insertion (ref len = 1) or clear deletion (alt len == 1)
       * are wholly contained within the STR or insert bps between the edges of STR and the first
         flanking bps
       * the derived repeat unit for the STR is at least twice as common as all possible alternative
         choices of repeat units at that length
       * insert full or partial copies of that repeat unit with no impurities

       There approach is conservative as there are other possible representations of indels that
       this simple paradigm would not catch, it only filters pure repeat indels when SNPs containing
       impure indels still have the possibility of confounding some of the signal of the STR they overlap,
       and it does not attempt to filter SNPs in complex regions with multiple STRs, among other limitations.
       A broader filtering choice would be to simply filter any SNPs which overlap an STR being studied.
       Still, this approach also has the potential to filter some SNPs that possibly should be kept
       in the finemapping analysis, such as SNPs which insert or delete partial repeats in the middle
       of an STR which might be more properly considered impurities, and SNPs which insert
       misaligned rotations of the repeat unit (e.g. inserting CAG at the beginning of an AGCAGCAGC repeat).

       This results in filtering 278,405 variants.
   
   .. details:: Implementation details 

       * Script at :code:`$UKB/finemapping/str_imp_snp_overlap.py/`
 
         .. details:: Contents

             .. literalinclude:: ../finemapping/str_imp_snp_overlap.py
                 :linenos:
                 :language: python

       * Outputs to :code:`$UKB/finemapping/str_imp_snp_overlaps/`



Running Finemapping Tools 
=========================

#. FINEMAP

   .. details:: Details

       Ran FINEMAP on all variants in each clump with p-value >= 0.05.
       Used default arguments, except :code:`--n-causal-snps 20`

