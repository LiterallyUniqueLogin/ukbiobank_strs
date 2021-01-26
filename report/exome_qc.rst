==================================
QC imputation vs exome-called STRs
==================================

Calling STRs with HipSTR
========================

#. Downloaded CRAMs for 40k samples (field IDs 23163 and 23164)

   .. details:: Details

       * Chose FE instead of SPB CRAMs because the SPB pipeline
         had issues and was down. Didn't use OQFE pipeline because it was
         not released at the time. 
       * Only downloaded 40k of the 50k files because that took me 2 months.
         Also why I haven't migrated to the OQFE files.
       * Files were downloaded in semirandom order, I believe prioritizing
         lower IDs numerically, but ID ordering is random so this is okay

   .. details:: Implementation details 

       * Files at :code:`$UKB/exome/fe_crams/`
       * Filenames: :code:`$ID_23163_0_0.cram, $ID_23163_0_0.cram.crai`
         and symlinked shorter versions of these names: :code:`$ID.cram,
         $ID.cram.crai`
       * Downloading was limited by max simultaneous downloads (10, set by UKB) and
         bandwidth (I believe also throttled by UKB)
       * Downloaded files on tscc-dm1 per instructions, supposedly faster
       * Downloaded with dask for parallelism. Script at
         :code:`$UKB/exome/fetch_fe_crams_local.py`.

         .. details:: Contents

             .. literalinclude:: ../exome/fetch_fe_crams_local.py
                 :linenos:
                 :language: python

       * Batch file at
         :code:`$UKB/exome/fe_cram.bulk`

         .. details:: Contents (truncated)

             .. literalinclude:: ../exome/fe_cram.bulk
                 :lines: 1-10



#. Build hg38 exome STR locus BED

   .. details:: BED contents

      .. csv-table:: SNPSTR Exome STRS hg38 BED
             :file: ../side_analyses/exome_strs/snpstr_exome_strs_38.bed
             :delim: tab

   .. details:: Details

   .. details:: Implementation details

       * File at: :code:`$UKB/side_analyses/exome_strs/snpstr_exome_strs_38.bed`
         Produced via

         .. details:: code

             .. code-block:: bash

                 python make_snpstr_bed.py

       * Helper script :code:`make_snpstr_bed.py`

         .. literalinclude:: ../side_analyses/exome_strs/make_snpstr_bed.py
             :language: python

       * Also produced :code:`snpstr_exome_str_ids.txt, snpstr_exome_str_cleaned_ids.txt`
         Just the list of IDs, with the latter having the '/' in some STR names replaced by '_'
       * Also produced files containing the calls for those STRs in the SNPSTR hg19 panel
         :code:`snpstr_panel_exome_calls_hg19_chrNN.txt`.

         .. details:: Code

             .. code-block:: bash

                for chr in $(seq 1 22) ; do
                     bcftools query -i ID=@$UKB/side_analyses/exome_strs/snpstr_exome_str_ids.txt \
                         -f '%ID %REF %ALT [%GT:]\n' \
                         $UKB/snpstr/1kg.snp.str.chr${chr}.vcf.gz \
                         > snpstr_panel_exome_calls_hg19_chr${chr}.txt &
                 done 

             
   #. Ran HipSTR


Filtering HipSTR calls
======================
