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

   .. details:: Methodology

       * Cull the STRs from the SNPSTR reference VCFs and make an hg19 BED file
       * LiftOver the BED to hg38 - nothing fails to lift over
       * Intersect the new hg38 BED with UKB's target capture region
         ( https://biobank.ctsu.ox.ac.uk/showcase/refer.cgi?id=3801 )

   .. details:: Details

       * 630 exome STRs
       * Per discussion with Melissa: Exome sequencing generally covers
         the boundaries of the target file and a bit more. If there are
         massive expansions in the target region then this boundary guarantee
         will fail. But small STR expansions should generally be captured by
         exome sequencing, so it's reasonable to include STRs that are at/near
         the boundaries of the capture region. Just subject any use of the STRs
         to standard filtering.


   .. details:: Implementation details

       * File at: :code:`$UKB/side_analyses/exome_strs/snpstr_exome_strs_38.bed`
         Produced via

         .. details:: code

             .. code-block:: bash

                 python make_snpstr_bed.py
                 #bedPlus argument means that the fields 4+ get carried over directly instead of being interpreted as standard bed fields
                 $UKB/utilities/liftOver/liftOver -bedPlus=3 snpstr_strs_19.bed $UKB/utilities/liftOver/hg19ToHg38.over.chain.gz snpstr_strs_38_unsorted.bed unmapped_19to38_snpstr_strs.bed
                 for chr in $(seq 1 22); do echo "chr$chr" >> chr.names ; done
                 bedtools sort -i snpstr_strs_38_unsorted.bed -g chr.name > snpstr_strs_38.bed

                 # see here https://biobank.ctsu.ox.ac.uk/showcase/refer.cgi?id=3801
                 wget  -nd  biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/xgen_plus_spikein.b38.bed ; mv xgen_plus_spikein.b38.bed exome_38.bed
                 awk '{ print "chr" $0 ; }' exome_38.bed > exome_38_chr.bed

                 bedtools intersect -u -f 1 -a snpstr_strs_38.bed -b exome_38_chr.bed > snpstr_exome_strs_38.bed

       * Helper script :code:`make_snpstr_bed.py`

         .. details:: code

             .. literalinclude:: ../side_analyses/exome_strs/make_snpstr_bed.py
                 :language: python

       * Overlapping checks:

          * STR_2876 isn't included because it only partially overlaps the exome region
          * While STR_1214 is included because it fully overalps the exome region

       * SNPSTR IDs do not match up with HipSTR reference IDs. I cannot find the
         script I wrote which showed this, but I'm confident of this.
         Shubham says they matched initially when SNPSTR ref was published
         and then HipSTR changed the IDs out from under him.
         
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

   .. details:: Details

   .. details:: Implementation Details


Filtering HipSTR calls
======================

Comparing HipSTR calls to Imputed Calls
=======================================
