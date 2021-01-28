==================================
QC imputation vs exome-called STRs
==================================

Calling STRs with HipSTR
========================

#. Downloaded CRAMs for 40k samples

   .. details:: Text

       Individual-level CRAM files and indicies produced from the UKB Whole Exome
       Sequencing data by the Functionally Equivalent pipeline were downloaded
       from the showcase fields 23163 and 23164. Due to UKB throttling of download
       speeds only 40k of 50k available CRMAs were downloaded
       (the selection of samples to download was random).   

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

   .. details:: Text
           
       Start and end coorindates and repeat periods for STRs in the SNPSTR reference
       panel [SNPSTR]_ were acquired through personal communication with Mr. Saini.
       These fields were arranged into a BED file meeting HipSTR's specification
       [Hipstr_BED_specification]_. :code:`LiftOver` [LiftOver]_ was used to convert
       this BED from hg19 to hg38 coordinates; no STRs were lost during conversion.
       The UKB exome capture region was downloaded from
       `<https://biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/xgen_plus_spikein.b38.bed>`_ .
       :code:`BEDTools intersect` [BEDTools]_ with the flag :code:`-f 1` was used
       to intersect the hg38 SNPSTR BED file and the exome capture BED file,
       producing a BED file of STRs from the SNPSTR reference panel which were
       completely contained in the exome capture region. This contained 630
       out of 445720 STRs

   .. details:: Details

       * Per discussion with Melissa: Exome sequencing generally covers
         the boundaries of the target file and a bit more. If there are
         massive expansions in the target region then this boundary guarantee
         will fail. But small STR expansions should generally be captured by
         exome sequencing, so it's reasonable to include STRs that are at/near
         the boundaries of the capture region. Just subject any use of the STRs
         to standard filtering.

   .. details:: BED contents

      .. csv-table:: SNPSTR Exome STRS hg38 BED
             :file: ../side_analyses/exome_strs/snpstr_exome_strs_38.bed
             :delim: tab

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

   .. details:: Text TODO
       
       :code:`HipSTR` [HipSTR]_ was invoked 
       to call individual loci in batches of 4k individuals at a time, this being
       close to the maximum number of open files permitted per process in the computational
       environment. Batching of individuals was random. In addition to the 
       standard flags, ``HipSTR`` was run with the flags ``--bam-libs CRAM_FILES --bam-samps CRAM_FIELS
       --regions BED --max-reads 10000000`` where ``CRAM_FILES`` is a comma separated list
       of paths to the CRAM files included in the batch and BED is a path to
       the hg38 BED of exome STRs produced above. ``--max-reads`` was set to this limit
       to effectively remote the arbitrary depth cap that HipSTR imposes. HipSTR successfully
       called 587 of the 630 loci. MergeSTR from
       TRTools [TRTools]_ was
       used with the ``--trim`` flag to combine the VCFs for each batch at each locus
       into a single VCF for each locus, and BCFtools [BCFtools]_ was used to concatenate
       the loci VCFs into a single exome VCF.

   .. details:: TODO Details

       * Batching was done in order of increasing sample number,
         as sample numbering is random.

   .. details:: TODO Implementation Details

       * Filtered IDs: ``$UKB/exome/fe_cram_str_calls/filtered_ids.txt``

         .. details:: Contents

             .. literalinclude:: ../exome/fe_cram_str_calls/filtered_ids.txt

       * Calling mergeSTR ``$UKB/exome/fe_cram_str_calls/call_mergeSTR.py``
        
         .. details:: Code 

             .. literalinclude:: ../exome/fe_cram_str_calls/call_mergeSTR.py
                 :language: python

       * Concatenating with bcftools

         .. details:: Code

             .. code:: bash

                 cd $UKB/exome/fe_cram_str_calls
                 awk '{ print "merged_vcfs/" $1 ".merged.vcf" }' unfiltered_ids.txt > merged_filenames.txt
                 bcftools concat -f merged_filenames.txt -o exome_calls -O z --threads 1


Filtering HipSTR calls
======================

#. Pre filtering figures

   .. details:: Figures

       .. image:: ../exome/fe_cram_str_calls/filtering/pre_filtering_q-chrom-callnum.png

       .. image:: ../exome/fe_cram_str_calls/filtering/pre_filtering_q-sample-callrate.png

       .. image:: ../exome/fe_cram_str_calls/filtering/pre_filtering_q-diffref-histogram.png

       .. image:: ../exome/fe_cram_str_calls/filtering/pre_filtering_q-diffref-bias.png

       .. image:: ../exome/fe_cram_str_calls/filtering/pre_filtering_q-quality-per-call.png

       .. image:: ../exome/fe_cram_str_calls/filtering/pre_filtering_q-quality-per-locus.png

       .. image:: ../exome/fe_cram_str_calls/filtering/pre_filtering_q-quality-per-sample.png

   .. details:: TODO Thoughts

   .. details:: Code

       .. code:: bash

           cd $UKB/exome/fe_cram_str_calls/filtering
           ./call_qcSTR.sh \
               $(pwd)/../exome_calls.vcf.gz \
               $(pwd)/../pre_filtering_q

       .. literalinclude:: ../exome/fe_cram_str_calls/filtering/call_qcSTR.sh
           :language: bash

#. Running DumpSTR

   .. details:: Filtering Reports

       Per locus ``$UKB/exome/fe_cram_str_calls/filtering/filtered_exome_calls.loclog.tab``
       
       .. literalinclude:: ../exome/fe_cram_str_calls/filtering/filtered_exome_calls.loclog.tab

       Per sample (head) ``$UKB/exome/fe_cram_str_calls/filtering/filtered_exome_calls.samplog.tab``

       .. literalinclude:: ../exome/fe_cram_str_calls/filtering/filtered_exome_calls.samplog.tab
           :lines: 1-3

   .. details:: Code

       ``$UKB/exome/fe_cram_str_calls/filtering/call_dumpSTR.sh`` 

       .. literalinclude:: ../exome/fe_cram_str_calls/filtering/call_dumpSTR.sh
           :language: bash

#. Post filtering figures

   .. details:: Figures

       .. image:: ../exome/fe_cram_str_calls/filtering/post_filtering_q-chrom-callnum.png

       .. image:: ../exome/fe_cram_str_calls/filtering/post_filtering_q-sample-callrate.png

       .. image:: ../exome/fe_cram_str_calls/filtering/post_filtering_q-diffref-histogram.png

       .. image:: ../exome/fe_cram_str_calls/filtering/post_filtering_q-diffref-bias.png

       .. image:: ../exome/fe_cram_str_calls/filtering/post_filtering_q-quality-per-call.png

       .. image:: ../exome/fe_cram_str_calls/filtering/post_filtering_q-quality-per-locus.png

       .. image:: ../exome/fe_cram_str_calls/filtering/post_filtering_q-quality-per-sample.png

   .. details:: TODO Thoughts

   .. details:: Code 

       .. code:: bash

           cd $UKB/exome/fe_cram_str_calls/filtering
           ./call_qcSTR.sh \
               $(pwd)/../exome_calls.vcf.gz \
               $(pwd)/../pre_filtering_q


Comparing HipSTR calls to Imputed Calls
=======================================


