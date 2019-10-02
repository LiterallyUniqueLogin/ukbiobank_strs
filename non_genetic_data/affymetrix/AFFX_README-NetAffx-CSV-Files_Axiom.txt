README for NetAffx CSV annotation files for Axiom Genotyping Products.
 
Copyright 2008-2011, Affymetrix Inc.
All Rights Reserved

The contents of this CSV file are covered by the terms of use or 
license located at http://www.affymetrix.com/site/terms.affx

Array name:     Axiom product line arrays
Organism:       Non-specific

This README provides a guide to the contents of the CSV files 
containing annotations for the Affymetrix Axiom Genotyping arrays. 

Contents
--------

I. General Notes
   A. Informative Header
   B. File format conventions
   C. Annotation overview

II. CSV file: Column Descriptions

III. References

I. General Notes
-----------------

I.A. Informative header

   The CSV file contains an informative header consisting of a set of
   key=value attributes describing the data within the file. This
   includes details about the array design and the NetAffx release on
   which the data in the file are based. The keys do not contain
   spaces, but the values may.

   For compatibility with Genotyping Console (GTC) and Affymetrix 
   Power Tools (APT) software, each attribute-containing line 
   of the informative header starts with '#%' while descriptive comment lines
   of the header start with '##'.

   The first non-commented row of the CSV file is a column header,
   containing the primary column titles.

   Description of header keys:

   General information:
      create_date     - date the CSV file was created
      chip_type       - official designation of the chip (a.k.a. array_type)
      lib_set_name    - array name portion of library file name
      lib_set_version - array version portion of library file name

   Genome information: 
      genome-species                - genus and species of organism
      genome-version                - version used during array annotation
      genome-version-ucsc           - UCSC version tag for genome assembly
      genome-version-ncbi           - NCBI version tag for genome assembly
      genome-version-create_date    - genome assembly release date

   NetAffx information:
      netaffx-annotation-date                    - NetAffx release date
      netaffx-annotation-netaffx-build           - NetAffx release number
      netaffx-annotation-tabular-format-version  - version of the file format
      netaffx-annotation-docgen-method           - specific class of the "docgen" software used to generate the annotation file
      netaffx-annotation-docgen-version          - version of the "docgen" software used to generate the annotation file

   Source database information:
      dbsnp-date                    - date of dbSNP version release
      dbsnp-version                 - dbSNP version used during array annotation
      dgv-date                      - date of the "Database of Genomic Variants" version release
      dgv-version                   - "Database of Genomic Variants" release used for array annotation
      hapmap-date                   - date of HapMap version release
      hapmap-version                - HapMap version used during array annotation
 
I.B. File format conventions

   The CSV files follow the conventions of other NetAffx tabular data 
   files as described at:   
 
   http://www.affymetrix.com/support/technical/manual/taf_manual.affx
 
   Missing values in the CSV annotations are represented as "---".
 
I.C. Annotation overview
   
   SNPs selected for Affymetrix arrays are from two sources:
   dbSNP (http://www.ncbi.nlm.nih.gov/projects/SNP/) and the 1000 Genomes
   Project (http://www.1000genomes.org/).

   1. Chromosomal positions.

   For SNPs in dbSNP and the 1000 Genomes Project, SNP chromosomal positions 
   are obtained from the dbSNP database and 1000 genomes data sources. 
   On rare occasions, a SNP selected at design time does not uniquely map to 
   the reference genome at time of annotation update. No chromosome positions 
   are reported in those cases.
   
   Base positions are reported by 1-base index. A single base location is reported
   for each SNP, regardless of the actual length of the polymorphism. For
   single-nucleotide replacement SNPs (i.e. single base transitions and
   transversions) the position reported is the reference-sequence base that is
   affected by the polymorphism.  For indel alleles (insertions and deletions
   with respect to the reference sequence), the position depends upon whether
   the reference genome contains the insertion or the deletion. 
   
   If the reference genome represents the single base insertion (i.e. the indel
   base is present in the genomic reference sequence), then the base reported
   is the base that is affected by the polymorphism (i.e. the base that is
   deleted in the non-reference case). 
   If the insertion (present in the genomic reference sequence) is multi-base
   (i.e. >1 base long), then the position reported is the first base [5' base]
   of the insertion allele on the plus strand of the reference genome.
   
   If the reference genome represents the deletion (i.e. the indel bases are
   not present in the genomic reference sequence) then the position reported
   is the base to the left [5'] of the insertion site on the plus strand of
   the reference genome. This is true for both single and multi-base
   deletions. i.e., the insertion site is the inter-base location immediately
   3' of the reported base position on the plus strand of the reference genome.
 
   2. Relationship to surrounding genes.

   SNP chromosomal positions were compared with the annotated gene structure 
   of Ensembl and RefSeq genes, and the relationship is annotated 
   as UTR-5, CDS, UTR-3, and intron when the SNP falls into regions of 5'UTR, 
   CDS, 3'UTR, and intron respectively. When a SNP is not within any genes, 
   it is annotated as upstream (when the SNP is upstream of the 5' end of 
   the transcript, on the transcript's strand) or downstream (when the SNP is 
   downstream of the 3' end of the transcript, on the transcript's strand) 
   relative to the neighboring transcripts. For upstream and downstream SNPs, 
   the distance from the SNP to the gene is reported. When a SNP is within an exon, 
   but there is no reported CDS for the transcript, then the value 'exon' is reported.
   When a SNP is within 2 bases from a splice site, then the value
   'splice-site' is reported. When a SNP is within a CDS but does not result in a 
   residue change, 'synon' is reported. When a SNP is within a CDS and may result
   in a residue change, 'missense' is reported. When a SNP is within a CDS and may
   result in a premature stop codon, 'nonsense' is reported. When 
   an indel is not a factor of 3 and may cause a frameshift, 'frameshift' is reported.
   When an indel is a factor of 3 and may cause an amino acid insert or deletion,
   'aa-indel' is reported.
   Only one value is reported per SNP-transcript relationship, according to the 
   following hierarchy:
    splice-site
    intron
    nonsense
    missense
    synon
    frameshift
    aa-indel
    CDS
    UTR-3
    UTR-5
    exon
    upstream 
    downstream
   e.g. Since 'splice-site' and 'intron' are not mutually exclusive, only 
   'splice-site' is reported. For each transcript data source (Ensembl 
   and RefSeq), all genic (i.e. non-upstream, non-downstream) relationships
   for all transcripts are reported. If a genic relationship is reportable 
   for a given data source, then upstream and downstream relationships are
   not also reported. For upstream and downstream relationships, only the
   2 closest transcripts 5' and 3' on the genomic sequence are reported,
   for each transcript data source. (Thus, it is possible to be downstream
   of 2 transcripts from the same data source, one 5' to the SNP, and one 
   3' on the genomic sequence.)

   3. Genetic maps.

   This annotation is to provide a rough estimate on SNP genetic
   distances to the p-telomere. Those estimates may be used as seed
   input for linkage analysis programs like MERLIN
   (http://www.sph.umich.edu/csg/abecasis/Merlin/). As a requirement
   of MERLIN, every SNP has to have a unique genetic distance. SNP
   genetic distances were extrapolated from three experimentally
   obtained genetic maps: deCODE map, Marshfield map, and SLM1 map.  

      1). deCODE genetic map was built by genotyping 5,136
          microsatellite markers for 146 families (1), and is available
          through Nature Genetics (see reference).  

      2). Marshfield genetic map (2) was based on CEPH family genotypes
          for 7,740 microsatellite markers, and is available from
          http://research.marshfieldclinic.org/genetics/Map_Markers/maps/IndexMapFrames.html 

      3). SLM1 (SNP Linkage Map) map was generated from unpublished data
          from Affymetrix and Dr. Aravinda Chakravarti group at Johns
          Hopkins University. It was based on genotypes for 2,022
          microsatellite markers  and 6,205 SNPs.   

   Physical locations of markers used in each genetic map were
   obtained from the UCSC database (ftp://genome.ucsc.edu). Markers
   are removed when their genetic order is opposite of their physical
   order. When physically neighboring markers share the same genetic
   distance, only the one with the largest physical position was kept
   in order that no two SNPs have exactly the same genetic
   distance. Physical positions of SNPs and physical locations of the
   markers in the cleaned genetic maps were compared to infer genetic
   distances for SNPs. We assume that genetic distance changes
   linearly with physical distance between any two neighboring markers
   in each genetic map. 
 
 
II. CSV file: Column Descriptions
---------------------------------

   Note: Missing annotation values are represented as "---".

   1. Probe Set ID
        
      Unique identifier for the probe set.
        
   2. dbSNP RS ID
        
      SNP identifier from dbSNP. 
        
   3. Chromosome (1-22, X, Y, MT, for Human.)

      Other values allowed for non-human, e.g. >22 (for bovine etc), 
      W, Z (for avian sex chromosomes), C (for Chloroplast).
      See also "Annotation overview".

   4. Physical Position 

      1-based SNP chromosomal position. See also "Annotation 
      overview". A positive integer, unless the SNP is unmapped, 
      in which case the value is "---".
        
   5. Strand (+|-)
        
      The strand (of the reference genome) where the SNP is 
      mapped. This information may not be the same as that of 
      dbSNP for the same SNP. For the Axiom product line, the 
      annotations are fixed on the plus strand of the reference genome.
      Thus, for Axiom arrays, this value should always be "+", unless the
      SNP is unmapped, in which case the value is "---". 
      See also "Allele A", wherein strand impacts the abstract allele naming.
   
   6. ChrX pseudo-autosomal region 1 (PAR1)
 
      Value is 1 if the SNP is located in PAR1; otherwise value is 0.
      This p-arm pseudo-autosomal region has the following chromosomal
      coordinates in humans: (may not be available or annotated in other species)
      
      per UCSC, hg18: (human PAR1 and PAR2)
      chrY:1-2709520 and chrY:57443438-57772954
      chrX:1-2709520 and chrX:154584238-154913754 
      http://genome.ucsc.edu/cgi-bin/hgGateway?clade=mammal&org=Human&db=hg18
     
      per UCSC, hg19: (human PAR1 and PAR2)
      chrY:10001-2649520 and chrY:59034050-59363566
      chrX:60001-2699520 and chrX:154931044-155260560
      http://genome.ucsc.edu/cgi-bin/hgGateway?clade=mammal&org=Human&db=hg19 
 
   7. Cytoband

      The chromosome band seen on Giemsa-stained chromosomes.
      Value is the band number where the SNP probe set maps.
 
   8. Flank
        
      Allele and Flanking probe sequence of the SNP, usually 71 nucleotides
      or more in size. The two alleles of the SNP are placed within square brackets.
      All sequences are listed 5' to 3' on the forward strand of the reference genome. 
      Alleles are reported as on the forward strand of the genome.  
  
   9. Allele A
        
      At array (or underlying database) design time, Affymetrix follows the following
      naming convention to assign allele nucleotide bases to the "Abstract" allele 
      codes "A" and "B":
      
      1. SNPs are fixed on the forward strand of the design-time reference genome;      
      2. For AT or CG SNPs (SNP alleles are A/T or C/G), the alleles are named in
         alphabetical order (A and C are the "A" alleles, in these cases);      
      3. For non-AT and non-CG SNPs, allele A is A or T, allele B is C or G;      
      4. For indels, allele A is "-", allele B is the insertion.      
      5. For multi-base alleles, the alleles are named in alphabetical order.
         (For [AGT/TTA], AGT would be "Allele A". For [GGT/TTA], GGT would be "Allele A". );
      
      The alleles are named according to the design-time reference genome build,
      and are fixed through subsequent genome lifts. 
      So if an allele-specific SNP (A/T or G/C) is lifted to the opposite strand
      on a new genome, and now annotated with respect to that new genome, it will
      now be in violation of the naming convention, had it been originally
      designed on that new genome. 
      This effect is a consequence of the fact that in the 2-colour assay system,
      the complementary A and T bases use the same colour-channel for non-A/T and
      non-C/G SNPs (as do the complementary C and G bases, for the other
      colour-channel). The fact that the 2 bases that share a colour channel are
      *complementary* is relevant. Thus, when a non-A/T and non-C/G SNP changes
      strands (as a result of a local inversion in a new genome re-assembly)
      the allele codes remain obedient to the naming convention. However, for
      A/T and C/G SNPs, the colour-channel cannot distinguish the 2 alleles, so
      separate allele-specific oligonucleotides are used (where the allele is
      the last 5' base of the oligo). In the event of a strand change as above,
      the CDF file still ties a specific allele code (e.g. "B") to a specific
      base (e.g. "T") in the context of its adjacent bases on the specific
      oligo at a specific x,y location on the array, regardless of the strand
      that oligo now maps to. Hence, it will now be in violation of the naming
      convention, had it been originally designed on that new genome.
  
   10. Allele B
        
      See definition in "Allele A" above.
        
   11. Associated Gene
        
      Values are a list of genes which the SNP is associated to 
      (separated by ///). Values for each gene are: "transcript 
      accession // SNP-gene relationship // distance (value 0 if 
      within the gene) // UniGene Cluster ID // gene name or symbol 
      // NCBI Gene ID // GenBank description". The SNP could be 
      within the gene region or be upstream or downstream of the 
      genes. See "Annotation overview" for details.
  
   12. Genetic Map
        
      See "Annotation overview" for details on calculation of genetic 
      distances. Values listed are "sex-averaged genetic
      distance calculated from genetic map // ID of first
      marker used for calculation // ID of second marker used
      for calculation // TSC (the SNP consortium) ID of SNP
      number 1 // TSC ID of SNP number 2 // source" and delimited
      by "///". Note for calculating genetic distance from SLM1, pairs marker number 1/
      marker number 2 or marker number 2/SNP number 1 or SNP number 1/SNP 
      number 2 may be used.
 
   13. Microsatellite
        
      Values are "marker upstream of the SNP // upstream // 
      distance from marker to SNP /// marker downstream of SNP // 
      downstream // distance from SNP to marker" or "marker 
      containing the SNP // within // 0 /// marker upstream of 
      the SNP // upstream // distance from marker to SNP /// 
      marker downstream of SNP // downstream // distance from SNP 
      to marker".
 
   14. Allele Frequencies
        
      Allele frequencies reported here were obtained from 
      Affymetrix' experimental data on different sets of 
      individuals. 
      
      Values are a list of allele frequencies for various populations
      (separated by "///").  Values for each population are 
      "A allele frequency // B allele frequency // Population".
      Populations are sorted alphabetically. If no allele frequencies are 
      available for a given population, then the population will not be 
      reported in the field. Thus, for some SNPs there may be fewer 
      populations reported than for others. The number of populations 
      reported in the following "allele frequency" and "number of 
      individuals" columns is expected to be internally consistent between 
      the currently described column and these subsequent columns.
 
   15. Heterozygous Allele Frequencies
        
      Heterozygosity was calculated from the above allele 
      frequency data. 
      
      Values are a list of heterozygosity frequencies for various 
      populations (separated by "///").  Values for each population
      are "Heterozygosity frequency // Population".
 
   16. Number of individuals
        
      Values are a list of the number of individuals 
      from various populations (separated by "///") used for 
      calculating the above allele frequencies.
      
      For each population the values are "Number of individuals // Population". 
      For the Bovine WG SNP 1.0 Axiom array, the values are numbers of samples.
 
   17. In Hapmap
 
      This indicates if the corresponding SNP is included in the HapMap project.
      Possible values are 'YES' if it is part of HapMap or '---' if either unknown or
      not included in the HapMap.
 
   18. Strand Versus dbSNP
 
      This value indicates the relationship of the annotated strand (always
      the plus strand, in the case of Axiom products) versus dbSNP. 
      Possible values here are 'same', 'reverse' or '---'.  '---' indicates that
      the corresponding SNP is not in dbSNP.
 
   19. Probe Count
 
      The values are the total number of probes in the probeset. Expect default to 
      be 2 for allele-specific oligos, 1 for all others. Else, if values outside [1,2]
      exist, actual probe replicate count is given.
 
   20. ChrX pseudo-autosomal region 2
 
      Value is 1 if the SNP is located in PAR2; otherwise value is 0.
      This q-arm pseudo-autosomal region has the following chromosomal
      coordinates in humans: (may not be available or annotated in other species)
      
      per UCSC, hg18: (human PAR1 and PAR2)
      chrY:1-2709520 and chrY:57443438-57772954
      chrX:1-2709520 and chrX:154584238-154913754 
      http://genome.ucsc.edu/cgi-bin/hgGateway?clade=mammal&org=Human&db=hg18
      
      per UCSC, hg19: (human PAR1 and PAR2)
      chrY:10001-2649520 and chrY:59034050-59363566
      chrX:60001-2699520 and chrX:154931044-155260560
      http://genome.ucsc.edu/cgi-bin/hgGateway?clade=mammal&org=Human&db=hg19 
 
   21. Minor Allele
 
      Values are: minor allele // population. The value is the nucleotide base
      found in either the "Allele A" or "Allele B" field that corresponds to 
      the lowest allele frequency value found in the "Allele Frequencies" field.
      If the values in the "Allele Frequencies" field are equal, then the 
      nucleotide base for "Allele A" is used.
 
   22. Minor Allele Frequency
 
      Values are: minor allele frequency // population.
 
   23. OMIM

      This track furnishes OMIM and Morbid Map IDs and their
      respective gene titles. 
      
      Overlaps of OMIM (TM) genes with the probeset. Values are "OMIM id 
      // Disease title // morbid map id // variant-transcript relationship"

      Disease titles surrounded by square brackets, "[]", indicate nondiseases -- 
      mainly genetic variations that lead to abnormal lab test values. 
      Braces, "{}", indicate mutations that contribute to susceptibility to
      multifactorial disorders (eg. diabetes), or to susceptibility to infection.
      Question mark, "?", indicates an unconfirmed or possibly spurious mapping.
      The number in parentheses after the name of each disorder indicates: 
      (1) the disorder was positioned by mapping of the wildtype gene; 
      (2) the disease phenotype itself was mapped; 
      (3) the molecular basis of the disorder is known; 
      (4) the disorder is a chromosome deletion or duplication syndrome.
      For more information, please read http://omim.org/help/faq#1.6 
 
      This database contains information from the Online Mendelian
      Inheritance in Man (OMIM (TM)) database, which has been obtained under a
      license from the Johns Hopkins University.  This database/product does
      not represent the entire, unmodified OMIM(TM) database, which is
      available in its entirety at www.ncbi.nlm.nih.gov/omim/.

   24. Affy SNP ID

      Unique Affymetrix identifier for the actual SNP (i.e. not the probeset ID).
      Format is Affx-[integer]. The Affx- prefix denotes an Axiom product line 
      ID, whereas a bare integer would designate a WGSA (e.g SNP 6.0) product line ID.

   25. Position End

      1-based SNP chromosomal position. See also "Annotation
      overview". A positive integer, unless the SNP is unmapped,
      in which case the value is "---".

   26. dbSNP Loctype

      Value is 1 if the Insert is in the reference genome. In this case, the single nucleotide 
      on the flanking sequence is substituted with a nucleotide sequence with length greater than one.
      Value is 2 when exactly one nucleotide in the reference genome is replaced. This is a restricted
      case of Loctype 1.
      Value is 3 if the Deletion is in the reference genome. In this case, the SNP site is deleted
      from the reference. 
      More information can be found at http://www.ncbi.nlm.nih.gov/SNP/specs/alignment_types.htm

   27. Annotation Notes

      Format: "descriptor // value" and delimited by "///" where descriptors include the following cases:

      "degenerate" : a variant with (indistinguishable from) at least one other putatively different SNP or indel
      "dbsnp retired" : variant's rsid has retired from this dbSNP build
      "dbsnp merged" : variant's rsid has merged with another
      "strand-flip" : variant has undergone a strand-flip when mapped to this genome build
      "ambiguous position" : variant's genomic position cannot be unambiguously determined

   28. Ref Allele and Alt Allele

      The reference allele and alternative alleles are specified according to the current reference genome build. 
      The value of Ref Allele could be "-", which indicates an insertion after the specified position. 
      Otherwise it is the sequence (one or more bases) of the allele of the marker which matches the current 
      reference genome. The value of Alt Allele could be "-", which indicates that the variant is a deletion with 
      respect to the current genome build. Otherwise it is the sequence (one or more bases) of the allele of the 
      marker which does not match the current reference genome. If neither allele of the marker matches the 
      current genome build sequence, then the value of the Ref Allele is set to "." and the value of the 
      Alt Allele is set to <allele_1>/<allele_2>, where <allele_1> and <allele_2> are the alleles of the marker. 
      If the current genome build position of the marker is unknown then the value of the Ref Allele and 
      Alt Allele are set to "---" to denote missing information.

   29. Biomedical

      Biomedical annotations this marker is a member of.

      Format: "category // value" and delimited by "///" where category indicates the Biomedical entity.


III. References
---------------
 
   1. Kong, A., Gudbjartsson, D.F., Sainz, J., Jonsdottir, G.M.,
      Gudjonsson, S.A., Richardsson, B., Sigurdardottir, S., Barnard,
      J., Hallbeck, B., Masson, G., Shlien, A., Palsson, S.T., Frigge,
      M.L., Thorgeirsson, T.E., Gulcher, J.R., Stefansson, K. (2002) A
      high-resolution recombination map of the human genome. Nat
      Genet., 31, 241-247. 

   2. Broman, K.W., Murray, J.C., Sheffield, V.C., White, R.L., Weber,
      J.L. (1998) Comprehensive human genetic maps: individual and
      sex-specific variation in recombination. Am J Hum Genet., 63,
      861-869.  

   3. This database/product contains information obtained from the
      Online Mendelian Inheritance in Man (OMIM (TM)) database, which has
      been obtained through a license from the Johns Hopkins University.
      This database/product does not represent the entire, unmodified OMIM(TM) 
      database, which is available in its entirety at www.ncbi.nlm.nih.gov/omim/.
      Online Mendelian Inheritance in Man, OMIM (TM). McKusick-Nathans 
      Institute of Genetic Medicine, Johns Hopkins University
      (Baltimore, MD) and National Center for Biotechnology
      Information, National Library of Medicine (Bethesda, MD).
