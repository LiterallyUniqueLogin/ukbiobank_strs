Any STR x phenotype association that has p-value <= 1e-10 and FINEMAP posterior probability of causality >= 0.8
phenotype
chrom
start_pos - the start bp of the STR, 1-based
repeat_unit - canonical repeat sequence as preliminarily inferred by the paper algorithm from the reference STR sequence and the repeat period stated in the SNPSTR reference panel
white_brit_allele_frequencies - frequencies of each allele (by dosage) among the tested population
multiallelicness - amongst the population being tested for association, this is the fraction of total allelic dosage at this locus taken up by all alleles but the two most common alleles
association_p_value
direction_of_association - "+" if an increase in STR length causes an increase in phenotype, "-" if the reverse
pcausal - FINEMAP's posterior probability of causality
relation_to_gene - if this STR is transcribed, for each transcript what is the GENCODE gene type of that transcript (i.e. protein coding, lncRNA, etc.) and what is the GENCODE feature type of the region the STR is in (i.e. intron, exon, etc.)
transcribed - each transcript this STR is in (if any), including the transcript's GENCODE gene type (i.e. protein coding, lncRNA, etc.) and the transcript's transcript support level (1-5 or missing)
