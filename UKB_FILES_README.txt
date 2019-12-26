Important links:
UKBiobank data showcase website: https://biobank.ndph.ox.ac.uk/showcase/
Document describing how to download data: https://biobank.ndph.ox.ac.uk/~bbdatan/Accessing_UKB_data_v2.0.pdf
Paper describing the UK biobank microarray and imputed data: www.nature.com/articles/s41586-018-0579-z
Bycroft, C., Freeman, C., Petkova, D. et al. The UK Biobank resource with deep phenotyping and genomic data. Nature 562, 203–209 (2018) doi:10.1038/s41586-018-0579-z

The contents of this directory (each directory has its own README):
ukb_utilities - the utilities used to download all the data from the data showcase
main_dataset - contains the phenotypes we have access to through UKB
	Doesn't consist of bulk data (microarray data, exomes, imaging, etc.)
microarray - the 650k snp genotypes gathered by microarray chips
array_imputed - the 90m genotypes imputed from the microarray chips
exome - the UKB whole exome sequencing data. NOTE: this has been recalled
   and is only expected to be rereleased Spring 2020
misc_data - miscellaneous small pieces of data associated with the larger datasets
	in the above folders, but located and downloaded separately
