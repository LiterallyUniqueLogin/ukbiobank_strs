class PhenotypeDescription:
    def __init__(
            self,
            data_field_id,
            *,
            unit,
            min_val = None,
            max_val = None,
            min_omit = None,
            max_omit = None,
            categorical_covars = [],
            previous_STR_findings = [], # format ("chr:pos", 'url')
            previous_SNP_findings = [], # unsure of format yet
            exciting_STR_hits= []) :# format "chr:pos"
        self.data_field_id = data_field_id
        self.unit = unit
        self.min_val = min_val
        self.max_val = max_val
        self.min_omit = min_omit
        self.max_omit = max_omit
        self.categorical_covars = categorical_covars.copy()
        self.previous_STR_findings = previous_STR_findings.copy()
        self.previous_SNP_findings = previous_SNP_findings.copy()
        self.exciting_STR_hits = exciting_STR_hits.copy()

pheno_descs = {
    # ---------------- Continuous phenotypes ------------------------
    'height': PhenotypeDescription(
        '50',
        unit = 'cm',
        previous_STR_findings = [("3:53128363", 'https://www.nature.com/articles/s41588-019-0521-9')]
    ),
    # blood pressure - combine 4079, 94, 4080, 93,  with covars 36,37 - all category 100011
    # --- Spirometry ---
    # could study FEV and FVC z-scores (20257, 20256)
    # but those are just z-scores for those values regressed on
    # age, height and gender among lifelong non-smokers with no history of lung disease
    # we can compute that ourselves, so I don't see any reason to study these separately
    # maybe we should study lung function separately in smokers and nonsmokers
    # probably using direct question 20116 over derived field 20160
    # Spirometry: device number 3132, device id 42
    # lots of different devices used with only a few thousand per device
    # so it doesn't make sense to add that as a categorical variable
    # could add it as a random component if I was using a mixed effect model
    # 3062Forced vital capacity (FVC)
    # 3063Forced expiratory volume in 1-second (FEV1)
    # should I remove 3088,89,90? or 3159?

    # ------------------------------ Binary Phenotypes --------------------------------
    'afib_and_flutter_I48': PhenotypeDescription(
        '131350',
        unit = 'binary_date_first_reported',
    )
}

# haematology phenotypes https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100081
# 'crit' values are percentage of blood volume occupied by that cell type
# overview here https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=1453
# https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/haematology.pdf
# despite the above docs, no negative values for any of these phenotypes
haematological_phenotypes = {
    'white_blood_cell_count': PhenotypeDescription(
        '30000',
        unit = '10^9 cells/L',
        exciting_STR_hits = ['17:38177325']
    ),
    'red_blood_cell_count': PhenotypeDescription(
        '30010',
        unit = '10^12 cells/L',
        exciting_STR_hits = ['3:141097026', '15:76246773']
    ),
    'haemoglobin_concentration': PhenotypeDescription(
        '30020',
        unit = 'g/dL',
        exciting_STR_hits = ['15:76246773'],
        previous_STR_findings = [("16:88800373", 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889333/')]
    ),
    'haematocrit': PhenotypeDescription(
        '30030',
        unit = '%',
        exciting_STR_hits = ['15:76246773'],
        previous_STR_findings = [("16:88800373", 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889333/')]
    ),
    # here corpuscular means of red blood cells
    'mean_corpuscular_volume': PhenotypeDescription(
        '30040',
        unit = '10^-15 L',
        exciting_STR_hits = ['22:32897468']
    ),
    'mean_corpuscular_haemoglobin': PhenotypeDescription(
        '30050',
        unit = '10^-12 g',
        exciting_STR_hits = ['22:32897468'],
        previous_STR_findings = [("16:88800373", 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889333/')]
    ),
    'mean_corpuscular_haemoglobin_concentration': PhenotypeDescription(
        '30060',
        unit = 'g/dL',
        previous_STR_findings = [("16:88800373", 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889333/')]
    ),
    # my (possibly incorrect) understanding is that this is the width of the image
    # the machine took relative to its receptive field and is not an innate
    # property of this cell type
    'red_blood_cell_distribution_width': PhenotypeDescription(
        '30070',
        unit = '%',
        exciting_STR_hits = ['1:43429968']
    ),
    'platelet_count': PhenotypeDescription(
        '30080',
        unit = '10^9 cells/L',
        exciting_STR_hits = ['11:119077000', '17:27842010']
    ),
    'platelet_crit': PhenotypeDescription(
        '30090',
        unit = '%',
        exciting_STR_hits = ['1:204527033', '4:7040815']
    ),
    'mean_platelet_volume': PhenotypeDescription(
        '30100',
        unit = '10^-15 L',
        exciting_STR_hits = ['17:27842010', '17:44033915']
    ),
    'platelet_distribution_width': PhenotypeDescription(
        '30110',
        unit = '%',
        exciting_STR_hits = ['22:43385874']
    ),
    'lymphocyte_count': PhenotypeDescription(
        '30120',
        unit = '10^9 cells/L',
    ),
    # skipping monocytes
    'neutrophil_count': PhenotypeDescription(
        '30140',
        unit = '10^9 cells/L',
    ),
    'eosinophil_count': PhenotypeDescription(
        '30150',
        unit = '10^9 cells/L',
        exciting_STR_hits = ['3:51119324', '10:26721459']
    ),
    # skipping basophils
    # skipping nucleated rbc measurements
    'lymphocyte_percent': PhenotypeDescription(
        '30180',
        unit = '%',
    ),
    # skipping monocytes
    'neutrophil_percent': PhenotypeDescription(
        '30200',
        unit = '%',
    ),
    'eosinophil_percent': PhenotypeDescription(
        '30210',
        unit = '%',
        exciting_STR_hits = ['8:130623382', '19:16424879']
    ),
    # skipping basophils
    # skipping nucleated rbc measurements
    'reticulocyte_percent': PhenotypeDescription(
        '30240',
        unit = '%',
        previous_STR_findings = [("16:88800373", 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889333/')]
    ),
    'reticulocyte_count': PhenotypeDescription(
        '30250',
        unit = '10^12 cells/L',
        previous_STR_findings = [("16:88800373", 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889333/')]
    ),
    'mean_reticulocyte_volume': PhenotypeDescription(
        '30260',
        unit = '10^-15 L'
    ),
    # I think this is referring to spherical rbcs?
    'mean_sphered_cell_volume': PhenotypeDescription(
        '30270',
        unit = '10^-15 L',
        exciting_STR_hits = ['17:44084292']
    ),
    'immature_reticulocyte_fraction': PhenotypeDescription(
        '30280',
        unit = 'fraction'
    ),
    'high_light_scatter_reticulocyte_percentage': PhenotypeDescription(
        '30290',
        unit = '%'
    ),
    'high_light_scatter_reticulocyte_count': PhenotypeDescription(
        '30300',
        unit = '10^12 cells/L'
    )
}

for pheno, args in haematological_phenotypes.items():
    args.categorical_covars.append((pheno + '_device_id', str(int(args.data_field_id) + 3)))

pheno_descs.update(haematological_phenotypes)

# blood biochemistry phenotypes https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=17518
# and related QC fields https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=18518
# See related QC fields _2 for aliquot and _5, _6 for missingness reasons
# Due to the companion doc, I'm including aliquot as a covariate (see code below)
# (section 7.2)
# companion documents:
# This one for understanding the dilution/timing issues and correcting:
# https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=5636
# https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/biomarker_issues.pdf
# This one for understanding the limits of the machine's analyzing range
# https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=1227
# https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/serum_biochemistry.pdf
# units of U are probably the enzyme unit https://en.wikipedia.org/wiki/Enzyme_unit
serum_biomarkers = {
    'albumin': PhenotypeDescription(
        '30600',
        unit = 'g/L',
        min_val = 15,
        max_val = 60,
        min_omit=14,
        max_omit=2,
        exciting_STR_hits = ['4:177692134']
    ),
    'alkaline_phosphatase': PhenotypeDescription(
        '30610',
        unit = 'U/L',
        min_val=3,
        max_val=1500,
        min_omit=46,
        max_omit=3
    ),
    'alanine_aminotransferase': PhenotypeDescription(
        '30620',
        unit = 'U/L',
        min_val=3,
        max_val=500,
        min_omit=45,
        max_omit=22,
        exciting_STR_hits = ['2:227110021']
    ),

    # 1.8k samples discarded due to too high values (>2.5)
    'apolipoprotein_a': PhenotypeDescription(
        '30630',
        unit = 'g/L',
        min_val=0.4,
        max_val=2.5,
        min_omit=20,
        max_omit=1809,
        exciting_STR_hits = ['1:93858126']
    ),

    # 1.2k samples discarded (both < 0.4 and > 2)
    'apolipoprotein_b': PhenotypeDescription(
        '30640',
        unit = 'g/L',
        min_val=0.4,
        max_val=2,
        max_omit=339,
        min_omit=871
    ),
    'aspartate_aminotransferase': PhenotypeDescription(
        '30650',
        unit = 'U/L',
        min_val=3,
        max_val=1000,
        min_omit=11,
        max_omit=3,
        exciting_STR_hits = ['14:35662185']
    ),
    # skipping direct bilirubin 30660, 70k values are too low
    'urea': PhenotypeDescription(
        '30670',
        unit = 'mmol/L',
        min_val=0.8,
        max_val=50,
        min_omit=10,
        max_omit=0
    ),
    'calcium': PhenotypeDescription(
        '30680',
        unit = 'mmol/L',
        min_val=1,
        max_val=5,
        min_omit=61,
        max_omit=0,
        exciting_STR_hits = ['11:2982726', '6:34182033']
    ),
    'cholesterol': PhenotypeDescription(
        '30690',
        unit = 'mmol/L',
        min_val=0.5,
        max_val=18,
        min_omit=7,
        max_omit=0
    ),
    'creatinine': PhenotypeDescription(
        '30700',
        unit = 'umol/L',
        min_val=4,
        max_val=4420,
        min_omit=0,
        max_omit=0,
        exciting_STR_hits = ['2:15786517', '2:148599967']
    ),
    'c_reactive_protein': PhenotypeDescription(
        '30710',
        unit = 'mg/L',
        min_val=0.08,
        max_val=80,
        min_omit=119,
        max_omit=324
    ),
    'cystatin_c': PhenotypeDescription(
        '30720',
        unit = 'mg/L',
        min_val=0.1,
        max_val='(7.7-8.99)',
        min_omit=6,
        max_omit=9
    ),
    'gamma_glutamyltransferase': PhenotypeDescription(
        '30730',
        unit = 'U/L',
        min_val=5,
        max_val=1200,
        min_omit=22,
        max_omit=73
    ),
    'glucose': PhenotypeDescription(
        '30740',
        unit = 'mmol/L',
        min_val=0.6,
        max_val=45,
        min_omit=7,
        max_omit=1
    ),
    'hdl_cholesterol': PhenotypeDescription(
        '30760',
        unit = 'mmol/L',
        min_val=0.05,
        max_val=4.65,
        min_omit=5,
        max_omit=1,
        exciting_STR_hits = ['11:47732457']
    ),
    'igf_1': PhenotypeDescription(
        '30770',
        unit = 'nmol/L',
        min_val=1.30,
        max_val=195,
        min_omit=6,
        max_omit=0
    ),
    'ldl_cholesterol_direct': PhenotypeDescription(
        '30780',
        unit = 'mmol/L',
        min_val=0.26,
        max_val=10.3,
        min_omit=7,
        max_omit=0
    ),

#    #48k < 3.8 , 34k > 189
#    'lipoprotein_a': PhenotypeDescription(
#        '30790',
#        unit = 'nmol/L',
#        min_val=3.8,
#        max_val=189
#    ),
#
#    # 375k < 175
#    'oestradiol': PhenotypeDescription(
#        '30800',
#        unit = 'pmol/L'
#    ),
    'phosphate': PhenotypeDescription(
        '30810',
        unit = 'mmol/L',
        min_val=0.32,
        max_val=6.4,
        min_omit=12,
        max_omit=0
    ),

#    # 444k < 10, 2k > 120
#    'rheumatoid_factor': PhenotypeDescription(
#        '30820',
#        unit = 'IU/mL'
#    ),

    # 750 > 226-242
    'shbg': PhenotypeDescription(
        '30830',
        unit = 'nmol/L',
        min_val=0.33,
        max_val='(226-242)',
        min_omit=10,
        max_omit=745
    ),
    'total_bilirubin': PhenotypeDescription(
        '30840',
        unit = 'umol/L',
        previous_STR_findings = [("2:234668880", 'https://academic.oup.com/clinchem/article/54/5/851/5628770')],
        min_val=1,
        max_val=513,
        min_omit=18,
        max_omit=0
    ),

#    # 800 < 0.35,  41k > 55.48
#    'testosterone': PhenotypeDescription(
#        '30850',
#        unit = 'nmol/L'
#    ),
    'total_protein': PhenotypeDescription(
        '30860',
        unit = 'g/L',
        min_val=30,
        max_val=120,
        min_omit=16,
        max_omit=1
    ),
    'triglycerides': PhenotypeDescription(
        '30870',
        unit = 'mmol/L',
        min_val=0.1,
        max_val=11.3,
        min_omit=6,
        max_omit=136
    ),
    'urate': PhenotypeDescription(
        '30880',
        unit = 'umol/L',
        min_val='unknown',
        max_val='unknown',
        min_omit=135,
        max_omit=0
    ),

    # 2.6k < 10
    'vitamin_d': PhenotypeDescription(
        '30890',
        unit = 'nmol/L',
        min_val=10,
        max_val=375,
        min_omit=2655,
        max_omit=2
    )
}

for pheno, args in serum_biomarkers.items():
    args.categorical_covars.append((pheno + '_aliquot', str(int(args.data_field_id) + 2)))

# this doesn't have an aliquot for some reason
serum_biomarkers['glycated_haemoglobin'] = PhenotypeDescription(
    '30750',
    unit = 'mmol/mol',
    min_val='unknown',
    max_val='unknown',
    min_omit=190,
    max_omit=0
)

pheno_descs.update(serum_biomarkers)

def is_binary(phenotype):
    return 'binary' in pheno_descs[phenotype].unit

'''
phenotypes_in_use = {
    'eosinophil_count',
    'eosinophil_percent',
    'haematocrit',
    'haemoglobin_concentration',
    #'hdl_cholesterol',
    'ldl_cholesterol_direct',
    'lymphocyte_count',
    'lymphocyte_percent',
    'neutrophil_count',
    'neutrophil_percent',
    'platelet_count',
    'platelet_crit',
    'red_blood_cell_count',
    'total_bilirubin',
    'white_blood_cell_count',
}
'''
phenotypes_in_use = set(
    key for key in haematological_phenotypes if 'reticulocyte' not in key
)
phenotypes_in_use = phenotypes_in_use.union(serum_biomarkers)

phenotypes_in_use = sorted(phenotypes_in_use)

