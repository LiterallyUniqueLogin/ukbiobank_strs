import groovy.transform.TupleConstructor

@TupleConstructor()
class PhenotypeDescription {
    String data_field_id
    String unit
    String min_val
    String max_val
    String min_omit
    String max_omit
    List categorical_covars = []
    List previous_STR_findings = [] // format ("chr:pos", 'url')
    List previous_SNP_findings = [] // unsure of format yet
    List exciting_STR_hits = [] // format "chr:pos"
}

pheno_descs = [
    // ---------------- Continuous phenotypes ------------------------
    'height': new PhenotypeDescription(
        data_field_id: '50',
        unit : 'cm',
        previous_STR_findings : [new Tuple("3:53128363", 'https://www.nature.com/articles/s41588-019-0521-9')]
    ),
    // blood pressure - combine 4079, 94, 4080, 93,  with covars 36,37 - all category 100011
    // --- Spirometry ---
    // could study FEV and FVC z-scores (20257, 20256)
    // but those are just z-scores for those values regressed on
    // age, height and gender among lifelong non-smokers with no history of lung disease
    // we can compute that ourselves, so I don't see any reason to study these separately
    // maybe we should study lung function separately in smokers and nonsmokers
    // probably using direct question 20116 over derived field 20160
    // Spirometry: device number 3132, device id 42
    // lots of different devices used with only a few thousand per device
    // so it doesn't make sense to add that as a categorical variable
    // could add it as a random component if I was using a mixed effect model
    // 3062Forced vital capacity (FVC)
    // 3063Forced expiratory volume in 1-second (FEV1)
    // should I remove 3088,89,90? or 3159?

    // ------------------------------ Binary Phenotypes --------------------------------
    'afib_and_flutter_I48': new PhenotypeDescription(
        data_field_id: '131350',
        unit : 'binary_date_first_reported',
    ),
    'ovarian_dysfunction_E28': new PhenotypeDescription(
        data_field_id: '130736',
        unit : 'binary_date_first_reported',
    )
]

// haematology phenotypes https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100081
// 'crit' values are percentage of blood volume occupied by that cell type
// overview here https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=1453
// https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/haematology.pdf
// despite the above docs, no negative values for any of these phenotypes
haematological_phenotypes = [
    'white_blood_cell_count': new PhenotypeDescription(
        data_field_id: '30000',
        unit : '10^9 cells/L',
        exciting_STR_hits : ['17:38177325']
    ),
    'red_blood_cell_count': new PhenotypeDescription(
        data_field_id: '30010',
        unit : '10^12 cells/L',
        exciting_STR_hits : ['3:141097026', '15:76246773']
    ),
    'haemoglobin_concentration': new PhenotypeDescription(
        data_field_id: '30020',
        unit : 'g/dL',
        exciting_STR_hits : ['15:76246773'],
        previous_STR_findings : [new Tuple("16:88800373", 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889333/')]
    ),
    'haematocrit': new PhenotypeDescription(
        data_field_id: '30030',
        unit : '%',
        exciting_STR_hits : ['15:76246773'],
        previous_STR_findings : [new Tuple("16:88800373", 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889333/')]
    ),
    // here corpuscular means of red blood cells
    'mean_corpuscular_volume': new PhenotypeDescription(
        data_field_id: '30040',
        unit : '10^-15 L',
        exciting_STR_hits : ['22:32897468']
    ),
    'mean_corpuscular_haemoglobin': new PhenotypeDescription(
        data_field_id: '30050',
        unit : '10^-12 g',
        exciting_STR_hits : ['22:32897468'],
        previous_STR_findings : [new Tuple("16:88800373", 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889333/')]
    ),
    'mean_corpuscular_haemoglobin_concentration': new PhenotypeDescription(
        data_field_id: '30060',
        unit : 'g/dL',
        previous_STR_findings : [new Tuple("16:88800373", 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889333/')]
    ),
    // my (possibly incorrect) understanding is that this is the width of the image
    // the machine took relative to its receptive field and is not an innate
    // property of this cell type
    'red_blood_cell_distribution_width': new PhenotypeDescription(
        data_field_id: '30070',
        unit : '%',
        exciting_STR_hits : ['1:43429968']
    ),
    'platelet_count': new PhenotypeDescription(
        data_field_id: '30080',
        unit : '10^9 cells/L',
        exciting_STR_hits : ['11:119077000', '17:27842010']
    ),
    'platelet_crit': new PhenotypeDescription(
        data_field_id: '30090',
        unit : '%',
        exciting_STR_hits : ['1:204527033', '4:7040815']
    ),
    'mean_platelet_volume': new PhenotypeDescription(
        data_field_id: '30100',
        unit : '10^-15 L',
        exciting_STR_hits : ['17:27842010', '17:44033915']
    ),
    'platelet_distribution_width': new PhenotypeDescription(
        data_field_id: '30110',
        unit : '%',
        exciting_STR_hits : ['22:43385872']
    ),
    'lymphocyte_count': new PhenotypeDescription(
        data_field_id: '30120',
        unit : '10^9 cells/L',
    ),
    // skipping monocytes
    'neutrophil_count': new PhenotypeDescription(
        data_field_id: '30140',
        unit : '10^9 cells/L',
    ),
    'eosinophil_count': new PhenotypeDescription(
        data_field_id: '30150',
        unit : '10^9 cells/L',
        exciting_STR_hits : ['3:51119324', '10:26721459']
    ),
    // skipping basophils
    // skipping nucleated rbc measurements
    'lymphocyte_percent': new PhenotypeDescription(
        data_field_id: '30180',
        unit : '%',
    ),
    // skipping monocytes
    'neutrophil_percent': new PhenotypeDescription(
        data_field_id: '30200',
        unit : '%',
    ),
    'eosinophil_percent': new PhenotypeDescription(
        data_field_id: '30210',
        unit : '%',
        exciting_STR_hits : ['8:130623382', '19:16424879']
    ),
    // skipping basophils
    // skipping nucleated rbc measurements
    'reticulocyte_percent': new PhenotypeDescription(
        data_field_id: '30240',
        unit : '%',
        previous_STR_findings : [new Tuple("16:88800373", 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889333/')]
    ),
    'reticulocyte_count': new PhenotypeDescription(
        data_field_id: '30250',
        unit : '10^12 cells/L',
        previous_STR_findings : [new Tuple("16:88800373", 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5889333/')]
    ),
    'mean_reticulocyte_volume': new PhenotypeDescription(
        data_field_id: '30260',
        unit : '10^-15 L'
    ),
    // I think this is referring to spherical rbcs?
    'mean_sphered_cell_volume': new PhenotypeDescription(
        data_field_id: '30270',
        unit : '10^-15 L',
        exciting_STR_hits : ['17:44084292']
    ),
    'immature_reticulocyte_fraction': new PhenotypeDescription(
        data_field_id: '30280',
        unit : 'fraction'
    ),
    'high_light_scatter_reticulocyte_percentage': new PhenotypeDescription(
        data_field_id: '30290',
        unit : '%'
    ),
    'high_light_scatter_reticulocyte_count': new PhenotypeDescription(
        data_field_id: '30300',
        unit : '10^12 cells/L'
    )
]

haematological_phenotypes.each{ pheno, desc -> desc.categorical_covars.add(new Tuple(pheno + '_device_id', (desc.data_field_id.toInteger() + 3).toString())) }

pheno_descs.putAll(haematological_phenotypes)

// blood biochemistry phenotypes https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=17518
// and related QC fields https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=18518
// See related QC fields _2 for aliquot and _5, _6 for missingness reasons
// Due to the companion doc, I'm including aliquot as a covariate (see code below)
// (section 7.2)
// companion documents:
// This one for understanding the dilution/timing issues and correcting:
// https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=5636
// https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/biomarker_issues.pdf
// This one for understanding the limits of the machine's analyzing range
// https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=1227
// https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/serum_biochemistry.pdf
// units of U are probably the enzyme unit https://en.wikipedia.org/wiki/Enzyme_unit
serum_biomarkers = [
    'albumin': new PhenotypeDescription(
        data_field_id: '30600',
        unit : 'g/L',
        min_val : 15,
        max_val : 60,
        min_omit:14,
        max_omit:2,
        exciting_STR_hits : ['4:177692134']
    ),
    'alkaline_phosphatase': new PhenotypeDescription(
        data_field_id: '30610',
        unit : 'U/L',
        min_val:3,
        max_val:1500,
        min_omit:46,
        max_omit:3
    ),
    'alanine_aminotransferase': new PhenotypeDescription(
        data_field_id: '30620',
        unit : 'U/L',
        min_val:3,
        max_val:500,
        min_omit:45,
        max_omit:22,
        exciting_STR_hits : ['2:227110021']
    ),

    // 1.8k samples discarded due to too high values (>2.5)
    'apolipoprotein_a': new PhenotypeDescription(
        data_field_id: '30630',
        unit : 'g/L',
        min_val:0.4,
        max_val:2.5,
        min_omit:20,
        max_omit:1809,
        exciting_STR_hits : ['1:93858126']
    ),

    // 1.2k samples discarded (both < 0.4 and > 2)
    'apolipoprotein_b': new PhenotypeDescription(
        data_field_id: '30640',
        unit : 'g/L',
        min_val:0.4,
        max_val:2,
        max_omit:339,
        min_omit:871
    ),
    'aspartate_aminotransferase': new PhenotypeDescription(
        data_field_id: '30650',
        unit : 'U/L',
        min_val:3,
        max_val:1000,
        min_omit:11,
        max_omit:3,
        exciting_STR_hits : ['14:35662185']
    ),
    // skipping direct bilirubin 30660, 70k values are too low
    'urea': new PhenotypeDescription(
        data_field_id: '30670',
        unit : 'mmol/L',
        min_val:0.8,
        max_val:50,
        min_omit:10,
        max_omit:0
    ),
    'calcium': new PhenotypeDescription(
        data_field_id: '30680',
        unit : 'mmol/L',
        min_val:1,
        max_val:5,
        min_omit:61,
        max_omit:0,
        exciting_STR_hits : ['11:2982726', '6:34182033']
    ),
    'cholesterol': new PhenotypeDescription(
        data_field_id: '30690',
        unit : 'mmol/L',
        min_val:0.5,
        max_val:18,
        min_omit:7,
        max_omit:0
    ),
    'creatinine': new PhenotypeDescription(
        data_field_id: '30700',
        unit : 'umol/L',
        min_val:4,
        max_val:4420,
        min_omit:0,
        max_omit:0,
        exciting_STR_hits : ['2:15786517', '2:148599967']
    ),
    'c_reactive_protein': new PhenotypeDescription(
        data_field_id: '30710',
        unit : 'mg/L',
        min_val:0.08,
        max_val:80,
        min_omit:119,
        max_omit:324
    ),
    'cystatin_c': new PhenotypeDescription(
        data_field_id: '30720',
        unit : 'mg/L',
        min_val:0.1,
        max_val:'(7.7-8.99)',
        min_omit:6,
        max_omit:9
    ),
    'gamma_glutamyltransferase': new PhenotypeDescription(
        data_field_id: '30730',
        unit : 'U/L',
        min_val:5,
        max_val:1200,
        min_omit:22,
        max_omit:73
    ),
    'glucose': new PhenotypeDescription(
        data_field_id: '30740',
        unit : 'mmol/L',
        min_val:0.6,
        max_val:45,
        min_omit:7,
        max_omit:1
    ),
    'hdl_cholesterol': new PhenotypeDescription(
        data_field_id: '30760',
        unit : 'mmol/L',
        min_val:0.05,
        max_val:4.65,
        min_omit:5,
        max_omit:1,
        exciting_STR_hits : ['11:47332457']
    ),
    'igf_1': new PhenotypeDescription(
        data_field_id: '30770',
        unit : 'nmol/L',
        min_val:1.30,
        max_val:195,
        min_omit:6,
        max_omit:0
    ),
    'ldl_cholesterol_direct': new PhenotypeDescription(
        data_field_id: '30780',
        unit : 'mmol/L',
        min_val:0.26,
        max_val:10.3,
        min_omit:7,
        max_omit:0
    ),

//    //48k < 3.8 , 34k > 189
//    'lipoprotein_a': new PhenotypeDescription(
//        data_field_id: '30790',
//        unit : 'nmol/L',
//        min_val:3.8,
//        max_val:189
//    ),
//
//    // 375k < 175
//    'oestradiol': new PhenotypeDescription(
//        data_field_id: '30800',
//        unit : 'pmol/L'
//    ),
    'phosphate': new PhenotypeDescription(
        data_field_id: '30810',
        unit : 'mmol/L',
        min_val:0.32,
        max_val:6.4,
        min_omit:12,
        max_omit:0
    ),

//    // 444k < 10, 2k > 120
//    'rheumatoid_factor': new PhenotypeDescription(
//        data_field_id: '30820',
//        unit : 'IU/mL'
//    ),

    // 750 > 226-242
    'shbg': new PhenotypeDescription(
        data_field_id: '30830',
        unit : 'nmol/L',
        min_val:0.33,
        max_val:'(226-242)',
        min_omit:10,
        max_omit:745
    ),
    'total_bilirubin': new PhenotypeDescription(
        data_field_id: '30840',
        unit : 'umol/L',
        previous_STR_findings : [new Tuple("2:234668880", 'https://academic.oup.com/clinchem/article/54/5/851/5628770')],
        min_val:1,
        max_val:513,
        min_omit:18,
        max_omit:0
    ),

//    // 800 < 0.35,  41k > 55.48
//    'testosterone': new PhenotypeDescription(
//        data_field_id: '30850',
//        unit : 'nmol/L'
//    ),
    'total_protein': new PhenotypeDescription(
        data_field_id: '30860',
        unit : 'g/L',
        min_val:30,
        max_val:120,
        min_omit:16,
        max_omit:1
    ),
    'triglycerides': new PhenotypeDescription(
        data_field_id: '30870',
        unit : 'mmol/L',
        min_val:0.1,
        max_val:11.3,
        min_omit:6,
        max_omit:136
    ),
    'urate': new PhenotypeDescription(
        data_field_id: '30880',
        unit : 'umol/L',
        min_val:'unknown',
        max_val:'unknown',
        min_omit:135,
        max_omit:0
    ),

    // 2.6k < 10
    'vitamin_d': new PhenotypeDescription(
        data_field_id: '30890',
        unit : 'nmol/L',
        min_val:10,
        max_val:375,
        min_omit:2655,
        max_omit:2
    )
]

serum_biomarkers.each{ pheno, desc -> desc.categorical_covars.add(new Tuple(pheno + '_aliquot', (desc.data_field_id.toInteger() + 2).toString())) }

// this doesn't have an aliquot for some reason
serum_biomarkers['glycated_haemoglobin'] = new PhenotypeDescription(
    data_field_id: '30750',
    unit : 'mmol/mol',
    min_val:'unknown',
    max_val:'unknown',
    min_omit:190,
    max_omit:0
)

pheno_descs.putAll(serum_biomarkers)

def is_binary(phenotype) {
    return pheno_descs[phenotype].unit =~ /binary/
}

/**
phenotypes_in_use = {
    'eosinophil_count',
    'eosinophil_percent',
    'haematocrit',
    'haemoglobin_concentration',
    //'hdl_cholesterol',
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
*/

phenotypes_in_use = haematological_phenotypes.keySet().findAll( { !(it =~ /reticulocyte/) } )
phenotypes_in_use.addAll(serum_biomarkers.keySet())

phenotypes_in_use = phenotypes_in_use.sort()

