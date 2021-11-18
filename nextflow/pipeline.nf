nextflow.enable.dsl=2

@Grab("com.opencsv:opencsv:5.5.2")
import com.opencsv.CSVReaderHeaderAwareBuilder
import com.opencsv.CSVParserBuilder

ukb = System.getenv('UKB')

//must be launched from $UKB

//load pheno_descs
//phenotypes_in_use
//is_binary(phenotype)
evaluate(new File("nextflow/phenotypes.groovy"))

phenotypes_in_use = ['red_blood_cell_count']

/**
Passing files around:

Processes output path qualified outputs with the same
path as the spot on the computer that I want them stored (relative paths, assumiing
$UKB as current working directory) and will
publish them with publishDir '<loc>', mode: 'link'
Either 'pattern' or 'saveAs' can be used to only publish some output files.
Note: output overwritting won't happen when using nextflow resume

Process will take paths as inputs for the files they depend on.
The scripts they execute will use those, instead of loading from hardcoded absolute paths
Input files that have already been generated and that probably won't change
will be noted in comments but not marked as inputs and will still be loaded
as hardcoded absolute paths. This should be for any file that's loaded
as a hardcoded path
*/

process testerator {
    input:
        val foo

    output:
        path "${bar}.txt"
        path "${bar}.tab"
        stdout emit: stdout
    
    when:
        println foo
        bar = "asheep"
        true

    script:
    """
    echo foo > ${bar}.txt
    echo bar > ${bar}.tab
    """
}

/**
workflow {
    for (line : file("$ukb/signals/regions/white_blood_cell_count.tab").readLines()[1..-1]) {
        println line
    }
}
*/

process summarize_finemap_output {
    label 'publish'
    time '30m'

    input:
        val phenotype
        val n_signals_no_strs
        path signal_dirs name 'dir*/*'

    output:
        path "${out}/all_STR_contribs.tab", emit: all_STR_contribs
        path "${out}/summary.txt"
        path "${out}/best_STR_contribs.tab"
        path "${out}/best_STR_ranks.tab"
        path "${out}/avg_causal_count.png"
        path "${out}/str_contrib_fractions.png"
        path "${out}/single_str_contrib_fractions.png"
        path "${out}/str_rank.png"

    script:
    out = "finemapping__finemap_results__${phenotype}__summary"
    ("mkdir ${out} && "
     "${launchDir}/finemapping/summarize_finemap_output.py " +
     "${out} " +
     "${phenotype} " +
     "${n_signals_no_strs} " +
     "${signal_dirs.inject('') { d1, d2 -> d1 + " " + d2 }} " // list out dirs
    )
}

process produce_summary_table {
    label 'publish'
    time '5h'
    memory '75GB'

    //Also depends on side_analyses/str_annotations/*
    //That may be updated if
    //* we process a new version of gencode data
    //* we change STR reference panel
    //Reads snpstr/info_field/chr*.vcf.gz
    //Will need to update if we update the STR reference panel
    input: 
        val phenotype
        path my_STR_results
        path all_STR_contribs

    output:
        path "${out_prefix}.tab"
        path "${out_prefix}_REAMDE.txt"
        stdout emit: stdout

    script:
    out_prefix = "finemapping__summary__${phenotype}_table"
    return 'foo'
    STRs = ''
    URLs = ''
    first = true
    pheno_descs[phenotype].previous_STR_findings.each{ STR, url -> 
        if (!first) {
            STRs += ' '
            URLs += ' '
        }
        first = false
        STRs += STR
        URLs += url
    }

    exciting_STRs = ''
    first = true
    pheno_descs[phenotype].exciting_STR_hits.each{ STR ->
        if (!first) {
            exciting_STRs += ' '
        }
        first = false
        exciting_STRs += STR
    }

    command = (
        "${launchDir}/finemapping/collate_strong_associations.py " +
        "${phenotype} " +
        "${out_prefix} " + 
        "'${pheno_descs[phenotype].unit}' " +
        "${my_STR_results} " +
        "${all_STR_contribs} "
    )
    if (STRs) {
        command += ' --previous-STR-findings ' + STRs
        command += ' --previous-STR-finding-URLs ' + URLs
    }
    if (exciting_STRs) {
        command += ' --cool-loci ' + exciting_STRs
    }
    println command
    command
}

workflow {

    // grab regions to finemap per phenotype. Assume regions files already exist
    pheno_to_regions = [:]
    pheno_to_n_not_strs = [:]
    phenotypes_in_use.each { phenotype -> 
        reader = new CSVReaderHeaderAwareBuilder(
            new FileReader("$ukb/signals/regions/${phenotype}.tab")
        ).withCSVParser(
            new CSVParserBuilder().withSeparator('\t' as char).build()
        ).build()
        pheno_to_regions[phenotype] = []
        pheno_to_n_not_strs[phenotype] = 0
        while (row = reader.readMap()) {
            if (row.any_strs == 'False') {
                pheno_to_n_not_strs[phenotype] += 1
            } else {
                pheno_to_regions[phenotype] << "${row.chrom}_${row.start}_${row.end}"
            }
        }
    }

    summarize_finemap_output(
        'a',
        'b',
        'c',
        Channel.fromList(phenotypes_in_use),
        Channel.fromList(phenotypes_in_use.collect { pheno_to_n_not_strs[it] } ),
        Channel.fromList(phenotypes_in_use.collect { pheno_to_regions[it] } )
    )

    // assume association results already exist
    produce_summary_table(
        Channel.fromList(phenotypes_in_use),
        Channel.fromList(phenotypes_in_use.collect { "$ukb/association/results/$it/my_str/results.tab" } ),
        summarize_finemap_output.out.all_STR_contribs
    )
    produce_summary_table.out.stdout.view()
}

