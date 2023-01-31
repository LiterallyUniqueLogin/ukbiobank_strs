version 1.0

import "tasks.wdl"

workflow main {

  call tasks.plot_locus as plot_SLC2A2_locus {
    phenotype_name = "total_bilirubin"
    unit = ""

  }
}
