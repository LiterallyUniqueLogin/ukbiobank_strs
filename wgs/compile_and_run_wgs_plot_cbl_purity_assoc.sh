#!/bin/bash

# run from $UKB
java -jar utilities/dxCompiler-2.10.7.jar compile \
	workflow/nexus_wdl/plot_cbl_purity_assoc.wdl \
	-verbose \
	-archive \
	-waitOnUpload \
	-separateOutputs \
	-project UKB_Test \
	-folder /imputed_strs_paper/workflow/app && \
dx run \
	/imputed_strs_paper/workflow/app/plot_cbl_purity_assoc \
	-y \
	--folder /imputed_strs_paper/wgs/20230425_comparison_2
