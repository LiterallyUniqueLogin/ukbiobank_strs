#!/bin/bash

# run from $UKB
java -jar utilities/dxCompiler-2.11.2.jar compile \
	workflow/nexus_wdl/plot_paper_wgs_assocs.wdl \
	-verbose \
	-archive \
	-waitOnUpload \
	-separateOutputs \
	-project UKB_Test \
	-folder /imputed_strs_paper/workflow/app && \
dx run \
	/imputed_strs_paper/workflow/app/plot_paper_wgs_assocs \
	-y \
	--folder /imputed_strs_paper/wgs/20230425_paper_assocs
