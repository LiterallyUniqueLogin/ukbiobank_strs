#!/bin/bash

# run from $UKB
java -jar utilities/dxCompiler-2.10.7.jar compile \
	workflow/nexus_wdl/wgs_compare_calls.wdl \
	-verbose \
	-archive \
	-waitOnUpload \
	-separateOutputs \
	-project UKB_Test \
	-folder /imputed_strs_paper/workflow/app && \
dx run \
	/imputed_strs_paper/workflow/app/compare_calls_w \
	-y \
	--folder /imputed_strs_paper/wgs/20230425_comparison
