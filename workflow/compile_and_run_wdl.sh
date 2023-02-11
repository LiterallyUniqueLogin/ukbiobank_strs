#!/bin/bash

# run from $UKB
java -jar workflow/dxCompiler-2.10.7.jar compile \
	workflow/research_analysis_platform_workflow.wdl \
	-verbose \
	-archive \
	-waitOnUpload \
	-separateOutputs \
	-project UKB_Test \
	-folder /imputed_strs_paper/workflow/app && \
dx run /imputed_strs_paper/workflow/app/rap_main -y
