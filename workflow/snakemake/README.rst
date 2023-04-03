Project Execution
=================

This project uses Snakemake to manage the workflow of the analyses
in this project. Specifically, Snakemake describes the analyses in this project
as a set of discrete tasks with named input and output files, infers the dependencies
of those tasks on each other and provides a mechanism for executing those tasks
on the cluster without reexecuting already completed worked.

The snakefile proivdes a decent look at the organization of all the analyses in the project

Configuring Snakemake
---------------------

I am using the config.yaml for setting the default argument values
to the ``snakemake`` command for running on the San Diego Supercomputer Center's Expanse cluster.
See that file for more details. This uses the snakemake profile functionality described in
`profile configuration <https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`_ That file
also uses the status-sacct.sh script for monitoring SLURM jobs.

To execute a Snakemake target with that configuration, use
   
.. code:: shell

    cd $UKB 
    snakemake \
        --profile workflow -rp <target>

Note that snakemake must be run from the root directory of the project

Each rule in the snakefile needs a resources section 
with the ``time`` parameter. The rules may override the default values for the parameters
``threads`` (1) and ``mem_gb`` (2). ``time``
should be ``'HH:MM:SS'``. The quotes are *not optional*. 

After a successful run, the contents of the output file will look like::

    Building DAG of jobs...
    Using shell: /usr/bin/bash
    Provided cores: 1 (use --cores to define parallelism)
    Rules claiming more threads will be scaled down.
    Job counts:
           count   jobs
           1       test
           1
    Select jobs to execute...
       
    [Fri Feb 26 17:36:09 2021]
    rule test:
        output: out.txt
        jobid: 0
           
     <stdout of the script executed by this rule>
    [Sat Feb 27 09:41:32 2021]
    Finished job 0.
    1 of 1 steps (100%) done

If using the shell: option in a rule in the snakefile, the current working
directory of the shell is set to the current working directory where snakemake
was invoked (not where the snakefile is located). Similarly (per
`here <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts>`_)
input and output paths are also relative to the invocation directory and not
the snakefile directory. Note: the target of a script: option is relative to
the snakefile and not the current working directory of the snakemake command!
For this reason, I'm using shell: and not script: .

Phenotype Configuration
-----------------------
The phenotypes.py file is used to define relevant aspects of the phenotypes that
were analyzed in this project or other side analyses. This is referred to directly
by the snakefile

Project Environment
===================
`envs/` contains yaml files describing the conda environments used in this project.
software_versions_paper_submission.txt contains a manually curated list of all pieces
of software used for the first paper submission. This was created in part with the
the conda_pkg_versions.sh script

