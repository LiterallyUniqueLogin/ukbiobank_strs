Project Execution
=================

This project uses Snakemake (TODO citation) to manage the workflow of the analyses
in this project. Specifically, Snakemake describes the analyses in this project
as a set of discrete tasks with named input and output files, infers the dependencies
of those tasks on each other and provides a mechanism for executing those tasks
on the cluster without reexecuting already completed worked.

The snakefile proivdes a decent look at the organization of the project:

.. details:: snakefile contents

    Location: ``$UKB/workflow/snakefile``

    .. literalinclude:: ../workflow/snakefile
        :language: python 


Configuring Snakemake
---------------------

.. details:: Configuring Snakemake

   I am using two separate files to configure Snakemake:
   ``$UKB/workflow/profile/config.yaml`` for setting the default argument values
   to the ``snakemake`` command and ``$UKB/workflow/cluster_config.yaml`` for
   setting per-rule values of the ``--cluster`` option which are passed on to
   ``sbatch``. The former uses
   `profile configuration <https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`_
   and the latter uses
   `cluster configuration <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#snakefiles-cluster-configuration>`_.

   To execute a Snakemake target with both of those configurations, use
   
   .. code:: shell

       cd $UKB 
       snakemake \
           --cluster-config workflow/cluster_config.yaml \
           --profile workflow/profile <target>

   I am not using per-rule ``resources`` in the snakefile. ``resources`` is a
   mechanism for allowing me to enforce limitations on how many jobs of specific types
   are run at any given time. Specifically, it is a map from arbitrary names
   to int values, where running ``snakemake --resources <resoure_name>=<limit>
   <target>`` ensures that jobs with no more than ``<limit>`` total of the named
   resource, as specified in the snakefile, are running at any given time. It 
   does not influence how jobs are scheduled by the cluster or what resources are
   allocated to them, and is not what I need.

   Similarly, I am not using  the ``threads`` section of each rule in the snakefile.
   That indicates how many threads the rule intends to use and Snakemake caps the
   number of jobs running at any given time so that the total used threads doesn't
   exceed ``--jobs``. Again, this doesn't tell the scheduler how many threads I
   want allocated per job.

   Each rule in the snakefile needs a corresponding entry in the cluster config
   with the parameters ``tasks-per-node``, ``time`` and ``output``. ``time``
   should be ``'HH:MM:SS'``. The quotes are *not optional*. ``output`` should be
   ``$UKB/<some_path>/output/%x_%j.out`` (``%x_%j`` will be replaced by
   ``<job_name>_<job_id>``). This corresponds to the stdout/err of the job run by the
   cluster. This is (and must be) different than the output section
   of each rule in the snakefile - that corresponds to the outputs that that rule
   is intended to create.

   After a successful run, the contents of the cluster output file will look like::

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

