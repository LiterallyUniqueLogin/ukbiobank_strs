Environment Setup
=================

.. details:: Environment setup

    Conda packages
    --------------

    Environment file

    .. literalinclude:: ../ukb.env.yml

    Package versions

    .. 
        .. command-output:: list=$(conda list) ; echo "$list" | head -n 3 | tail -n 1 ; for pkg in $(grep -m1 -A50000 dependencies ../ukb.env.yml | tail -n +2 | awk '{print $2}') ; do echo "$list" | grep $pkg ; done
            :shell:

    Other software dependencies
    ---------------------------

    .. command-output:: ../utilities/plink2 --version

    .. command-output:: ../utilities/hipstr/hipstr/HipSTR --version
