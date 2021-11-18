nextflow.enable.dsl=2

UKB = System.getenv()['UKB']
str_imputation_run_name = 'first_pass'

process helloWorld {
    output:
        stdout emit: result
        path "hello_world.txt", emit: hello_world

    """
    echo 'Hello, world!'
    echo 'Hello, world!' > hello_world.txt
    """
}

// requires association to be done
process susie_load_gts {

    time { task.attempt == 1 ? '30m': '47h30m' }
    memory '50GB'

    input:
        val str_imputation_run_name
        val phenotype
        val chrom
        val start
        val end
    output:
        
    script:
        """
        finemapping/susie_load_gts.py $str_imputation_run_name \
        $phenotype $chrom $start $end \
        """
}

workflow {
    helloWorld()
    helloWorld.out[0].view()
}
