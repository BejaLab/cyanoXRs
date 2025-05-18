import yaml
from os.path import basename

rule all_pgap:
    input:
        expand("analysis/pgap/{genome}/annot.gbk", genome = genomes)

rule copy_seq:
    input:
        "genomes/{genome}.fna"
    output:
        "analysis/pgap/{genome}.fna"
    shell:
        "cp {input} {output}"

rule copy_submol:
    input:
        "workflow/resources/cyanobacterium.yaml"
    output:
        "analysis/pgap/{genome}_submol.yaml"
    shell:
        "cp {input} {output}"

rule make_yaml:
    input:
        fasta  = "analysis/pgap/{genome}.fna",
        submol = "analysis/pgap/{genome}_submol.yaml"
    output:
        "analysis/pgap/{genome}.yaml"
    wildcard_constraints:
        genome = '.*(?<!_submol)'
    run:
        sets = dict(
            fasta  = { "class": "File", "location": basename(input.fasta)  },
            submol = { "class": "File", "location": basename(input.submol) }
        )
        with open(output[0], 'w') as fd:
            yaml.dump(sets, fd)

rule pgap:
    input:
        yaml = "analysis/pgap/{genome}.yaml",
        submol = "analysis/pgap/{genome}_submol.yaml"
    output:
        gbk = "analysis/pgap/{genome}/annot.gbk"
    params:
        dirname = directory("analysis/pgap/{genome}/")
    threads:
        8
    shell:
        """
        rm -fr {params.dirname}
        pgap.py --no-internet --no-self-update --debug --ignore-all-errors --report-usage-true --cpu {threads} -o {params.dirname} {input.yaml}
        """
