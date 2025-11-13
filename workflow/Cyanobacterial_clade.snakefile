# Concatenate sequences
rule cat_sequences:
    input:
        "input/cyanobacteria.fasta",
        "input/ingroup.fasta"
    output:
        "analysis/cyanobacteria/sequences.fasta"
    params:
        m = 240,
        x = lambda w: 'X{10}'
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit seq -gim {params.m} {input} | seqkit grep -vdsp '{params.x}' | tr : _ > {output}"

# Cluster identical sequences
rule cdhit:
    input:
        "analysis/cyanobacteria/sequences.fasta"
    output:
        fasta = "analysis/cyanobacteria/sequences.fasta.cdhit",
        clstr = "analysis/cyanobacteria/sequences.fasta.cdhit.clstr"
    params:
        c = 1
    conda:
        "envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output.fasta} -c {params.c} -d 0"

# Align with mafft
rule mafft:
    input:
        "analysis/cyanobacteria/sequences.fasta.cdhit"
    output:
        "analysis/cyanobacteria/sequences.fasta.mafft"
    conda:
        "envs/mafft.yaml"
    shell:
        "mafft --auto --reorder {input} > {output}"

# Trim with trimal
rule trimal:
    input:
        "analysis/cyanobacteria/sequences.fasta.mafft"
    output:
        "analysis/cyanobacteria/sequences.fasta.trimal"
    conda:
        "envs/trimal.yaml"
    shell:
        "trimal -in {input} -out {output} -automated1"

# Run phylogeny with RAxML
rule raxml:
    input:
        "analysis/cyanobacteria/sequences.fasta.trimal"
    output:
        "analysis/cyanobacteria/RAxML_info.txt",
        "analysis/cyanobacteria/RAxML_bipartitions.txt",
        "analysis/cyanobacteria/RAxML_bestTree.txt",
        "analysis/cyanobacteria/RAxML_bipartitionsBranchLabels.txt",
        "analysis/cyanobacteria/RAxML_bootstrap.txt"
    params:
        model = "PROTGAMMAAUTO",
        seed = 123,
        bootstrap = 1000
    conda:
        "envs/raxml.yaml"
    threads:
        20
    shell:
        "raxmlHPC-PTHREADS-SSE3 -f a -p {params.seed} -x {params.seed} -# {params.bootstrap} -m {params.model} -T {threads} -s {input} -n txt -w $(dirname $(realpath {output}))"

# Plot the tree of the cyanobacterial XRs
rule plot_clade_tree:
    input:
        tree = "analysis/cyanobacteria/RAxML_bipartitions.txt",
        metadata = "input/metadata.xlsx",
        refs = "input/ingroup.fasta"
    output:
        image = "output/Cyanobacterial_XR_clade.svg",
        jtree = "output/Cyanobacterial_XR_clade.jtree"
    params:
        ingroup = [ "GR" ]
    conda:
        "envs/r-plot.yaml"
    script:
        "scripts/plot_clade_tree.R"
