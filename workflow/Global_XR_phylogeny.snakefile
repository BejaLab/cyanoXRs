# Extract a single ingroup protein
rule extract_protein:
    input:
        "input/ingroup.fasta"
    output:
        "analysis/refs/proteins/{protein}.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit grep -p {wildcards.protein} {input} | seqkit seq -i -o {output}"

# tblastn the databases with SrXR
rule tblastn:
    input:
        query = "analysis/refs/proteins/XR.faa",
        ndb = "databases/{database}/{base}.fna.ndb"
    output:
        "analysis/tblastn/{database}/{base}.txt"
    params:
        db = lambda w, input: os.path.splitext(input.ndb)[0]
    conda:
        "envs/search.yaml"
    threads:
        20
    shell:
        "tblastn -num_threads {threads} -db {params.db} -query {input.query} -outfmt 6 -out {output} -max_target_seqs 1000000000"

# Extract matching scaffolds
rule extract_scaffolds:
    input:
        tblastn = "analysis/tblastn/{database}/{base}.txt",
        ndb = "databases/{database}/{base}.fna.ndb"
    output:
        "analysis/tblastn/{database}/{base}.fna"
    params:
        db = lambda w, input: os.path.splitext(input.ndb)[0],
        max_e = 1e-10,
        min_len = 200
    conda:
        "envs/search.yaml"
    shell:
        "awk -v E={params.max_e} -v L={params.min_len} '$11<E && $4>=L' {input.tblastn} | cut -f2 | sort -u | blastdbcmd -db {params.db} -entry_batch - > {output}"

# Collect scaffold headers with genome and taxonomy information
rule scaffold_headers:
    input:
        expand("analysis/tblastn/{db}/{base}.fna", zip, db = databases, base = bases)
    output:
        "analysis/tblastn/scaffold_headers.txt"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit seq -n {input} | awk '!_[$1];{{_[$1]=1}}' > {output}"

# Predict genes in the scaffolds
rule prodigal_scaffolds:
    input:
        "analysis/tblastn/{database}/{base}.fna"
    output:
        gff = "analysis/tblastn/{database}/{base}.gff",
        faa = "analysis/tblastn/{database}/{base}.faa"
    conda:
        "envs/prodigal.yaml"
    shell:
        "prodigal -i {input} -f gff -o {output.gff} -a {output.faa} -p meta"

# Blast the predicted proteins against the list of prokaryotic rhodopsin families 
rule blast_rhods:
    input:
        query = "analysis/tblastn/{database}/{base}.faa",
        db = "analysis/rhodopsin_list/Rhodopsin_list.faa",
        pdb = "analysis/rhodopsin_list/Rhodopsin_list.faa.pdb"
    output:
        "analysis/tblastn/{database}/{base}_rhodopsin_list.txt"
    params:
        outfmt = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle',
        max_target_seqs = 5,
        evalue = 1e-10
    conda:
        "envs/search.yaml"
    threads:
        20
    shell:
        "blastp -query {input.query} -db {input.db} -out {output} -outfmt '6 {params.outfmt}' -max_target_seqs {params.max_target_seqs} -evalue {params.evalue} -num_threads {threads}"

# Extract the XRs based on best blast matches
rule blast_rhods_extract_xrs:
    input:
        blast = "analysis/tblastn/{database}/{base}_rhodopsin_list.txt",
        faa = "analysis/tblastn/{database}/{base}.faa"
    output:
        "analysis/tblastn/{database}/{base}_xrs.faa"
    params:
        clade = '/P:XR/',
        id = 40,
        evalue = 1e-30
    conda:
        "envs/kits.yaml"
    shell:
        "awk -vc={params.clade} -F\\\\t '!_[$1] && $13~c && $11<{params.evalue} && $3>={params.id} {{print$1}}; {{_[$1]=1}}' {input.blast} | seqkit grep -f- {input.faa} -o {output}"

# Combine the ingroups, outgroups and the collected XRs
rule cat_proteins:
    input:
        "input/ingroup.fasta",
        "input/outgroups.fasta",
        expand("analysis/tblastn/{database}/{base}_xrs.faa", zip, database = databases, base = bases)
    output:
        "analysis/proteins/big.faa"
    conda:
        "envs/kits.yaml"
    shell:
        "seqkit rmdup -o {output} {input}"

# Cluster the proteins
rule cdhit_xrs:
    input:
        "analysis/proteins/big.faa"
    output:
        fasta = "analysis/proteins/big_cdhit.faa",
        clstr = "analysis/proteins/big_cdhit.faa.clstr"
    params:
        c = 0.95
    conda:
        "envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output.fasta} -c {params.c} -d 0"

# Select the cluster representatives to guarantee that ingroup refs are present
rule cdhit_select_proteins:
    input:
        refs = "input/ingroup.fasta",
        fasta = "analysis/proteins/big.faa",
        clstr = "analysis/proteins/big_cdhit.faa.clstr"
    output:
        "analysis/proteins/big_cdhit_selected.faa"
    conda:
        "envs/python.yaml"
    script:
        "scripts/cdhit_select.py"

# Align
rule mafft_xrs:
    input:
        "analysis/proteins/big_cdhit_selected.faa"
    output:
        "analysis/proteins/big_cdhit.mafft"
    conda:
        "envs/mafft.yaml"
    threads:
        20
    shell:
        "mafft --thread {threads} --reorder --auto {input} > {output}"

# Trim
rule trimal_xrs:
    input:
        "analysis/proteins/big_cdhit.mafft"
    output:
        "analysis/proteins/big_cdhit.trimal"
    params:
        gt = 0.9
    conda:
        "envs/trimal.yaml"
    shell:
        "trimal -in {input} -out {output} -gt {params.gt} -keepheader"

# Do iqtree2 phylogeny
rule gene_iqtree:
    input:
        "analysis/proteins/big_cdhit.trimal"
    output:
        expand("analysis/proteins/big_cdhit.{ext}", ext = [ "treefile", "model.gz", "iqtree", "ckp.gz" ])
    params:
        seed = 123,
        B = 1000,
        prefix = "analysis/proteins/big_cdhit"
    threads:
        4
    conda:
        "envs/iqtree.yaml"
    shell:
        "iqtree2 -s {input} --prefix {params.prefix} -redo --alrt 1000 -B {params.B} --seed {params.seed} -T {threads} -msub nuclear"

# Plot the tree
rule plot_global_tree:
    input:
        scaffolds = "analysis/tblastn/scaffold_headers.txt",
        tree = "analysis/proteins/big_cdhit.treefile",
        synonyms = "workflow/resources/phyla_synonyms.txt",
        ref_taxonomy = "input/ingroup_taxonomy.txt",
        outgroups = "input/outgroups.fasta"
    output:
        image = "output/All_prokaryotic_XRs.svg",
        jtree = "output/All_prokaryotic_XRs.jtree"
    conda:
        "envs/r-plot.yaml"
    script:
        "scripts/plot_global_tree.R"
