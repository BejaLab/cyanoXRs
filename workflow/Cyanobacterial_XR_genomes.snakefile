# Taxonomic classification of the genomes
rule gtdbtk:
    input:
        "genomes"
    output:
        directory("analysis/gtdbtk")
    params:
        data_path = "/data/gtdbtk/release214/214.1/auxillary_files/gtdbtk_r214_data/release214"
    conda:
        "envs/gtdbtk.yaml"
    shell:
        "GTDBTK_DATA_PATH='{params.data_path}' gtdbtk classify_wf --genome_dir {input} --out_dir {output} --skip_ani_screen"

# Copy a genome fasta to analysis/
rule copy_genome:
    input:
        "genomes/{genome}.fna"
    output:
        "analysis/genomes/{genome}.fna"
    shell:
        "cp {input} {output}"

# Make blast database for a genome
rule makeblast_nucl:
    input:
        "analysis/genomes/{genome}.fna"
    output:
        "analysis/genomes/{genome}.fna.ndb"
    conda:
        "envs/search.yaml"
    shell:
        "makeblastdb -in {input} -dbtype nucl"

# Predict genes in a XR-containing cyanobacterial genome
rule prodigal:
    input:
        "analysis/genomes/{genome}.fna"
    output:
        gff = "analysis/genomes/{genome}.gff",
        faa = "analysis/genomes/{genome}.faa"
    conda:
        "envs/prodigal.yaml"
    shell:
        "prodigal -i {input} -f gff -o {output.gff} -a {output.faa}"

# Download genes for blast searches
rule dload_gene_reps:
    output:
        "analysis/genes/{gene}.fasta"
    params:
        query = lambda w: '+AND+'.join(genes[w.gene])
    shell:
        "wget -qO {output} 'https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query={params.query}'"

# blast selected genes against genome proteins
rule blastp:
    input:
        query = "analysis/genes/{gene}.fasta",
        db = "analysis/genomes/{genome}.faa",
        pdb = "analysis/genomes/{genome}.faa.pdb"
    output:
        "analysis/blastp/{gene}/{genome}.txt"
    params:
        fmt = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'
    conda:
        "envs/search.yaml"
    shell:
        "blastp -db {input.db} -query {input.query} -outfmt '6 {params.fmt}' -out {output}"

# Extract matches proteins
rule blastp_fasta:
    input:
        matches = expand("analysis/blastp/{{gene}}/{genome}.txt", genome = genomes),
        fasta = expand("analysis/genomes/{genome}.faa", genome = genomes)
    output:
        "output/{gene}_homologs.fasta"
    params:
        e = 1e-10
    conda:
        "envs/kits.yaml"
    shell:
        "awk '$11<{params.e}' {input.matches} | cut -f2 | seqkit grep -f- {input.fasta} | seqkit rmdup -o {output}"
