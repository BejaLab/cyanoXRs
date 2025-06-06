# Create list of cyanobacterial scaffolds in a database
rule taxon_contigs:
    input:
        ndb = "databases/{database}/{base}.fna.ndb"
    output:
        "analysis/cyano_genomes/{database}/{base}.txt"
    params:
        db = lambda w, input: os.path.splitext(input.ndb)[0],
        taxa = 'c__Cyanobacteriia;'
    conda:
        "envs/search.yaml"
    shell:
        "blastdbcmd -db {params.db} -entry all -outfmt '%a %t' | grep -E '{params.taxa}' > {output}"

# Fetch the cyanobacterial scaffolds, translate and hmmsearch
rule hmmsearch:
    input:
        hmm = "analysis/pfam/{profile}.hmm",
        ndb = "databases/{database}/{base}.fna.ndb",
        taxon = "analysis/cyano_genomes/{database}/{base}.txt"
    output:
        "analysis/hmmsearch/{database}/{base}_{profile}.txt"
    params:
        db = lambda w, input: os.path.splitext(input.ndb)[0],
        min = 50*3,
        max = 10000*3
    conda:
        "envs/search.yaml"
    threads:
        4
    shell:
        """
        cut -f1 -d' ' {input.taxon} | blastdbcmd -entry_batch - -db {params.db} | \
            sed 's/:/%3A/g' | getorf -find 0 -filter -minsize {params.min} -maxsize {params.max} | sed 's/%3A/:/g' | \
            hmmsearch --cpu {threads} --cut_ga --tblout {output} -o /dev/null {input.hmm} -
        """

# Assign colors to taxa (for consistency across databases)
rule taxa_colors:
    input:
        expand("analysis/cyano_genomes/{database}/{base}.txt", zip, database = databases, base = bases)
    output:
        "analysis/cyano_genomes/{rank}.csv"
    conda:
        "envs/r-plot.yaml"
    script:
        "scripts/taxa_colors.R"

# Extract the scaffolds matched with hmmsearch
rule extract_contigs:
    input:
        tbl = "analysis/hmmsearch/{database}/{base}_{profile}.txt",
        ndb = "databases/{database}/{base}.fna.ndb"
    output:
        "analysis/hmmsearch/{database}/{base}_{profile}.fna"
    params:
        db = lambda w, input: os.path.splitext(input.ndb)[0]
    conda:
        "envs/search.yaml"
    shell:
        "grep -v '^#' {input.tbl} | cut -f1 -d' ' | sed -E 's/_[0-9]+$//' | sort -u | blastdbcmd -db {params.db} -entry_batch - > {output}"

# Extract the matching ORF translations based on the hmmsearch matched
rule extract_contig_proteins:
    input:
        fna = "analysis/hmmsearch/{database}/{base}_{profile}.fna",
        tbl = "analysis/hmmsearch/{database}/{base}_{profile}.txt"
    output:
        "analysis/hmmsearch/{database}/{base}_{profile}.faa"
    conda:
        "envs/python.yaml"
    script:
        "scripts/extract_contig_proteins.py"

# Blast the extracted proteins against the rhodopsin list
rule blast_contig_proteins:
    input:
        query = "analysis/hmmsearch/{database}/{base}_{profile}.faa",
        db = "analysis/rhodopsin_list/Rhodopsin_list.faa",
        pdb = "analysis/rhodopsin_list/Rhodopsin_list.faa.pdb"
    output:
        "analysis/hmmsearch/{database}/{base}_{profile}_rhodopsins.blast"
    params:
        outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle',
        max = 10
    threads:
        10
    shell:
        "blastp -num_threads {threads} -max_target_seqs {params.max} -query {input.query} -db {input.db} -out {output} -outfmt '{params.outfmt}'"

# Plot the distribution of the rhodopsins in cyanobacteria
rule plot_families:
    input:
        blast = "analysis/hmmsearch/{{database}}/{{base}}_{profile}_rhodopsins.blast".format(profile = profile),
        scaffolds = "analysis/cyano_genomes/{{database}}/{{base}}.txt".format(profile = profile),
        hmmsearch = "analysis/hmmsearch/{{database}}/{{base}}_{profile}.txt".format(profile = profile),
        taxa = "analysis/cyano_genomes/order.csv"
    output:
        plot = "output/Rhodopsins_in_cyanobacteria_{database}-{base}.svg",
        counts = "output/Rhodopsins_in_cyanobacteria_{database}-{base}.csv"
    params:
        min_id = 50,
        evalue = 1e-10
    conda:
        "envs/r-upset.yaml"
    script:
        "scripts/plot_families.R"
