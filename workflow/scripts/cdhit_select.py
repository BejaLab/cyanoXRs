from Bio import SeqIO
import re

clstr_file = snakemake.input['clstr']
fasta_file = snakemake.input['fasta']
refs_file = snakemake.input['refs']

output_file = str(snakemake.output)

def parse_cluster_file(cluster_file, ref_ids):
    representatives = {}
    with open(cluster_file, 'r') as file:
        for line in file:
            if line.startswith(">Cluster"):
                current_cluster = line
            else:
                search = re.search("(\\d+)aa, >(.+)[.][.][.] (.+)", line)
                seq_id = search.group(2)
                pct = search.group(3)
                if seq_id in ref_ids or pct == "*" and current_cluster not in representatives:
                    representatives[current_cluster] = seq_id
    return representatives.values()

refs = SeqIO.to_dict(SeqIO.parse(refs_file, 'fasta'))
ref_ids = set(refs.keys())
reps = parse_cluster_file(clstr_file, ref_ids)
with open(output_file, 'w') as file:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id in reps:
            SeqIO.write(record, file, 'fasta')
