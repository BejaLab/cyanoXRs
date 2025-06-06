from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation
import csv
import re

tbl_file = snakemake.input['tbl']
fna_file = snakemake.input['fna']
output_file = str(snakemake.output)

coords_re = re.compile(r'\[([0-9]+) *- *([0-9]+)\]')
def get_location(desc):
    coord1, coord2 = coords_re.match(desc).groups()
    coord1 = int(coord1)
    coord2 = int(coord2)
    if coord1 < coord2:
        start = coord1 - 1
        stop = coord2
        strand = 1
    else:
        stop = coord1
        start = coord2 - 1
        strand = -1
    return SimpleLocation(start, stop, strand = strand)

suffix_re = re.compile(r'_[0-9]+$')
def get_record(file):
    for line in file:
        if not line.startswith('#'):
            orf = line.split()[0]
            scaffold, orf_num = orf.rsplit('_', 1)
            description = line.split(maxsplit = 18)[18]
            loc = get_location(description)
            yield orf, scaffold, loc

records = SeqIO.to_dict(SeqIO.parse(fna_file, 'fasta'))
orfs = {}

with open(tbl_file) as file, open(output_file, 'w') as out:
    for orf, scaffold, loc in get_record(file):
        if orf not in orfs:
            record = records[scaffold]
            feature = SeqFeature(location = loc, type = "CDS")
            translation = feature.translate(record, cds = False).seq.rstrip('*')
            strand = '+' if loc.strand > 0 else '-'
            out.write(f">{orf} cds:{scaffold}:{loc.start}-{loc.end}({strand}) {record.description}\n{translation}\n")
            orfs[orf] = True
