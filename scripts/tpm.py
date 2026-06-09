#how to run this script:
#python tpm.py --counts counts.tsv --gtf genes.gtf --out tpm.tsv
#try: --counts /home/cvillazon/results/mRNA_SNAKEPipes_out/featureCounts/counts_edit.tsv
#gtf: --gtf /home/seisbio/data/snakePipes/GRCm39/annotation/genes.gtf
#out: --out /home/cvillazon/results/tpm_counts.csv


import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser(
    description="Convert raw counts to TPM using gene lengths from a GTF file.")
parser.add_argument("--counts", help="Counts file")
parser.add_argument("--gtf", help="GTF file")
parser.add_argument("--out", help="Output file")

args = parser.parse_args()

counts_file = args.counts
gtf_file = args.gtf
out_file = args.out

# -------------------------
# Load counts
# -------------------------
counts = pd.read_csv(counts_file, sep="\t", index_col=0)

# -------------------------
# Extract gene lengths from GTF
# -------------------------
gene_lengths = {}

with open(gtf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if fields[2] != "exon":
            continue

        start = int(fields[3])
        end = int(fields[4])
        length = end - start + 1

        attrs = fields[8]
        gene_id = re.search('gene_id "([^"]+)"', attrs).group(1)

        gene_lengths.setdefault(gene_id, 0)
        gene_lengths[gene_id] += length

gene_lengths = pd.Series(gene_lengths)

# -------------------------
# Align genes
# -------------------------
common = counts.index.intersection(gene_lengths.index)

counts = counts.loc[common]
gene_lengths = gene_lengths.loc[common]

# -------------------------
# Compute TPM
# -------------------------
rpk = counts.div(gene_lengths / 1000, axis=0)
tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6

# -------------------------
# Save
# -------------------------
tpm.to_csv(out_file)