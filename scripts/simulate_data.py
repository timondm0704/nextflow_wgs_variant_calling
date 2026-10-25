#!/usr/bin/env python3
"""Simulate a small population-genomics test dataset for the pipeline.

Creates, from one reference sequence:
  1. A set of SNPs. Some are "clinal": their alternate-allele frequency
     rises from south to north across the sampling sites. The rest have the
     same frequency everywhere.
  2. Diploid individuals sampled at several sites along a latitude gradient.
     Each individual gets a genotype (0/0, 0/1 or 1/1) at every SNP.
  3. Paired-end Illumina-style reads (2 x 100 bp) for each individual, drawn
     from both of its chromosome copies, at a chosen depth of coverage.
     One individual is deliberately sequenced at low coverage, so the
     pipeline's 12x coverage check has something to catch.
  4. A truth VCF with every individual's real genotypes, so the variant
     calls can be checked later.
  5. data/samplesheet.csv linking each individual to its FASTQ files and
     sampling metadata (latitude, longitude, depth, date).

Everything is fake but reproducible (fixed random seed).

Usage (from the project folder):
  python3 scripts/simulate_data.py
"""
import csv
import gzip
import os
import random

SEED = 42
REFERENCE = "data/reference/SHH_locus_hg38.fa"
READS_DIR = "data/reads"
TRUTH_VCF = "data/truth/truth_variants.vcf"
SAMPLESHEET = "data/samplesheet.csv"

READ_LEN = 100          # bases per read
INSERT_MEAN = 300       # DNA fragment length
INSERT_SD = 30
ERROR_RATE = 0.002      # sequencing substitution errors
TAIL_FRAC = 0.20        # reads with a low-quality 3' tail
N_SNPS = 40
N_CLINAL = 10           # SNPs whose frequency changes along the gradient
DEFAULT_COVERAGE = 15
LOW_COVERAGE = 6        # one individual below MacroMar's 12x minimum

# Four sampling sites along a south-to-north gradient (made-up coordinates)
SITES = [
    {"site": "S1", "latitude": 41.30, "longitude": -70.60, "depth_m": 12, "clinal_freq": 0.10},
    {"site": "S2", "latitude": 42.40, "longitude": -70.90, "depth_m": 8,  "clinal_freq": 0.35},
    {"site": "S3", "latitude": 43.60, "longitude": -70.20, "depth_m": 15, "clinal_freq": 0.65},
    {"site": "S4", "latitude": 44.40, "longitude": -68.20, "depth_m": 20, "clinal_freq": 0.90},
]
INDIVIDUALS_PER_SITE = 2
DATES = ["2025-06-10", "2025-06-24", "2025-07-08", "2025-07-22"]

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    break
                name = line[1:].split()[0]
            else:
                seq.append(line.strip().upper())
    return name, "".join(seq)


def add_errors(seq, rng):
    return "".join(
        rng.choice([b for b in "ACGT" if b != c]) if rng.random() < ERROR_RATE else c
        for c in seq
    )


def qualities(rng):
    q = [rng.randint(32, 40) for _ in range(READ_LEN)]
    if rng.random() < TAIL_FRAC:
        tail = rng.randint(10, 30)
        q[-tail:] = [rng.randint(2, 12) for _ in range(tail)]
    return "".join(chr(x + 33) for x in q)


def main():
    rng = random.Random(SEED)
    chrom, ref = read_fasta(REFERENCE)
    L = len(ref)

    # 1. SNPs: positions away from the ends, at least 50 bp apart
    positions = []
    while len(positions) < N_SNPS:
        p = rng.randint(400, L - 400)
        if all(abs(p - x) >= 50 for x in positions):
            positions.append(p)
    positions.sort()
    clinal = set(rng.sample(positions, N_CLINAL))
    snps = []
    for p in positions:  # p is 1-based
        ref_base = ref[p - 1]
        alt_base = rng.choice([b for b in "ACGT" if b != ref_base])
        flat_freq = round(rng.uniform(0.1, 0.5), 2)
        snps.append({"pos": p, "ref": ref_base, "alt": alt_base,
                     "clinal": p in clinal, "flat_freq": flat_freq})

    # 2. Individuals and genotypes
    individuals = []
    for si, site in enumerate(SITES):
        for k in range(1, INDIVIDUALS_PER_SITE + 1):
            sid = f"{site['site']}_ind{k}"
            hap = [[], []]  # alt allele (1) or ref (0) on each chromosome copy
            for s in snps:
                f = site["clinal_freq"] if s["clinal"] else s["flat_freq"]
                hap[0].append(1 if rng.random() < f else 0)
                hap[1].append(1 if rng.random() < f else 0)
            individuals.append({"id": sid, "site": site, "date": DATES[si], "hap": hap})
    individuals[-1]["low_cov"] = True  # last individual gets low coverage

    # 3. Reads
    os.makedirs(READS_DIR, exist_ok=True)
    for ind in individuals:
        cov = LOW_COVERAGE if ind.get("low_cov") else DEFAULT_COVERAGE
        ind["coverage"] = cov
        haplotypes = []
        for h in ind["hap"]:
            seq = list(ref)
            for s, allele in zip(snps, h):
                if allele:
                    seq[s["pos"] - 1] = s["alt"]
            haplotypes.append("".join(seq))

        n_pairs = int(cov * L / (2 * READ_LEN))
        r1_path = os.path.join(READS_DIR, f"{ind['id']}_R1.fastq.gz")
        r2_path = os.path.join(READS_DIR, f"{ind['id']}_R2.fastq.gz")
        ind["fastq_1"], ind["fastq_2"] = r1_path, r2_path
        with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
            for i in range(1, n_pairs + 1):
                hap_seq = rng.choice(haplotypes)
                ins = max(READ_LEN + 20, int(rng.gauss(INSERT_MEAN, INSERT_SD)))
                start = rng.randint(0, L - ins)
                frag = hap_seq[start:start + ins]
                if rng.random() < 0.5:
                    frag = revcomp(frag)
                read1 = add_errors(frag[:READ_LEN], rng)
                read2 = add_errors(revcomp(frag)[:READ_LEN], rng)
                name = f"{ind['id']}_pair{i}"
                r1.write(f"@{name}/1\n{read1}\n+\n{qualities(rng)}\n")
                r2.write(f"@{name}/2\n{read2}\n+\n{qualities(rng)}\n")

    # 4. Truth VCF (only SNPs where at least one individual carries the alt allele)
    os.makedirs(os.path.dirname(TRUTH_VCF), exist_ok=True)
    with open(TRUTH_VCF, "w") as v:
        v.write("##fileformat=VCFv4.2\n")
        v.write(f"##contig=<ID={chrom},length={L}>\n")
        v.write('##INFO=<ID=CLINAL,Number=0,Type=Flag,Description="Frequency changes along the latitude gradient">\n')
        v.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        v.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(i["id"] for i in individuals) + "\n")
        for j, s in enumerate(snps):
            gts = [f"{min(i['hap'][0][j], i['hap'][1][j])}/{max(i['hap'][0][j], i['hap'][1][j])}"
                   for i in individuals]
            if all(g == "0/0" for g in gts):
                continue
            info = "CLINAL" if s["clinal"] else "."
            v.write(f"{chrom}\t{s['pos']}\tsnp{j + 1}\t{s['ref']}\t{s['alt']}\t.\tPASS\t{info}\tGT\t"
                    + "\t".join(gts) + "\n")

    # 5. Sample sheet
    with open(SAMPLESHEET, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["sample_id", "site", "latitude", "longitude", "depth_m",
                    "collection_date", "fastq_1", "fastq_2"])
        for ind in individuals:
            s = ind["site"]
            w.writerow([ind["id"], s["site"], s["latitude"], s["longitude"], s["depth_m"],
                        ind["date"], ind["fastq_1"], ind["fastq_2"]])

    print(f"Reference: {chrom} ({L:,} bp)")
    print(f"SNPs simulated: {len(snps)} ({N_CLINAL} clinal)")
    for ind in individuals:
        print(f"  {ind['id']}: {ind['coverage']}x")
    print(f"Wrote {SAMPLESHEET}, {TRUTH_VCF} and reads in {READS_DIR}/")


if __name__ == "__main__":
    main()
