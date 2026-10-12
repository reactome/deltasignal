"""specs/021 pilot: decryptM curves -> one row per (protein, site, experiment).

decryptM (Zecha et al., Science 2023; PRIDE PXD037285, CC0) reports one fitted
dose-response curve per phosphopeptidoform. This maps each phospho group to a
residue position on the canonical UniProt sequence, so it can be joined with
Reactome's modified residues (accession, coordinate).

A peptide is placed on a REVIEWED (Swiss-Prot) canonical sequence, because
Reactome coordinates are on those: the first accession in the peptide's
protein list (isoform suffix dropped) that is reviewed and contains the
peptide exactly once. decryptM's leading protein is often an unreviewed
TrEMBL fragment (MAPK3 T202/Y204 leads with B3KR49 at 88/90), so the leading
protein alone would misplace sites. A peptide with no such accession is
skipped and counted, not guessed.

Usage: decryptm_sites.py FASTA OUT_TSV CURVES_TXT [CURVES_TXT ...]
"""
import csv
import collections
import re
import sys

MOD = re.compile(r"\(([a-z]+)\)")


def read_fasta(path, reviewed_only=True):
    seqs, acc, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if acc:
                seqs[acc] = "".join(buf)
            acc = None
            if line.startswith(">sp|") or not reviewed_only:
                acc = line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def phospho_offsets(modified):
    """(plain sequence, [0-based offsets of phosphorylated residues])."""
    s = modified.strip("_")
    plain, offs, i = [], [], 0
    while i < len(s):
        m = MOD.match(s, i)
        if m:
            if m.group(1) == "ph" and plain:
                offs.append(len(plain) - 1)
            i = m.end()
            continue
        plain.append(s[i])
        i += 1
    return "".join(plain), offs


def sites(rows, seqs, stats):
    for r in rows:
        plain, offs = phospho_offsets(r["Modified sequence"])
        if not offs:
            continue
        accs = (r.get("Leading proteins") or "").split(";") + (r.get("Proteins") or "").split(";")
        acc, start = None, -1
        for a in dict.fromkeys(x.split("-")[0] for x in accs if x):
            seq = seqs.get(a)
            if not seq:
                continue
            i = seq.find(plain)
            if i >= 0 and seq.find(plain, i + 1) < 0:
                acc, start = a, i
                break
        if acc is None:
            stats["no reviewed sequence places it uniquely"] += 1
            continue
        stats["placed"] += 1
        exp = r["Experiment"]
        for o in offs:
            yield {"uniprot": acc, "position": start + o + 1, "residue": plain[o],
                   "gene": r.get("Gene names", ""), "experiment": exp,
                   "regulation": r.get("Regulation", "") or "not",
                   "effect_size": r.get("Curve effect size", ""), "pEC50": r.get("pEC50", ""),
                   "n_sites_on_peptide": len(offs)}


def main():
    fasta, out, *curves = sys.argv[1:]
    seqs = read_fasta(fasta)
    stats = collections.Counter()
    cols = ["uniprot", "position", "residue", "gene", "experiment", "regulation",
            "effect_size", "pEC50", "n_sites_on_peptide"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for path in curves:
            csv.field_size_limit(sys.maxsize)
            for row in sites(csv.DictReader(open(path), delimiter="\t"), seqs, stats):
                w.writerow(row)
    print(dict(stats), file=sys.stderr)


if __name__ == "__main__":
    main()
