"""Step 1: build the Arg389 WT sequence from the UniProt reference.

UniProt P08588 carries Gly at position 389 (the minor-allele form).
This script substitutes G389R and writes a new FASTA; the original is untouched.

Run from the repository root:
    python3 scripts/01_make_arg389_wt.py
"""
IN_PATH = "data/raw/ADRB1_human_P08588.fasta"
OUT_PATH = "data/processed/ADRB1_human_Arg389_WT.fasta"
POS = 389  # 1-based position in the UniProt sequence

seq = "".join(l.strip() for l in open(IN_PATH) if not l.startswith(">"))
assert seq[POS - 1] == "G", "expected Gly at position %d" % POS  # python is 0-based

wt = seq[:POS - 1] + "R" + seq[POS:]
with open(OUT_PATH, "w") as f:
    f.write(">ADRB1_human_P08588_G389R (Arg389 WT, edited from UniProt)\n")
    for i in range(0, len(wt), 60):
        f.write(wt[i:i + 60] + "\n")

print("length:", len(wt))
print("residue %d:" % POS, wt[POS - 1])
print("residues %d-%d:" % (POS - 5, POS + 5), wt[POS - 6:POS + 5])
print("differences from original:", sum(a != b for a, b in zip(seq, wt)))
