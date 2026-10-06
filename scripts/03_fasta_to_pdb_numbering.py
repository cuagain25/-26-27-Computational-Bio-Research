"""Step 3: convert a turkey FASTA index to the PDB residue number (2VT4, Chain A).

Residues that have coordinates (CA atoms) in Chain A are aligned to the construct
sequence, so gaps from unresolved regions are handled correctly.
Requires Biopython (tested with 1.85).

Run from the repository root:
    python3 scripts/03_fasta_to_pdb_numbering.py
"""
from Bio import Align
from Bio.SeqUtils import seq1

TURKEY_PATH = "data/raw/2VT4_turkey.fasta"
PDB_PATH = "data/raw/2VT4.pdb"
CHAIN = "A"
FASTA_INDEX = 295  # 1-based index in the turkey construct sequence (from step 2)

turkey = "".join(l.strip() for l in open(TURKEY_PATH) if not l.startswith(">"))

# Chain residues that actually have coordinates (CA atoms)
obs, seen = [], set()
for l in open(PDB_PATH):
    if l.startswith("ATOM") and l[21] == CHAIN and l[12:16].strip() == "CA":
        num = l[22:26].strip() + l[26].strip()  # residue number + insertion code
        if num in seen:  # skip alternate locations
            continue
        seen.add(num)
        obs.append((num, l[17:20]))
obs_seq = "".join(seq1(r) for _, r in obs)
print("chain %s residues with coordinates:" % CHAIN, len(obs_seq),
      "| PDB numbers", obs[0][0], "to", obs[-1][0])

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -1
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0

aln = aligner.align(turkey, obs_seq)[0]
t2o = {}
for (ts, te), (os_, oe) in zip(*aln.aligned):
    for k in range(te - ts):
        t2o[ts + k] = os_ + k

mism = sum(turkey[t] != obs_seq[o] for t, o in t2o.items())
print("matched:", len(t2o), "| mismatches:", mism)

print("\nFASTA index | residue | PDB number")
center = FASTA_INDEX - 1
for t in range(center - 6, center + 5):
    if t in t2o:
        print("%11d | %7s | %s" % (t + 1, turkey[t], obs[t2o[t]][0]))
    else:
        print("%11d | %7s | missing (no coordinates)" % (t + 1, turkey[t]))
