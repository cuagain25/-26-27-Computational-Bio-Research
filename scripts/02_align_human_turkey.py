"""Step 2: align human Arg389 WT to the 2VT4 (turkey) construct sequence.

Output: human 389 -> turkey FASTA index (1-based), identity, gap-free run, window view.
Requires Biopython (tested with 1.85).

Run from the repository root:
    python3 scripts/02_align_human_turkey.py
"""
from Bio import Align
from Bio.Align import substitution_matrices

HUMAN_PATH = "data/processed/ADRB1_human_Arg389_WT.fasta"
TURKEY_PATH = "data/raw/2VT4_turkey.fasta"
POS = 389  # 1-based human position


def read_fasta(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


human = read_fasta(HUMAN_PATH)
turkey = read_fasta(TURKEY_PATH)

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0

aln = aligner.align(human, turkey)[0]
print("score:", aln.score)

# human index (0-based) -> turkey index (0-based) for aligned residue pairs
h2t = {}
for (hs, he), (ts, te) in zip(*aln.aligned):
    for k in range(he - hs):
        h2t[hs + k] = ts + k

same = sum(human[h] == turkey[t] for h, t in h2t.items())
print("aligned pairs:", len(h2t), "| identical:", same,
      "| identity: %.1f%%" % (100 * same / len(h2t)))

i = POS - 1
t = h2t.get(i)
if t is None:
    print("human %d falls in a gap (no turkey residue)" % POS)
else:
    print("human %d:" % POS, human[i], "-> turkey FASTA index:", t + 1, "residue:", turkey[t])

    lo = hi = i
    while (lo - 1) in h2t and h2t[lo - 1] == h2t[lo] - 1:
        lo -= 1
    while (hi + 1) in h2t and h2t[hi + 1] == h2t[hi] + 1:
        hi += 1
    print("gap-free run: human %d-%d <-> turkey %d-%d"
          % (lo + 1, hi + 1, h2t[lo] + 1, h2t[hi] + 1))

    hrow = mid = trow = ""
    for h in range(i - 10, i + 10):
        hrow += human[h]
        if h in h2t:
            c = turkey[h2t[h]]
            trow += c
            mid += "|" if c == human[h] else "."
        else:
            trow += "-"
            mid += " "
    print("human  %d-%d:" % (i - 9, i + 10), hrow)
    print("               ", mid)
    print("turkey         :", trow)
