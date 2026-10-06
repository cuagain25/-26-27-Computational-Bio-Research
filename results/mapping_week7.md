# Week 7: Human ADRB1 Arg389 -> Turkey 2VT4 residue mapping

**Decision: GO** (mapping confirmed 2026-10-05)

## Mapping record

| Sequence / Structure | Residue |
| --- | --- |
| Human ADRB1 Arg389 WT (UniProt P08588, G389R edited) | Arg389 (R) |
| Turkey sequence (FASTA index, 1-313) | Arg295 (R) |
| 2VT4 Chain A (PDB residue number) | Arg355 (R) |

295 is the FASTA sequence index; 355 is the PDB residue number.

## Alignment window (human 379-399 vs turkey)

```
human  379-399: RSPDFRKAFQRLLCCARRAAR
                |||||||||.|||...|.|.|
turkey         : RSPDFRKAFKRLLAFPRKADR
```

## Numbering table (Chain A)

```
FASTA index | residue | PDB number
        290 |       R | 350
        291 |       K | 351
        292 |       A | 352
        293 |       F | 353
        294 |       K | 354
        295 |       R | 355
        296 |       L | 356
        297 |       L | 357
        298 |       A | 358
        299 |       F | missing (no coordinates)
        300 |       P | missing (no coordinates)
```

## Methods

- Sequences: human ADRB1 UniProt P08588 (477 aa); turkey beta1AR construct from RCSB FASTA for 2VT4 (313 aa, includes C-terminal His-tag).
- Arg389 WT: UniProt sequence with G389R (1 difference, verified). The UniProt reference carries Gly389.
- Alignment: Biopython 1.85 PairwiseAligner, global, BLOSUM62, gap open -10 / extend -0.5, free end gaps. 313 aligned pairs, 76.4% identity; human 389 <-> turkey 295; gap-free run human 313-407 <-> turkey 219-313 (tail end overlaps the His-tag).
- PDB numbering: Chain A CA atoms in 2VT4.pdb aligned to the construct sequence: 274 residues with coordinates (PDB 40-358), 0 mismatches; FASTA 295 <-> PDB 355.
- Scripts: `scripts/01_make_arg389_wt.py`, `02_align_human_turkey.py`, `03_fasta_to_pdb_numbering.py`.

## Limitations

- 2VT4 is a turkey beta1AR construct with stabilizing mutations, not the human WT protein.
- Chain A coordinates end at residue 358; residue 355 is 3 residues from the end of the resolved region (end of helix 8, flexible).
- Residue 358 carries the stabilizing C358A mutation.
- Sequence is well conserved around 389 (379-387 identical, 388 Q/K, 389 R/R, 390-391 identical) but diverges after ~392.

## Allele frequency context (verify on dbSNP before citing; record access date)

- rs1801253 (c.1165G>C): reference allele G encodes Gly, alt allele C encodes Arg (codon-level inference). ALFA global: G=0.29, C=0.71.
- Arg389 allele frequency 69.49% in a Southeastern European cohort: Katsarou MS et al., Front Genet 2018;9:560, doi:10.3389/fgene.2018.00560.
