# pos_files

Position-and-frequency tables used by PosMix to spike in synthetic somatic SNVs.

## Layout

```
pos_files/
  train/          # HG001-003 pairs (manuscript training set)
    12/ 13/ 21/ 23/ 31/ 32/
      pf_<platform>_reads<N>M.txt   # exact VAF assignments used in the paper
      only*.txt                     # germline site lists used when generating VAFs
  test/           # HG004-005 pairs (manuscript test set)
    45/pf45_low.txt
    54/pf54_low.txt
  templates/      # shared low-VAF templates (pf*_low.txt) for all pairs
```

Pair IDs encode donor/background reference samples, e.g. `12` = HG001 reads mixed into HG002 background; `45` = HG004 into HG005.

## Format

Tab-delimited columns (no header), typically:

`CHROM  POS  P  REF  ALT  GT  f  b  GT2  group  AF`

`P` / `AF` are the target allele frequencies for PosMix.
