# mitotyping - vendored reference code

**Reference only. Not to be edited, and not part of this repo's pipeline.**

This is the mitoPPS reference implementation from Monzel et al. (2025), vendored into
`myc_human_validation` via the `myc_mouse` repo. It is here to be read - to check what
the original mitoPPS calculation does - not to be run, sourced, or modified.

## Chain of provenance

Upstream (see the original `SOURCE.md` in this directory, preserved unedited):

- Repo: https://github.com/annamonzel/mitotyping (branch `new_code`)
- Commit: `79b8e11`
- Cloned into `myc_mouse`: 2026-02-19
- Monzel et al. (2025), doi: 10.1101/2025.02.03.635951

Into this repo:

- Via: `myc_mouse` at `/Users/gs/G/data/MK_myc_2022/myc_mouse`
- Branch: `paper-final`
- Commit (dereferenced): `6a9c7dd513800a2a433934314a87d161ce98caa2`
- Snapshot date: 2026-08-27
- Method: `git archive` of `external/mitotyping` at that commit, so this is the
  committed state exactly. 48 files.

## What was and was not carried over

Only files **tracked** in `myc_mouse` came across. Untracked working-tree cruft was
dropped by construction: `.Rproj.user/`, `.DS_Store`, the zero-byte Google Drive `Icon`
files, and a 4.5 MB rendered `1_RNAseq_normalization_prep.html`. The tree is 1.8 MB here
against 11 MB in the source working directory. The `.Rmd` source of that HTML is
present.

## Two things that look like rule violations and are not

Flagged so a later reader does not treat them as leaks. Both were reviewed and kept
deliberately at snapshot time.

1. **`Code/Figure2/Fig2A-B_MouseData.R`** has "Mouse" in its filename and sits outside
   `data/from_myc_mouse/`. It is upstream Monzel et al. code; renaming it would corrupt
   the vendored tree. The repo rule exists to stop mouse *data* entering human analysis,
   and this is untouched third-party reference code that nothing here runs.
2. **`Data/MitoCarta/OriginalData/readme.txt`** is a 202-byte download-instruction
   pointer ("Place HumanMitoCarta3_0.xls and MouseMitoCarta3_0.xls in this folder"). It
   contains no gene data and is not a MitoCarta sheet. **No MitoCarta data file of any
   kind is present in this repo.**

The `Data/` subtree here is almost entirely such `readme.txt` pointers. The two real
data files are `Data/Fibroblasts/OriginalData/RNAseq_meta.csv.zip` (27 KB) and
`Data/Fibroblasts/ProcessedData/Lifespan_Study_data_seahorse.csv` (1.3 MB), both human
fibroblast data from the upstream lifespan study, both unused here.

## If you need to change something

Do not edit in place. Re-vendor from a newer upstream commit and update this file.
