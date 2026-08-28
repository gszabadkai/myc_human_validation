---
date: 2026-08-28
status: G1 count criterion discharged; correlation criterion OPEN
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 2, 6, 7.1, 8, 10)
decides:
  - D2 - RESOLVED, M-a Felsher stays primary
  - D7 - NEW, PROPOSED, awaiting confirmation before script 09
next-action: resolve D4 (Lee et al. 2017), then G2
---

# G1 result and consequent decisions

Gate G1, the MYC signature overlap audit (plan section 6). Run 2026-08-28 via
`scripts/04_snapshot_human_genesets.R`. Results in `results/g1_overlap_audit.rds`
and `outputs/tables/g1_overlap_audit.csv`.

---

## 1. Verdict

**The count criterion passes for all three estimators. The mitochondrial
circularity G1 exists to detect is absent.** A MYC-to-OXPHOS result in this
cohort will not be true by construction.

**G1 is NOT fully discharged.** Plan section 6 fails the signature on either of
two criteria: falling below ~50 genes, or losing its correlation structure. Only
the first is answered here. The second needs expression data from script 01. Do
not cite G1 as closed until that runs.

## 2. Strip result (full 1,136-gene MitoCarta inventory)

| Estimator | Before | Removed | Remaining | Threshold 50 |
|---|---|---|---|---|
| FELSHER | 67 | 6 | **61** | PASSES |
| COLLECTRI_MYC_ALL | 886 | 75 | 811 | PASSES |
| COLLECTRI_MYC_STIM | 810 | 71 | 739 | PASSES |

## 3. Mitochondrial overlap - clean, and cleaner than expected

| Estimator | MitoCarta ALL | OXPHOS umbrella (169) | OXPHOS subunits (102) | mtDNA (13) |
|---|---|---|---|---|
| FELSHER | 6 (9.0%) | **0** | **0** | 0 |
| COLLECTRI_MYC_ALL | 75 (8.5%) | 6 | 3 | 1 |
| COLLECTRI_MYC_STIM | 71 (8.8%) | 6 | 3 | 1 |

The Felsher overlap with the plan's primary OXPHOS measure is zero. Not small -
zero. The six mitochondrial genes sit nowhere near the electron transport chain:

```
APEX1    mtDNA maintenance | mtDNA repair
KARS1    mt-tRNA synthetases
MARS2    mt-tRNA synthetases
MTHFD2   Folate and 1-C metabolism
NOA1     mtRNA granules | Mitochondrial ribosome assembly
PUS1     mt-tRNA modifications
```

## 4. Proliferation entanglement - real, and NOT fixed by stripping

Overlap with the Hallmark proliferation sets, against chance (hypergeometric,
background 19,247 - the shared MitoCarta Sheet-3 vocabulary):

| Estimator | Version | Set | n | Expected | Fold | p |
|---|---|---|---|---|---|---|
| FELSHER | raw | E2F_TARGETS | 6 | 0.7 | 8.6x | 6.9e-5 |
| FELSHER | raw | G2M_CHECKPOINT | 5 | 0.7 | 7.2x | 6.6e-4 |
| FELSHER | stripped | E2F_TARGETS | 5 | 0.6 | 7.9x | 4.3e-4 |
| FELSHER | stripped | G2M_CHECKPOINT | 5 | 0.6 | 7.9x | 4.3e-4 |
| COLLECTRI_MYC_ALL | raw | E2F_TARGETS | 58 | 9.2 | 6.3x | 1.2e-30 |
| COLLECTRI_MYC_ALL | raw | G2M_CHECKPOINT | 49 | 9.2 | 5.3x | 1.7e-22 |

Stripping MitoCarta genes does not reduce it - the fold enrichment is flat or
slightly higher afterwards, because the denominator shrinks. This is a separate
problem from mitochondrial circularity and the strip was never going to fix it.

**The enrichment itself is expected biology, not a defect.** MYC drives
proliferation. The problem is narrower and structural, and it is stated in D7.

---

## 5. D2 - RESOLVED. M-a (Felsher, MitoCarta-stripped) stays primary

The plan specified M-a primary with mandatory M-b/M-c concordance, and asked for
confirmation or a flip. **Confirmed, not flipped.**

Reasoning: the stripped signature clears the threshold with room (61 vs 50), and
on the dimension G1 actually tests it is the *cleaner* of the two estimators -
zero OXPHOS-subunit overlap against CollecTRI's three. The cross-species
continuity argument with the mouse arm therefore costs nothing in circularity.

M-b (CollecTRI) and M-c (8q24 GISTIC) remain mandatory concordance checks. The
rule in plan section 7.1 is unchanged: if M-a and M-b disagree, report both and
treat the claim as unsupported.

**Which CollecTRI variant is M-b:** use COLLECTRI_MYC_ALL (886 targets, 811
stripped). The signed `mor` handles direction at scoring time, so discarding the
inhibitory edges up front would throw away information the ULM/VIPER model uses.
Note 92 MYC edges are flagged as both stimulatory and inhibitory; they are
retained in ALL and also counted as stimulatory in STIM.

## 6. D7 - NEW. The proliferation covariate overlaps the MYC exposure

**Status: PROPOSED. Must be confirmed before script 09 is built, and before any
model is fitted.**

Plan section 8 specifies a proliferation covariate (Hallmark E2F_TARGETS +
G2M_CHECKPOINT GSVA) and plan section 10 puts it in the Block C model:

```
PRIME ~ MYC * OXPHOS + purity + leukocyte_fraction + proliferation
        + PAM50 + TP53_status + plate
```

Gene-level sharing between the exposure and that covariate:

| Estimator (stripped) | In proliferation covariate | Share |
|---|---|---|
| FELSHER | 9 of 61 | **14.8%** |
| COLLECTRI_MYC_ALL | 80 of 811 | 9.9% |
| COLLECTRI_MYC_STIM | 78 of 739 | 10.6% |

This is not confounding. It is **overlapping measurement**: roughly one in seven
genes constituting the MYC score also constitutes the covariate meant to adjust
it. Adjusting for proliferation therefore partially adjusts away MYC itself,
biasing the `MYC:OXPHOS` coefficient toward the null. No sample size fixes this;
n ~1000 makes it more precise, not less biased.

**Proposed handling, all three:**

1. **Primary model unchanged.** Keep the covariate as pre-specified. Changing the
   primary model on the basis of a gate result is exactly the drift this arm's
   pre-registration exists to prevent.
2. **Pre-specified sensitivity:** refit without the proliferation term and report
   both coefficients side by side. Pre-specified now, before any outcome data is
   seen, so it is a sensitivity analysis and not a post-hoc rescue.
3. **Disjoint covariate:** construct a second proliferation score from E2F/G2M
   genes *not* in the MYC estimator, and report the interaction under it. This is
   the clean comparison and it costs one extra GSVA run.

If the interaction survives all three, the result is robust to the construction.
If it appears only without the covariate, that is informative and must be
reported as such, not selected.

## 7. Caveat to carry into the specificity battery

`MTHFD2` is in MitoCarta's **Folate and 1-C metabolism** pathway, and one-carbon
/ glycine cleavage is one of the pre-specified *pathway negatives* (plan section
2). MTHFD2 is also a canonical MYC target.

The strip removes it from the estimator, so the estimator is clean. The caveat is
interpretive: if the one-carbon negative later shows signal, the honest reading
is "MYC drives one-carbon metabolism" - long established - rather than "the
OXPHOS coupling is non-specific". The same caveat applies more weakly to the
mitoribosome negative via `NOA1`, and to any mt-tRNA-synthetase-adjacent set via
`KARS1` / `MARS2`. **State this in Methods.**

## 8. Symbol harmonisation - it mattered

Five of the 67 Felsher symbols are deprecated and match nothing in a
19,247-gene background. One changes the answer: `KARS` resolves to `KARS1`, which
IS in MitoCarta. A naive symbol join reports the Felsher-MitoCarta overlap as 5
when it is 6.

```
CD3EAP -> POLR1G     CIRH1A -> UTP4      VARS -> VARS1
METTL13 -> EEF1AKNMT  KARS -> KARS1  (mitochondrial)
```

Eight symbols across all sets were unresolvable and are kept under their original
names: `HSP90AA2P`, `MKRN4P`, `PTTG3P` (pseudogenes), `DLEU1` (lncRNA),
`MIR17HG`. All non-coding or pseudogenes, none of which could be in MitoCarta, so
nothing is lost. `MIR17HG` is the miR-17-92 host gene - the miRNA arm of the
Menegollo MB2_UF fork - appearing inside CollecTRI's MYC regulon.

## 9. G3 - PASS

`github.com/gszabadkai/Menegollo_Bentham` (private) carries
`Input and output data/TCGA analysis/TCGA.all.biclusters.RNAseq.Rdata`: one
data.frame, 849 x 213, keyed on TCGA patient barcodes.

`MB1.forkscale`, `MB2.forkscale`, `MB3.forkscale` are **complete, zero NAs across
all 849 samples**, alongside `.forkscale.log`, thresholded `.fork_0.4`/`.fork_0.8`
calls, `MB{1,2,3}.fork` membership, PAM50, purity, proliferation and histology.

Continuous forkscale is what Block D asks for, so no MCbiclust re-derivation is
needed and the fidelity risk in plan section 13 does not apply. **ED2 is cheap.**
This also supplies `MB1_forkscale` for the F3 incremental-value stop gate.

---

## 10. Infrastructure settled this session

Three provenanced inputs, each with a README (commit `f63f14d`):

| Input | Location | Note |
|---|---|---|
| Human MitoCarta 3.0 | `data/mitocarta_human/` | Broad workbook, MD5-verified, 1,136 genes / 149 pathways |
| CollecTRI | `data/collectri_human/` | OmniPath web service snapshot, SHA-256 recorded |
| Felsher signature | `data/genesets_from_library_human/` | 67 genes, tag-pinned from library `v1.0` |

**The library's `outputs/gmt/human/` tree was rejected and must stay rejected.**
It is mouse-native sets run through `mouse_to_human()`: human symbols, mouse
provenance. Its MitoCarta GMT carries 1,083 genes against the real 1,136, and its
Felsher set is a human-mouse-human round trip carrying 62 against 67, dropping
NPM1 and RPLP0 and adding the paralog artefact EIF5AL1. It is also gitignored in
the library, so the tag pins none of it. See `CLAUDE.md`.

`decoupleR::get_collectri()` does not work in this environment. OmnipathR 4.0.0
(Bioc 3.23) fails on an `ncbi_tax_id` join, after a cache wipe and with organism
passed explicitly. Not connectivity. Hence the snapshot. Both packages remain
installed and are still used for scoring.

## 11. Still open

- **D4 - Lee et al. 2017 scoping. BLOCKS G2.** Plan section 14 says do this before
  G2. If it pre-empts H1, H1 becomes a citation, H2 becomes the novel mechanism
  and H4 the novel consequence.
- **D5** - primary neoadjuvant cohort for H4. Plan recommends GSE25066 primary,
  BrighTNess as subtype-matched replication. Unconfirmed. Does not block G2.
- **D7** - above. Confirm before script 09.
- **G1 correlation criterion** - needs script 01.
- **Script 07 specificity panel source.** FAO, one-carbon, mitoribosome, TCA, ROS
  defence would naturally have come from the library's human GMT tree, which is
  now rejected. Source undecided. Building them from MitoCarta pathways directly
  is the obvious candidate and would keep provenance clean, but it has not been
  decided. Do not default to the library.
