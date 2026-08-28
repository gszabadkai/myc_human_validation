---
date: 2026-08-28
status: D7 RESOLVED
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 8, 10)
  - 2026-08-28_G1_result_and_decisions.md (section 6, where D7 was raised)
decides:
  - D7 - RESOLVED. Estimator-specific handling. M-a gets a disjoint proliferation
    covariate in the primary model; M-b keeps the standard one. All three
    specifications reported for both.
next-action: script 01 (discharges G1's correlation criterion)
---

# D7 - the proliferation covariate overlaps the MYC exposure

Raised by gate G1 (see that note, section 6). Resolved 2026-08-28, **before script
09 exists and before any model has been fitted.**

---

## 1. The problem

Plan section 8 specifies a proliferation covariate - Hallmark `E2F_TARGETS` +
`G2M_CHECKPOINT` GSVA - and section 10 puts it in the Block C model:

```
PRIME ~ MYC * OXPHOS + purity + leukocyte_fraction + proliferation
        + PAM50 + TP53_status + plate
```

G1 found the covariate and the exposure share genes:

| Estimator (MitoCarta-stripped) | shared with covariate | share of estimator |
|---|---|---|
| `FELSHER` (M-a, primary) | 9 of 61 | **14.8%** |
| `COLLECTRI_MYC_ALL` (M-b) | 80 of 811 | 9.9% |
| `COLLECTRI_MYC_STIM` | 78 of 739 | 10.6% |

**This is not confounding. It is overlapping measurement.** Roughly one in seven
genes constituting the MYC score also constitutes the covariate meant to adjust it,
so adjusting for proliferation partially adjusts away MYC itself and biases the
`MYC:OXPHOS` coefficient toward the null. No sample size fixes this - n ~1000 makes
it more precise, not less biased.

The enrichment itself is expected biology: MYC drives proliferation, and G1 confirmed
the overlap is far above chance (Felsher vs E2F 8.6x, p=6.9e-5). The problem is
structural, not biological.

## 2. The asymmetry that decided it

The obvious fix - rebuild the covariate from genes not in the estimator - looked
uniformly cheap. It is not. The cost is **wildly asymmetric between estimators**:

```
Covariate = union(HALLMARK_E2F_TARGETS 200, HALLMARK_G2M_CHECKPOINT 200)
          = 327 genes (73 shared between the two sets)
```

| Estimator | shared | cost to the EXPOSURE | cost to the COVARIATE | disjoint covariate |
|---|---|---|---|---|
| `FELSHER` | 9 | 14.8% | **2.8%** | 327 -> **318 genes** |
| `COLLECTRI_MYC_ALL` | 80 | 9.9% | **24.5%** | 327 -> **247 genes** |
| `COLLECTRI_MYC_STIM` | 78 | 10.6% | 23.9% | 327 -> 249 genes |

For **Felsher**, the disjoint covariate drops 9 genes from 327. Essentially free, and
what remains is still a perfectly good proliferation score.

For **CollecTRI**, it drops a quarter of the covariate - and exactly the wrong
quarter. The removed genes include `AURKA`, `BIRC5`, `CDK1`, `CCNB2`, `CDC20`,
`CHEK1`, `MYBL2`, `PCNA`. A proliferation score stripped of the canonical
proliferation machinery is a degraded construct, and adjusting for it would be worse
than not adjusting at all.

That difference is real biology rather than an artefact: CollecTRI's 811-gene MYC
regulon genuinely contains most of the cell-cycle machinery, because MYC drives the
cell cycle. Incidentally it is a further point in favour of the D2 decision - **M-a is
the less entangled estimator on this axis too**, as it was on mitochondrial
circularity.

## 3. Decision - estimator-specific

**M-a (`FELSHER`, primary): the disjoint covariate goes in the PRIMARY model.**
A covariate that eats 15% of the exposure is a specification error, not a
conservative choice, and here the fix costs 2.8% of the covariate.

**M-b (`COLLECTRI_MYC_ALL`, concordance check): keep the standard covariate.**
The disjoint version is not viable. Rely on the no-covariate sensitivity instead.

**Both: report all three specifications, always.** If the interaction survives all
three, the result is robust to the construction. If it appears only without the
covariate, that is informative and must be **reported as such, never selected**.

### The three specifications

| Label | Covariate | Applies to |
|---|---|---|
| **S1 primary (M-a)** | `PROLIF_DISJOINT` (318 genes) | M-a only |
| **S1 primary (M-b)** | `PROLIF_STD` (327 genes) | M-b only |
| **S2 sensitivity** | none - proliferation term dropped | both |
| **S3 sensitivity** | `PROLIF_STD` (327 genes) | both |

For M-a, S3 is the plan's original pre-specified model, retained so the amendment is
auditable against it. For M-b, S1 and S3 are the same fit and are reported once.

### Definitions, fixed

```
PROLIF_STD      = GSVA over union(HALLMARK_E2F_TARGETS, HALLMARK_G2M_CHECKPOINT)
                  = 327 genes

PROLIF_DISJOINT = PROLIF_STD minus the 9 genes shared with the MitoCarta-stripped
                  Felsher estimator = 318 genes
```

The nine removed genes, listed so this is reproducible without re-running G1:

```
CTPS1  DNMT1  HMGA1  NCL  NOP56  PRMT5  RANBP1  TFDP1  UCK2
```

(`HMGA1` is in both Hallmark sets; the other eight appear in one each.)

**`PROLIF_DISJOINT` is defined against the stripped Felsher set, not the raw 67.**
It is therefore specific to M-a and must never be used with M-b - a different
estimator implies a different disjoint set, and silently reusing this one would
reintroduce the circularity it exists to remove.

## 4. Why this is not drift

This is the **second** primary-model amendment made on 2026-08-28, after the
section 7.5 amendment under D5. That pattern deserves stating plainly rather than
leaving for a reviewer to notice.

Both amendments meet the same conditions, and the conditions are the point:

- raised by a **gate**, whose entire purpose is to find this class of problem before
  models are built;
- made **before script 09 exists** and before any model has been fitted;
- made with **no outcome data seen** - the audit is over gene-set membership only;
- they **improve measurement of the same construct** rather than changing the
  hypothesis, the endpoint, or the estimator;
- the original specification is **retained and reported** (S3), so the amendment is
  auditable against it.

What the pre-registration prohibits is swapping a specification *after* seeing that
the first one failed. Nothing here has been fitted. Recorded with a date so the
ordering is checkable.

## 5. Consequences for the build

- **Script 07** scores both `PROLIF_STD` and `PROLIF_DISJOINT`. One extra GSVA run,
  same VST input object, same cohort-relative rules.
- **Script 09** fits S1/S2/S3 and reports the `MYC:OXPHOS` coefficient from each side
  by side. The three-specification panel is part of the result, not an appendix.
- **Scale discipline unchanged**: both covariates are GSVA and therefore want
  log-scale VST input, `kcdf = "Gaussian"`. They share an input object with each
  other, never with mitoPPS.
- **Cohort-relativity unchanged**: score both covariates within each cohort's own
  run. Do not carry a disjoint gene list across cohorts and assume comparability of
  the resulting scores.

## 6. What this does not fix

Stripping MitoCarta genes did not reduce the proliferation enrichment (G1 note
section 4: fold enrichment was flat or slightly higher afterwards, because the
denominator shrank), and neither does this. **The Felsher signature remains
enriched for proliferation genes above chance, and that is expected biology.**

D7 removes the *circularity* between exposure and covariate. It does not, and
cannot, separate MYC activity from proliferation as biological constructs. If the
`MYC:OXPHOS` interaction is real, some of it will run through proliferation, and
that is a limitation to state in Methods rather than a defect to engineer away.
