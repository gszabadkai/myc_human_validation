---
date: 2026-08-29
status: DRAFT - awaiting sign-off before any code is written
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 2, 7, 8, 10)
  - 2026-08-28_D7_proliferation_covariate.md
  - 2026-08-28_D8_D9_missing_data_and_plate.md
  - 2026-08-28_specificity_panel_proposal.md
  - results/tcga_brca_mito_scores.rds, results/tcga_brca_priming.rds
decides: nothing yet. This is the build spec for scripts/09_interaction_models.R.
---

# Script 09 - build spec

Block C, the primary interaction test. Everything it needs now exists: 07 has
the exposures on both co-primary instruments plus 17 arms x 2,000 matched null
sets, 08 has the endpoints, and D7/D8/D9 fix the specification.

This document fixes the model space **before** any coefficient is computed. The
reason is arithmetic: the dimensions below cross to tens of thousands of
possible fits, and a script written without a fixed spine is a machine for
finding whichever one confirms.

---

## 1. The design principle: a spine, and one-at-a-time excursions

There are **two primary fits and no others**. Everything else varies exactly one
dimension from a spine and is reported alongside it, never selected between.

```
SPINE, one per co-primary instrument:

  PRIME ~ M_a * OXPHOS
        + purity + leukocyte_fraction + PROLIF_DISJOINT
        + PAM50 + TP53_status + plate_pooled

  OXPHOS = `OXPHOS subunits`, GSVA on VST      -> spine G
  OXPHOS = `OXPHOS subunits`, mitoPPS          -> spine M

  n = 938 complete cases (D8), lm, plate pooled to 28 levels (D9)
```

**The primary test is the `M_a:OXPHOS` coefficient, two-sided, on both spines.**
Per the co-primary rule fixed before any model was fitted: report both, claim
only what both support. An effect on one instrument alone is instrument-
dependent and is not a positive result.

### 1.1 Both exposures are z-scored, and the constants are frozen

`M_a` and each OXPHOS arm are centred and scaled **on the primary 938-patient
set**, and those same constants are reused in every stratum, sensitivity and
null fit. Two reasons, both structural:

- GSVA scores run roughly -1 to 1 while mitoPPS sits near 1 with a quite
  different spread, and the 18 arms differ in spread among themselves. Raw
  coefficients would not be comparable across arms or instruments, and the
  specificity forest - which is Panel a - would be meaningless.
- Freezing the constants means every coefficient in the paper is per SD of the
  same reference set, so a stratum estimate can be read against the full-cohort
  one.

Each null gene set is z-scored using **its own** mean and SD on the 938, because
it is a different variable. Not doing so would make the null a test of scale.

---

## 2. The excursion ladder

Each row changes one thing from the spine. Counts are fits, and both instruments
are always run.

| # | Dimension varied | Variants | Fits |
|---|---|---|---|
| 0 | **spine** | - | 2 |
| 1 | D7 proliferation spec, M-a | S2 (no term), S3 (`PROLIF_STD`) | 4 |
| 2 | MYC estimator | M-b with `PROLIF_STD` (S1=S3 for M-b), M-b with S2 | 4 |
| 3 | MYC instrument | M-c = `MYC_amp` (+2 GISTIC), available for all 938 | 2 |
| 4 | Missing data (D8) | M2 (n = 1,079), M3 (n = 1,095) | 4 |
| 5 | Plate (D9) | random intercept, primary arm only | 2 |
| 6 | Purity | ABSOLUTE > 0.7, **n = 270** | 2 |
| 7 | **Specificity battery** | 17 comparator arms in place of OXPHOS subunits | 34 |
| 8 | Endpoint negatives | BID, BAX, BCL2L11, BAK1 over BCL2L1 | 8 |
| 9 | Whole-family index | `PRIME_INDEX` (see 5) | 2 |
| 10 | **Limb-wise** (see 4) | BBC3, BCL2L1, and each negative's numerator | 12 |
| 11 | ER handling (see 6) | ER-adjusted PRIME + 4 negatives; within-ER PRIME | 14 |
| 12 | Strata (see 7) | TP53 x2, PIK3CA x2, PI3K-intact x2, PAM50 x4 | 20 |

**About 110 observed fits.** Plus the null, which is the expensive part:

```
17 arms x 2,000 draws x 2 instruments = 68,000 fits, spine specification only
```

---

## 3. The expression-matched null, and how a percentile is computed

Fixed here, before any coefficient exists.

For arm A on instrument I, let `b_obs` be the `M_a:A` coefficient from the spine
specification, and `b_null[1..2000]` the same coefficient refitted with each
matched draw substituted for A.

```
percentile(A) = 100 * mean(b_null < b_obs)
p_emp(A)      = 2 * min( mean(b_null <= b_obs), mean(b_null >= b_obs) )
```

Two-sided, matching plan section 10's two-sided primary test. With 2,000 draws
the smallest reportable value is 0.001; report it as `p < 0.001` rather than 0.

**The null models are the model.** Same covariates, same plate map, same 938
patients, same frozen z-scaling constants. A percentile is only calibrated if
the only thing that changed is the gene set (D9 consequence (a)).

`mtDNA-encoded OXPHOS` gets a coefficient and **no percentile** - all 13 of its
genes sit in the top expression ventile and no ventile-matched draw stands in
for them. Already settled; restated because a blank cell invites someone to
fill it.

### 3.1 The paired null - the comparative claim, tested

A one-set null answers "is OXPHOS extreme?". The claim is comparative -
"the respiratory chain, and not its own assembly machinery" - so the null has to
be on the **difference**, with both sets redrawn. The mouse did exactly this
(`myc_mouse/scripts/43`, `paired_null`).

```
d_obs    = b_obs(A) - b_obs(B)
d_null_i = b_null_i(A) - b_null_i(B),   i = 1..2000
```

This costs **no extra fitting** - the 68,000 coefficients above already contain
everything. It is valid because the draws for A and B are generated from
independent seeds, so their difference is a proper null for the difference of
two independently expression-matched sets.

Three pre-specified contrasts, all against `OXPHOS subunits`:

1. **vs `OXPHOS assembly factors`** - the primary. Same complexes, same
   umbrella, different function.
2. vs `Mitochondrial ribosome` - the growth/translation alternative.
3. vs `Nucleotide metabolism` - the proliferation alternative.

**Pre-stated, because the run of 07 makes it necessary.** `OXPHOS subunits` and
`OXPHOS assembly factors` correlate at rho = 0.875 as sample scores. So:

- interaction present for subunits and **absent** for assembly factors, despite
  that correlation -> a strong result;
- present for **both** -> uninformative about specificity, and must be reported
  as such rather than as double confirmation.

Same for the five per-complex arms, which correlate 0.81-0.97 with the claim arm
and with each other. **They are one compartment seen five ways, not five
independent confirmations**, and the text must say so. `CII subunits` is the one
outlier (0.56-0.67) and is a 4-gene set; if it behaves differently, small-set
noise is the first explanation, not biology.

---

## 4. Limb-wise fits - because PRIME is not the ratio it looks like

From script 08:

```
                    vs numerator   vs denominator   sd(numerator)
PRIME                    0.838         -0.128           1.081
BID/BCL2L1               0.766         -0.622           0.670
BAX/BCL2L1               0.721         -0.430           0.704
BCL2L11/BCL2L1           0.709         -0.679           0.585
BAK1/BCL2L1              0.743         -0.595           0.660
                                        sd(BCL2L1) =    0.559
```

`BBC3` varies about twice as much as `BCL2L1`, so **PRIME is close to a BBC3
readout while its four negative controls are genuine ratios.** A PRIME-versus-
negative difference could come from that structural asymmetry rather than from
biology.

So the spine is additionally refitted with `log2(BBC3)` alone and `log2(BCL2L1)`
alone as the endpoint, and with each negative's numerator alone. This turns a
hidden confound into a reported quantity: if the interaction lives entirely in
the `BBC3` limb, that is the PUMA result the mouse predicts; if it lives in the
`BCL2L1` limb, PRIME is measuring the guardian, not the trigger.

---

## 5. `PRIME_INDEX` is a different endpoint, not a robustness check

Also from 08: the index correlates with PRIME at **rho = 0.312** (equal-weight
variant 0.425, with-BCL2 0.29). The weighting explains it - `MCL1` is 75.6% of
the anti-apoptotic sum and `BBC3` is 4.8% of the pro-apoptotic sum - so as
specified the index is an `MCL1` versus `BCL2L11/BAX/BOK` contrast that happens
to contain two of PRIME's sixteen genes.

It stays exactly as pre-registered. What changes is the **label and the
inference**: it is a whole-family priming index, reported in its own right.
**Disagreement between it and PRIME is not a robustness failure**, because at
rho = 0.31 the two were never measuring the same thing. Stated now so it cannot
be read either way after the fact.

---

## 6. ER handling

08 found that the ER confound hits the endpoint negatives and not PRIME, and the
mechanism is now known:

```
BBC3 +0.231, BCL2L1 +0.259  ->  PRIME       +0.091   (cancels)
BID  -0.396                 ->  BID/BCL2L1  -0.452   (adds)
BAK1 -0.422                 ->  BAK1/BCL2L1 -0.474   (adds)
```

PRIME is ER-clean by coincidence of its two genes moving together on ER; BID and
BAK1 move opposite to the shared denominator, so their ratios accumulate ER
signal. **A PRIME-versus-negative difference could be an ER difference.**

`PAM50` is already in the spine and carries most of the ER information, so the
negatives are partly adjusted. On top of that, and pre-specified:

- an **ER-adjusted** variant (`+ er_call`) for PRIME and each negative -
  costs 46 patients, **n = 892**;
- a **within-ER-status** refit of PRIME (ER+ n = 691, ER- n = 201).

---

## 7. Strata, with their sizes stated before they are fitted

Within the 938:

```
TP53              mutant 325   wildtype 613      H3
PIK3CA_altered      TRUE 328      FALSE 610      H2
PI3K_pathway_intact TRUE 549      FALSE 389      H2, the better split
PTEN_altered        TRUE  83      FALSE 855      H2
PAM50   LumA 478  LumB 195  Basal 159  Her2 77  Normal 29
purity > 0.7        n = 270
```

Consequences, fixed now rather than discovered later:

- **`PTEN_altered` is not fitted as its own stratum** at n = 83. PTEN enters H2
  only through `PI3K_pathway_intact`, which is the pre-specified split.
- **`PAM50 = Normal` (n = 29) is not fitted.** Four PAM50 strata, not five.
  `Her2` at 77 is reported with its n against every estimate and read as
  indicative.
- **The purity-high sensitivity is n = 270.** For a three-term interaction that
  is weak, and a null there is uninformative rather than contradictory. Said
  now, because after the fact it would sound like an excuse.
- D9 applies within every stratum: re-pool plates at n < 10 **within the
  stratum**, after `droplevels()`.

---

## 8. FOXO3 and H2

`FOXO3_activity` correlates **-0.508 with purity** and **+0.374 with leukocyte
fraction** - the strongest confound anywhere in this pipeline, and H2 rests on
it. So, pre-specified:

- every FOXO3 model carries purity and leukocyte fraction, and
- **the purity-high refit is mandatory for any FOXO3 claim, not optional.** An
  unadjusted FOXO3 result is not reportable.

H2's own test - FOXO3 activity against OXPHOS within the PI3K-intact stratum -
is fitted here; its interpretation waits on the H1 result.

---

## 9. Multiplicity - what is corrected and what is not

- **The primary test is one coefficient per instrument.** Pre-registered, not
  corrected, and the co-primary rule is stricter than any correction: it must
  hold on both.
- **The specificity battery is not a family of hypotheses.** It is a set of
  controls whose *pattern* is the evidence, and the expression-matched null is
  what calibrates it. Bonferroni over 17 arms would be answering a question
  nobody asked.
- **Every sensitivity is reported, never selected.** The D7 three-specification
  panel, the D8 ladder and the D9 plate variants all appear side by side, as
  part of the result rather than as an appendix.
- No result is quoted from a single cell of this table without the pattern
  around it.

---

## 10. The decision rule for H1, fixed before the fits

H1 is supported only if **all** of the following hold:

1. `M_a:OXPHOS` on PRIME is **positive with a 95% CI excluding zero on both
   instruments**;
2. `OXPHOS subunits` beats its expression-matched null (`p_emp < 0.05`) on both;
3. the **paired** contrast against `OXPHOS assembly factors` beats its paired
   null on both;
4. the endpoint negatives do not show the interaction, or show it materially
   weaker, and the limb-wise fits show the effect is not carried by the
   `BCL2L1` limb alone.

Anything short of that is reported as what it is. The full pattern is reported
either way - plan section 15 already says what a well-formed failure looks like,
and it does not require a fifth hypothesis.

The `BUFFER` clause of H1 is Block B and belongs to script 10, not here.

---

## 11. Implementation notes that change the runtime by an order of magnitude

- **Do not call `lm()` with a formula inside the null loop.** Formula parsing and
  `model.frame` dominate and would turn a one-minute job into an hour. Build the
  model matrix once per (instrument, spec, analysis set); only two columns change
  between draws - the OXPHOS main effect and its interaction with `M_a`. Swap
  those two columns and use `.lm.fit`, extracting the interaction coefficient
  only. Standard errors are needed for the ~110 observed fits, not for the
  68,000 null fits.
- Assert once, loudly, that the swapped design matrix has the column order the
  coefficient index assumes. An off-by-one here silently reports the wrong
  coefficient 68,000 times.
- The null draw scores come from `results/mito_null/`, walked through
  `null_manifest`. Nothing is re-scored.
- `set.seed` is not needed - 09 draws nothing. All randomness lives in 07 and is
  already frozen in the cache.

### Outputs

```
results/block_c_models.rds
  $coefficients   tidy: endpoint, myc, arm, instrument, spec, missing, stratum,
                  n, term, estimate, se, ci_lo, ci_hi, p
  $null           arm, instrument, b_obs, null_median, percentile, p_emp
  $paired_null    arm_a, arm_b, d_obs, percentile, p_emp
  $decision       the four clauses of section 10, evaluated
outputs/tables/   the same as csv
```

No figure. Panel a is drawn in script 18 from `$coefficients` and `$null`.

---

## 12. What script 09 does not do

- Block B / `BUFFER` (script 10), `STATE` (11), outcome models (13).
- Nothing with the neoadjuvant cohorts, METABRIC or SCAN-B.
- No new scoring of any kind. Exposures and endpoints are consumed as built.
