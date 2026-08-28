---
date: 2026-08-28
status: D8 and D9 RESOLVED
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 8 and 10)
  - 2026-08-28_D7_proliferation_covariate.md (the other Block C specification decision)
  - scripts/03_build_covariate_table.R (where both were raised, unsigned)
decides:
  - D8 - missing data. COMPLETE CASES, n = 938, no imputation. Two named
    sensitivity refits, and a generalisability comparison of the 157 dropped.
  - D9 - plate. Fixed factor with plates of n < 10 pooled into `other`,
    41 levels -> 28. Random intercept as a sensitivity on the primary arm only.
next-action: script 08, then script 09
---

# D8 and D9 - the two Block C specification decisions

Both were raised in the header of `scripts/03_build_covariate_table.R` and left
unsigned. Resolved 2026-08-28, **before script 09 exists, before script 08
exists, and before any model has been fitted.** Neither decision touches the
hypotheses, the endpoint, or the exposure; both are about how the nuisance
structure of TCGA-BRCA is handled.

The numbers below are from `results/tcga_brca_covariates.rds` and are structural
- counts of who has which assay, and how samples were plated. No outcome data
was read to reach either decision.

---

## D8 - missing data

### 1. What is actually missing

```
purity              75 missing (6.8%)     PAM50   114 (10.4%)
leukocyte_fraction  16        (1.5%)      TP53     78  (7.1%)
plate                0        (0.0%)

complete cases on the Block C covariate set: 938 of 1,095 (85.7%)
```

**The missingness stacks, and that is the fact that decides this.** Of the 157
patients dropped, only 52 are missing one variable; 84 are missing two and 21
are missing three. So trimming the model to buy sample size does not work:

```
drop purity              -> 966 complete (+28)
drop leukocyte_fraction  -> 951 complete (+13)
drop PAM50               -> 949 complete (+11)
drop TP53_status         -> 938 complete  (+0)
drop plate               -> 938 complete  (+0)
```

Every TP53-missing patient is already missing something else. The 157 are one
coherent group of poorly-characterised samples, not a per-variable loss, and no
single covariate is the culprit that could be sacrificed.

### 2. Decision: complete cases, n = 938, no imputation

Three reasons, in order of strength.

**(a) Imputing PAM50 would be circular.** PAM50 is called from expression. The
MYC and OXPHOS exposures are scored from the same matrix. Imputing the covariate
that adjusts the exposure, from the exposure's own data, is not a repair - it
manufactures exactly the structure the adjustment is supposed to break. This
argument alone is sufficient and it is why multiple imputation is rejected here
rather than merely not chosen.

**(b) Complete-case regression is unbiased under the missingness that is
present.** A complete-case fit is unbiased when missingness is independent of
the outcome given the covariates in the model. Here it is driven by whether
ABSOLUTE, PAM50 and MC3 were run and callable - sample purity, coverage, and
assay availability - not by `BBC3` or `BCL2L1`. Critically, **RNA is complete
for all 1,095 patients, so neither the exposure nor the endpoint is ever
missing.** The entire loss is in covariates.

**(c) Plate cannot be imputed at all.** It is a design fact recorded in the
barcode, not a latent quantity. A missing plate would mean a missing barcode.

### 3. The sensitivity ladder, corrected

Script 03's header proposed "a pre-specified sensitivity refit on all patients
using only fully-observed covariates". **That sentence degenerates**: plate is
the only fully-observed covariate, so the fit it describes is

```
PRIME ~ MYC * OXPHOS + plate        n = 1,095
```

That is still worth running - it shows the interaction is not manufactured by
the complete-case restriction - but it must be named for what it is rather than
described as a covariate-adjusted model. A rung is added between:

| Fit | Covariates | n |
|---|---|---|
| **M1 primary** | full Block C set | 938 |
| **M2 sensitivity** | drop purity, PAM50, TP53; keep leukocyte fraction + plate | 1,079 |
| **M3 sensitivity** | plate only | 1,095 |

M2 costs 1.5% of the cohort and retains the stromal covariate that matters most
in breast. M3 is the unadjusted-except-batch bound.

### 4. What must be reported

- **n alongside every estimate.** Not in a footnote.
- **A comparison of the 157 dropped against the 938 kept, on M-a and on PRIME.**
  This is a **generalisability** statement, not a bias check, and it must be
  labelled that way. Complete-case analysis stays unbiased for the conditional
  model even if the two groups differ; what differing groups cost is the claim
  that the estimate describes all of TCGA-BRCA.
- **Whether missingness is concentrated in particular plates.** If it is, the
  complete-case restriction interacts with the batch structure and D9's pooling
  is doing more work than it appears to. One cross-tab, printed in script 09's
  QC block.

---

## D9 - plate

### 5. The structure

Within the 938 complete cases:

```
41 plate levels    median 15 samples, range 1 - 82
 4 plates of n = 1
14 plates with n < 10, covering 65 patients (7%)
plate x PAM50: 27% of the 205 cells empty
```

**The obvious alternative does not exist.** Tissue source site (barcode
positions 6-7) gives 36 levels against plate's 41 - barely coarser - and
**36 of the 41 plates span more than one TSS**. Plate and TSS are crossed, not
nested. There is no cleaner batch handle in the barcode, and substituting TSS
would trade a known processing batch for a differently-shaped one with the same
sparsity problem.

### 6. Decision: fixed factor, plates of n < 10 pooled into `other`

41 levels become 28, costing 27 degrees of freedom instead of 40. Every
surviving level has at least 10 observations, which is enough to estimate an
intercept; the four singleton plates stop consuming a parameter they cannot
inform.

**The threshold is set on design grounds and fixed before any model is fitted.**
A factor level needs enough observations to estimate a mean; 10 is the round
number at which that is true here. It is not tuned to any result, and the
alternative thresholds are recorded so the choice is auditable: n < 5 would pool
6 plates / 10 patients, n < 15 would pool 20 plates / 131 patients.

### 7. Why not a random intercept as primary, when it is the better model

It is the better model. A random intercept costs one variance parameter instead
of 27, shrinks small plates toward the grand mean rather than fitting them
exactly, and keeps the singletons. The reason it is not primary is compute, and
the number is not marginal:

> Script 09 fits the specificity battery against its expression-matched null:
> **17 arms x 2,000 draws = 34,000 null model fits per instrument**, before the
> three D7 specifications and the strata, and there are two co-primary
> instruments. With `lm` that is minutes. With `lmer` it is hours per instrument,
> and the null is not optional - it is what makes the battery a test rather than
> a ranking.

**So: random intercept as a pre-specified sensitivity on the primary arm only**
- OXPHOS subunits, both instruments, no battery, no null. If the `MYC:OXPHOS`
coefficient is stable across the two plate specifications, plate handling is not
load-bearing and that is one sentence in Methods. If it is not stable, that is a
finding about the data and must be reported rather than resolved by choosing.

### 8. Two consequences, not further choices

**(a) The same plate specification in the observed model and in all 34,000
nulls.** A percentile is only calibrated if the null models *are* the model. The
pooling map is computed once, on the primary complete-case analysis set, and
reused unchanged for every arm and every draw.

**(b) In a stratified refit, the same rule is re-applied within the stratum**,
after `droplevels()`. A plate with two observations inside the TP53-mutant
stratum is the same problem as a plate with two observations overall, and the
strata carry no null of their own, so nothing needs to stay comparable to the
full-cohort mapping. This is stated so that the two behaviours are not read as
an inconsistency.

---

## 9. Why neither of these is drift

Same test as D5 and D7, and it should be applied every time:

- both were **raised by the pipeline's own QC** (script 03's coverage report),
  which is what that report is for;
- both are settled **before script 09 exists** and before any model has been
  fitted;
- both were decided from **structural counts only** - who has which assay, how
  samples were plated - with **no outcome data read**;
- neither changes the hypothesis, the endpoint, the exposure, or the estimator.
  They change how a nuisance is handled;
- in both cases the **discarded alternative is retained and reported** (M2/M3
  for D8, the random intercept for D9), so the decision is auditable against
  what it displaced.

What the pre-registration prohibits is swapping a specification after seeing
that the first one failed. Nothing has been fitted. Dated so the ordering is
checkable.
