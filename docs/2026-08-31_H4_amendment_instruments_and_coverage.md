---
date: 2026-08-31
status: AMENDMENT - pre-data, written from script 12's marginals and nothing else
amends:
  - 2026-08-30_STATE_frozen_and_H4_buffer_declaration.md section 6.5 (two instruments)
relates-to:
  - 2026-08-28_D5_cohort_selection.md (cohorts, meta-analysis, power)
  - 2026-08-27_human_validation_plan.md (sections 2, 3, 7.5; H4)
  - data/neoadjuvant/README.md (provenance, landmines)
decides:
  - A1. H4 is GSVA-ONLY. The two-instrument rule is suspended FOR H4 ONLY, on
    feasibility grounds. mitoPPS does not exist in any of the three cohorts.
  - A2. Each cohort is scored on the genes it has, coverage reported beside
    every estimate, with the 3-cohort intersection as a pre-specified
    sensitivity. GSE25066 is excluded from the SPECIFICITY battery only.
  - C3. The score x treatment contrast against I-SPY2's control arm is
    estimable only within HRpos_HER2neg and TNBC.
  - A4. The F1 covariate set is subtype + treatment. Stage and grade DO NOT
    EXIST in the primary cohort; the plan's F1 model cannot be fitted as
    written. Cohort-specific covariates are a sensitivity only.
  - A5. Missing values: drop genes above 5% NA, impute the rest at the gene
    median within cohort. Declared before any score was computed.
next-action: script 13. Once it computes the first score-versus-pCR
  association, NOTHING in this note or in the 2026-08-30 declaration may be
  amended again.
---

# H4 amendment - instruments, coverage, and what the treatment contrast can reach

Written after `scripts/12_fetch_neoadjuvant_cohorts.R` was run on 2026-08-31 and
before any score has been computed in any of the three cohorts.

**Every fact below is exposure-side.** Gene names against gene names, scales,
and arm-by-subtype counts. Script 12 computes no score and no score-versus-pCR
association, by construction and by D5 section 2's standing rule. The pCR
marginals it prints were already published in D5. So nothing here could have
been tuned to a result, and this is the last moment at which any of it can
honestly be written.

---

## 0. What makes an amendment legitimate here, and what would not

An amendment is legitimate when it is forced by a **fact about the
instrument or the data** that was not knowable earlier, and is made before any
outcome association exists. It is illegitimate when it is prompted by a result,
however framed.

A1 and A2 are of the first kind: the instrument does not exist in these
cohorts, and the platforms do not carry the genes. Both were unknowable before
the files were parsed, and both were surfaced by guards written into script 12
before it was ever run.

**What would NOT be legitimate, and is not done here:** relaxing anything in the
2026-08-30 declaration's section 6.5 directions, section 6.6 informative-failure
criterion, or section 6.1 BUFFER definition. Those are outcome-side and they
stand exactly as written. See section 4.

---

## 1. AMENDMENT A1 - H4 is GSVA-only

### The fact

All three cohorts are deposited on the log2 scale:

```
cohort        min      max    median   frac_neg    frac_NA   genes  samples
GSE194040   2.335   21.250     7.504      0.000   1.50e-04   19134      988
GSE164458   0.000   16.431     1.561      0.000   0.00e+00   24613      482
GSE25066   -5.007   20.861     8.589   9.32e-05   0.00e+00   12502      508
```

GSVA is satisfied - log scale, `kcdf = "Gaussian"`. **mitoPPS is not available
in any of them.** It wants linear DESeq2-normalised counts (CLAUDE.md). Two of
these are microarrays and the third is deposited already logged; no linear
count matrix exists for any.

**Exponentiating does not recover it.** `2^x` gives a linear scale
arithmetically, but not the quantity mitoPPS was defined on. mitoPPS takes
ratios of pathway MEANS ACROSS GENES, which assumes a cross-gene comparability
that probe affinity does not provide, and I-SPY2's matrix was additionally
ComBat-adjusted on the log scale, so `2^x` there is a batch-corrected intensity
rather than an abundance.

### The amendment

> The 2026-08-30 declaration section 6.5's "on both instruments" is **suspended
> for H4 only**. H4's primary and secondary tests are computed on **GSVA
> alone**.

The two-instrument rule remains in force everywhere else in this arm. It is
what Blocks B, C and G were held to and nothing about H4 changes that.

### What this costs, stated rather than absorbed

H4's evidential standard is now **weaker than every other block's**. Blocks C, B
and G each required agreement between two instruments that demonstrably can and
do disagree - the surviving Block C finding cleared that bar, and the DepMap BIM
result failed it at n = 51. H4 has one instrument and therefore cannot clear it.

**A positive H4 must be reported as single-instrument** and must not be written
as though it met the same standard as the Block C result. A null H4 is likewise
single-instrument, and D5 section 6.3 already requires it to be reported as
"not powered to exclude a modest effect".

### The move that is forbidden, named now so it is not made later

Adding ssGSEA, a mean-z of the gene set, or any other level metric computed on
the same log matrix **does not restore two instruments.** The rule exists
because GSVA measures LEVEL and mitoPPS measures COMPOSITION, and they answer
different questions - which is precisely why the same genes can move opposite
ways on them (CLAUDE.md). Two level metrics on one matrix are one instrument
reported twice, and presenting them as two would be worse than admitting one,
because it would look like corroboration.

---

## 2. AMENDMENT A2 - how OXPHOS is scored when the platform is short of genes

### The fact

Set coverage, measured against each cohort's own gene universe:

```
set                       full  TCGA  I-SPY2  BrighTNess  Hatzis   int3
OXPHOS subunits             89    88      69          70      57     57
OXPHOS assembly factors     68    68      57          65      34     33
Felsher M-a (stripped)      61    61      57          58      44     44
```

`int3` is the 3-cohort intersection (11,676 genes). GSE25066 is the bottleneck:
its GPL96 probe collapse yields 12,502 genes, against 19,134 and 24,613.

### The amendment

> **Primary:** each cohort is scored on the genes that cohort has. The coverage
> fraction is reported beside every per-cohort estimate.
>
> **Sensitivity, pre-specified:** the same models on the 3-cohort intersection
> sets, so that a set-composition artefact is visible rather than assumed
> absent.
>
> **Specificity battery:** a cohort may contribute to it only if the
> ASSEMBLY-FACTOR control reaches **0.80** coverage in that cohort.
> **GSE25066 (0.50) therefore does not**, and contributes to the main effect
> only.

### Why not harmonise to the intersection as primary

Harmonising would let the cohort D5 already **demoted to third** - single-arm,
2011 array, no BBC3 - define the exposure in the n = 988 primary cohort,
discarding 12 of I-SPY2's 69 OXPHOS genes to accommodate it. That is the tail
wagging the dog.

GSVA is rank-based and cohort-relative. A set measured on 57 of 89 genes rather
than 69 is a noisier estimate of the same latent quantity, not a different
construct, provided the missing genes are not systematically one complex - which
the intersection sensitivity is there to check. Reporting coverage beside each
estimate keeps the heterogeneity visible in the meta-analysis rather than
hidden inside it.

### Why the specificity carve-out, and why 0.80

The assembly factors are **the primary specificity control** in this arm, not a
leftover: CLAUDE.md's standing convention records that in the mouse they sit at
percentile 50.2 of expression-matched sets while the subunits sit at 0.0, and
being the same complexes they exclude what a distant pathway cannot. A control
measured at half coverage cannot do that work.

**0.80 is not invented here.** It is `MIN_SET_FRAC` from
`scripts/14_depmap_dependency.R`, this arm's existing convention for the same
question. Applying it to the specificity claim - where a measurable control is
the whole point - rather than as a blanket floor that would eject an entire
cohort from the main test, is the narrower rule and the one tied to what is
actually broken.

Consequence to carry into the text: **the H4 specificity battery rests on
I-SPY2 (57/68) and BrighTNess (65/68).**

---

## 3. RECORDED CONSTRAINT C3 - what the treatment contrast can actually reach

Not an amendment. The 2026-08-30 declaration section 6.7 explicitly left the
`score x treatment` interaction to script 13's build spec. This fixes the
constraint that build spec must respect, so the claim cannot quietly widen.

### The fact

I-SPY2's arm-by-subtype table:

- The **179-patient Paclitaxel control arm is 94 HRpos_HER2neg + 85 TNBC and
  contains ZERO HER2+ patients.**
- The HER2-targeted arms (Trastuzumab, Pertuzumab, T-DM1, and Neratinib in part)
  are almost entirely HER2+; the others are almost entirely HER2-negative.

So arm and subtype are close to collinear across much of the design.

### The constraint

> The `score x treatment` contrast against I-SPY2's control arm is estimable
> **only within HRpos_HER2neg and TNBC**. It is not estimable in either HER2+
> stratum, at any n, because those strata contain no control patients.

Where the randomised treatment evidence therefore comes from:

| Cohort | Contributes to `score x treatment`? |
|---|---|
| GSE194040 | yes, within HRpos_HER2neg and TNBC only |
| GSE164458 | yes - randomised, 3 arms (A 237 / B 122 / C 123), TNBC only |
| GSE25066 | **no** - single-arm by design (D5 section 4) |

D5 preferred I-SPY2 over GSE25066 partly because it "can test `score x
treatment`". That remains true and is now bounded: it tests it in two of four
receptor subtypes.

---

## 3b. AMENDMENT A4 - the F1 covariate set, which the plan's model cannot have

**Added 2026-08-31 while writing script 13, before it was run.** Still
exposure-side: this is metadata availability, not a result.

### The fact

Plan Block F1 specifies

```
pCR ~ MYC * OXPHOS * BUFFER + stage + grade + ER_status + arm
```

**Stage and grade do not exist in the primary cohort.** GSE194040's deposited
phenotype is `patient id, tissue, hr, her2, mp, pcr, arm` and nothing else - no
stage, no grade, no age, no nodal status. What each cohort actually carries:

| | subtype | treatment | stage | grade | nodal | age |
|---|---|---|---|---|---|---|
| GSE194040 | yes (hr/her2) | yes, 14 arms | **no** | **no** | **no** | **no** |
| GSE164458 | constant TNBC | yes, 3 arms | **no** | **no** | yes | **no** |
| GSE25066 | yes (IHC) | constant | yes | yes | yes | yes |

The only covariate available in all three is **subtype**, and it is constant in
GSE164458 by design.

### The amendment

> **Primary, harmonised across all three:**
> `pCR ~ MYC * OXPHOS * BUFFER_c + subtype + treatment`, logistic. A term that
> is constant within a cohort is dropped from that cohort's fit and the drop is
> reported.
>
> **Cohort-specific sensitivity:** each cohort's model additionally adjusted for
> what that cohort has - nodal and ECOG in GSE164458, stage/grade/age in
> GSE25066. Reported alongside, never instead.
>
> The meta-analysis pools the **harmonised** estimates, because pooling
> estimates from differently-adjusted models is not a meta-analysis.

`mp` is deliberately NOT used in GSE194040. D5 section 4 flagged it as probably
MP1-vs-MP2 within I-SPY2's high-risk criterion and said "confirm before use".
It has not been confirmed, so it stays out of the primary and is reported only
as an unconfirmed sensitivity.

**What this costs:** residual confounding by stage and grade in the primary
cohort cannot be excluded, and neither can it be quantified there. Subtype
absorbs part of it - `subtype` is the strongest single predictor of pCR in
these cohorts - but not all. State it as a limitation; do not imply the
plan's covariate set was fitted.

## 3c. AMENDMENT A5 - missing values, declared before scoring

### The fact

GSE194040 carries 2,835 NA values: 1,141 of 19,134 genes and 478 of 988
samples are touched, but only **one gene** exceeds 5% missing and the worst gene
is 6.3%. The other two cohorts have none.

The damage is concentrated where it matters: **12 of the 69 OXPHOS-subunit genes
present in I-SPY2 carry at least one NA**, and 7 of 57 assembly factors. So
"drop any gene with any NA" would cut I-SPY2's OXPHOS coverage from 69 to 57 -
down to GSE25066's level - to remove 2,835 values out of 18.9 million.

### The amendment

> 1. Drop any gene with **more than 5% NA** across that cohort's samples
>    (removes one gene in GSE194040, none elsewhere).
> 2. Impute the remaining NAs with **that gene's median across samples**, within
>    cohort.
> 3. Report the number of values imputed per cohort.

Gene-median imputation is the neutral choice for a rank-based score: GSVA
estimates each gene's distribution across samples and then walks a per-sample
rank, so an imputed sample sits at that gene's central rank contribution rather
than being pushed to either tail. It is not a claim about the missing values;
it is the choice that moves the score least.

---

## 4. What does NOT change

Stated explicitly, because an amendment note is exactly where scope creeps.

- **`BUFFER_c = mean(z(log2 MCL1), z(log2 BCL2L1))`**, within cohort, both genes
  also reported separately, BCL2 excluded. Unchanged.
- **Directions**: `MYC:OXPHOS:BUFFER_c` predicted **NEGATIVE**, `MYC:OXPHOS`
  positive at mean buffering. Unchanged.
- **The informative-failure criterion** in coefficient terms - the OXPHOS slope
  at `MYC = +1 SD` negative with CI excluding zero WHILE the three-way is null.
  Unchanged, and still the thing that separates an inversion from a null.
- **MYC = M-a; OXPHOS = `OXPHOS subunits`.** Unchanged.
- **Within-cohort z always; scores never pooled; estimates meta-analysed.**
  Unchanged (D5 section 6.2, CLAUDE.md).
- **STATE** as frozen in script 11, tertile BUFFER and the level-3-vs-4 contrast
  included. Untouched.
- **Subtype adjustment mandatory**, met by receptor subtype because PAM50 exists
  only in GSE25066 - already recorded in the 2026-08-30 note section 5 and in
  script 12.
- **Report CIs, not p-values alone** (D5 section 6.3). Unchanged.

## 5. Two further facts recorded, neither requiring an amendment

- **BBC3 is absent from GSE25066.** GPL96 carries exactly one BBC3 probe and it
  multi-maps to `MIR3190`/`MIR3191`, so the declared collapse rule drops it.
  `PRIME` is therefore not computable in that cohort. This does not affect H4,
  whose endpoint is pCR, nor the BIM replication, which needs BCL2L11 - present
  in all three.
- **I-SPY2 carries 0.015% NA values** (about 2,800 of 18.9 million). GSVA needs
  an explicit rule for them; script 13 must state one rather than inherit a
  default.

## 6. The audit trail

1. `scripts/12_fetch_neoadjuvant_cohorts.R` written and committed **before** the
   data was downloaded, with the coverage and scale checks already in it.
2. Data downloaded; checksums recorded in `data/neoadjuvant/README.md`.
3. Script 12 run 2026-08-31. Marginals only. **No score computed, no
   score-versus-pCR association computed.**
4. This note.
5. Only then: script 13.

Three separate guards in script 12 stopped the run on real silent-corruption
bugs before any of this was read - a merged sample column, a connection-buffer
overflow, and GEO's ragged characteristic lines mislabelling 198 of GSE25066's
508 endpoints. All three are fixed and recorded in the git history. The
coverage and scale numbers above come from the run that followed.

**After script 13 computes its first score-versus-pCR association, nothing in
this note or in the 2026-08-30 declaration may be amended.** If something is
wrong in either, it has to be fixed now.
