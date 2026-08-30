---
date: 2026-08-31
status: RECORDED - result note. Section 6 is a correction, recorded beside the
  registered readings rather than folded into them.
relates-to:
  - 2026-08-30_STATE_frozen_and_H4_buffer_declaration.md (the declaration)
  - 2026-08-31_H4_amendment_instruments_and_coverage.md (A1, A2, A4, A5, C3)
  - 2026-08-28_D5_cohort_selection.md (cohorts, power, what is detectable)
  - 2026-08-27_human_validation_plan.md (sections 2, 15; Block F1)
decides:
  - H4 is NOT SUPPORTED. Every substantive declared reading fails.
  - It fails as a NULL, not in the informative direction. Falsification
    criterion 4 is therefore NOT met, and the count stays at three of four.
  - The one nominal result in the block is MCL1 in one cohort, not buffering,
    with I2 = 82% across the three. NOT PURSUED.
  - All four hypotheses have now been tested and none is supported.
next-action: the BIM replication - but read section 8, it cannot be run in these
  cohorts as declared.
---

# Block F1 - H4 is not supported

Run 2026-08-31 via `scripts/13_outcome_models.R`, to the declaration signed off
2026-08-30 and the amendment of 2026-08-31, both written before any score met a
pCR column. Results in `results/h4_outcome_models.rds`.

n = 988 / 482 / 470, with 319 / 236 / 95 pCR events. (GSE25066 enters at 470 of
508: 20 samples have no pCR call and 18 no callable receptor subtype.)

---

## 1. Every substantive reading fails

| | Test | Predicted | |
|---|---|---|---|
| **F1-a** | meta three-way `MYC:OXPHOS:BUFFER_c` | NEGATIVE, CI excludes 0 | **FAIL** |
| F1-b | same sign in the primary cohort | negative in GSE194040 | pass - **see section 6** |
| **F1-c** | meta two-way `MYC:OXPHOS` | POSITIVE, CI excludes 0 | **FAIL** |
| F1-d | specificity: not the assembly factors | null / weaker | pass, **vacuously** |
| F1-e | A2 sensitivity agrees (int3 scoring) | same sign as F1-a | **pass, genuine** |
| **F1-f** | secondary `STATE` level 4 vs level 3 | level 3 higher pCR | **FAIL** |
| **F1-g** | informative failure (criterion 4) | slope negative WHILE three-way null | **FAIL** |
| F1-h | composite coheres with its limbs (6.3) | reporting, not a prediction | pass, **vacuously** |

**Four rows pass and three of the four are vacuous.** Only F1-e is a real pass,
and it is a sensitivity check confirming the result is not a set-composition
artefact - not support for anything.

## 2. F1-a - the primary test, and why 0.052 is not "nearly significant"

```
                        estimate      95% CI            p
per-cohort scoring
  GSE194040  n=988      -0.0082   -0.239 to  0.222    0.944
  GSE164458  n=482      +0.4066   +0.104 to  0.710    0.0085
  GSE25066   n=470      +0.2399   -0.128 to  0.608    0.201

  meta, fixed effect    +0.1629   -0.001 to  0.327    0.052
  meta, random effects  +0.1952   -0.069 to  0.459    0.148
  Q = 4.77, p = 0.092   tau2 = 0.0315   I2 = 58.1%
```

The prediction was **NEGATIVE**. The pooled point estimate is **positive**.

**The fixed-effect p of 0.052 must not be reported as a near-miss.** At
I2 = 58% the fixed-effect model is the wrong summary; the random-effects
estimate is +0.195 with a CI comfortably spanning zero. And the heterogeneity is
structured in the least favourable way possible: **the primary cohort - n = 988,
randomised, the one D5 selected as primary - is at exactly zero**, with a CI of
+/- 0.23. The entire pooled signal comes from the smallest cohort.

The int3 sensitivity reproduces all of it (-0.005 / +0.383 / +0.240; meta FE
+0.154). **F1-e passes: this is not an artefact of gene-set composition.**

## 3. The one nominal result is MCL1, in one cohort, and it is not buffering

This is what declaration section 6.3's mandatory per-gene companion was for.

```
GSE164458      composite  +0.4066   +0.104 to 0.710   p 0.0085
               MCL1       +0.4383   +0.189 to 0.688   p 0.00057
               BCL2L1     +0.0099   -0.200 to 0.220   p 0.926
```

BrighTNess's entire result is **MCL1 alone**. BCL2L1 contributes nothing; the
composite is a slightly diluted copy of a one-gene effect.

And that gene does not behave the same way anywhere else:

```
              GSE194040   GSE164458   GSE25066    meta RE      I2      Q p
MCL1            -0.057      +0.438      -0.001   +0.120     82.2%   0.0036
BCL2L1          +0.062      +0.010      +0.215   +0.075      0.0%   0.497
```

**MCL1's heterogeneity is strong enough to reject a common effect outright**
(Q p = 0.0036). Two larger cohorts put it at zero; one small one puts it at
+0.44.

Meanwhile **BCL2L1 - this arm's a priori focus, the PRIME denominator, the term
G2 identified as the one Lee et al. 2017 never touch - is homogeneously null.**
I2 exactly 0, fixed and random effects identical at +0.075.

### 3.1 Why the composite was never measuring "buffering" here

`BUFFER_c` was specified as the mean of two z-scores on the stated ground that
two anti-apoptotic proteins both being high is more buffering than one being
high (declaration 6.2). That argument assumes the limbs are related. **They are
not:** Spearman between `z_MCL1` and `z_BCL2L1` is 0.019, 0.120 and -0.016 in
the three cohorts.

So `BUFFER_c` averages two variables carrying different information rather than
two readings of one latent capacity. The specification stands as declared -
`max()` would have been worse, and the reasoning was fixed before any fit - but
the construct is weaker than it looked, and section 6.3 is the only reason that
is visible rather than buried in a composite.

### 3.2 NOT PURSUED

G2 found MCL1 **basal-specific** while BCL2L1 is subtype-broad, and BrighTNess
is TNBC-only. An MCL1 signal appearing in the basal-enriched cohort is therefore
unsurprising, which makes the pattern *less* interesting rather than more.

Chasing it would be the fifth post-hoc hypothesis plan section 2 forbids. It
joins `Glycine metabolism` (Block C) and `Folate and 1-C` (Block B) on the
reported-and-not-pursued list.

## 4. The other failures

**F1-c, the two-way.** Meta FE +0.079 (-0.036 to +0.193), p 0.176; RE +0.121,
p 0.290, I2 71%. Predicted positive with a CI excluding zero. Null.

**F1-f, the frozen `STATE` contrast.** Level 4 versus level 3, direction fixed
by plan section 7.5 as level 3 having *higher* pCR:

```
GSE194040   +0.379   OR 1.46   p 0.146      cells 145 / 159
GSE164458   -0.016   OR 0.98   p 0.960      cells  74 /  93
GSE25066    +0.259   OR 1.30   p 0.483      cells  91 /  67
```

Two of three run **opposite** to the declared direction, and the cells were not
degenerate - the contrast had data to work with and went the other way.

**F1-g, the informative failure.** The OXPHOS slope within MYC-high, evaluated
at `MYC = +1 SD` in the primary cohort, is **-0.132 (-0.362 to +0.097),
p 0.259.** Criterion 4 requires that slope to be negative with a CI excluding
zero *while the three-way is null*. It is not. **This is a null, not an
inversion.**

**F1-d, specificity, passes vacuously.** The primary arm did not fire in the
qualifying cohorts either (GSE194040 +0.039 p 0.743; GSE164458 +0.255 p 0.110),
so "not the assembly factors" passes because nothing fired - exactly as Block
B's E2b and Block D's D-b/D-c did. It is not evidence of specificity. Note also
that in BrighTNess the control arm sits at 63% of the primary arm's magnitude in
the same direction; that is not a clean separation.

## 5. What the exposure does on its own

Reported because a three-way built on a score with no marginal signal is not
interpretable:

```
              MYC                    OXPHOS                BUFFER_c
GSE194040   +0.238  p 0.0011      -0.074  p 0.321      +0.107  p 0.270
GSE164458   +0.238  p 0.026       -0.066  p 0.538      +0.117  p 0.346
GSE25066    +0.708  p 4.8e-08     -0.120  p 0.367      -0.079  p 0.649
```

**MYC predicts pCR in all three cohorts; OXPHOS does nothing in any of them.**
The MYC association is not news - proliferative tumours respond to chemotherapy,
and MYC-high is proliferative. The exposure H4 is actually about has no marginal
effect anywhere.

The treatment stratification (descriptive, bounded by C3) is null on both sides:
control +0.003 (n = 179, 31 events), experimental +0.159 (n = 564, 175 events).

## 6. CORRECTION - F1-b is a rule I wrote badly

**Recorded beside the registered reading, not folded into it.**

F1-b is scored **pass** because GSE194040's three-way is -0.008, which is
negative. But the rule as I wrote it tests only the **sign of a point estimate,
with no magnitude requirement**, so a coefficient of -0.008 with a CI of
-0.239 to +0.222 satisfies it.

**A coefficient indistinguishable from zero cannot corroborate a direction.**
F1-b carries no information and should be read as uninformative rather than as
a pass. The registered table stands as registered; this is the correction beside
it, in the same pattern as the BCL2L11 declaration's section 5.1.

F1-h passes for a related reason - the composite does not fire in the meta, so
"the composite is not doing something neither limb does" is trivially satisfied.
It was declared as a reporting requirement rather than a prediction, so this is
not a defect, but it is not a pass in any meaningful sense either.

## 7. Where the falsification criteria stand

| Criterion | Status |
|---|---|
| `MYC:OXPHOS` on PRIME null in the cohort and every stratum | **MET** |
| No `BUFFER` enrichment in the quadrant at CNV or expression level | **MET** |
| FOXO3 regulon shows no relationship to OXPHOS in any stratum | **MET** |
| H4 fails in the **informative** direction | **NOT MET** - it failed as a null |

**Three of four, unchanged.** The fourth is unmet only because criterion 4
demands an active *inversion*, which is a stronger claim than a null. That is a
technicality cutting the wrong way: it records that the informative result did
not arrive, not that anything survived.

**All four hypotheses have now been tested and none is supported.** H1 falsified
on both clauses; H2 and H3 flat in every stratum; the escape reading unsupported
and its reversal not mitochondrial; Block G failed at adequate power; H4 null,
with its point estimate opposed to the prediction.

## 8. What this result is NOT, and what it does not license

- **It is single-instrument.** Amendment A1: mitoPPS does not exist in any of
  these cohorts. Blocks C, B and G each required agreement between two
  instruments that demonstrably disagree. **H4 could not clear that bar and must
  not be written as though it had.**
- **A null is not falsification.** D5 section 6.3, pre-registered: H4 is powered
  for a large conditional effect only - 86% at rr = 2.0 pooled, under 60% in any
  single cohort. This must be reported as "not powered to exclude a modest
  effect".
- **Residual confounding cannot be excluded** in the primary cohort. Amendment
  A4: stage and grade do not exist in GSE194040, so the plan's F1 covariate set
  was never fittable. Subtype absorbs part of it and not all.
- **It does not license the two-way as a fallback.** Plan Block F1 is explicit:
  dropping BUFFER buys about five percentage points of power and sacrifices the
  entire distinction from Lee et al. 2017. F1-c is reported because it was
  declared, not as a retreat.

## 9. The BIM replication cannot be run in these cohorts as declared

Recorded here because it was the assumed next step and it is blocked.

The Block C note section 9 declaration fixes, among other things: **"Both
co-primary instruments; the claim requires both, as in Block C."** Amendment A1
established that mitoPPS does not exist in any of the three neoadjuvant cohorts.
So the replication **as declared** cannot be run in them - the same fact that
forced A1 for H4 blocks it here.

A1 was an amendment made for H4 only, on feasibility grounds. The BIM
declaration's own instruction is harder: *"If it is not honoured, delete it
rather than amend it."*

The options, none of which is taken here:

1. **Run it where both instruments exist.** SCAN-B is RNA-seq with counts and is
   already on the plan's cohort list; METABRIC is an array and cannot carry
   mitoPPS either.
2. **Delete the commitment**, as its own text provides for.
3. Run it single-instrument in the neoadjuvant cohorts and report it explicitly
   as **not** the declared replication. This is the weakest option and would
   need saying so in every sentence that mentions it.

This is a decision for the author, not a default.
