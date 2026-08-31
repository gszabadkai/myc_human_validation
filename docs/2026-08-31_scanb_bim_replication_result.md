---
date: 2026-08-31
status: RESULT - the BIM replication FAILED. The declaration is honoured and the
        commitment is discharged.
relates-to:
  - 2026-08-29_block_c_result_H1_not_supported.md (section 9, the declaration)
  - 2026-08-31_scanb_bim_replication_declaration.md (every parameter, pre-data)
  - 2026-08-31_block_F1_result_H4_not_supported.md (section 9, the three options)
decides:
  - BCL2L11 does NOT replicate in SCAN-B. It is reported as a TCGA-specific
    observation and dropped, exactly as the declaration provides.
  - No further variants are tried. The arm has no remaining analysis.
next-action: the write-up. Plan section 15.
---

# The BIM replication failed, and it failed as a reversal

**`BCL2L11` does not replicate in SCAN-B.** The interaction was predicted
POSITIVE. In 3,143 patients it is **negative on both co-primary instruments,
nominally significant on both, and stable across all four specifications.**

```
                             GSVA                        mitoPPS
TCGA   (n =   938)   +0.051 (+0.015, +0.086) p 0.005   +0.068 (+0.031, +0.105) p 3.0e-4
SCAN-B (n = 3,143)   -0.036 (-0.066, -0.006) p 0.019   -0.033 (-0.065, -0.002) p 0.037
```

This is the informative failure rather than the uninformative one. The
declaration's failure condition applies verbatim and is discharged here:

> *"if `BCL2L11` does not replicate in the independent cohorts with the
> direction above, it is reported as a TCGA-specific observation and dropped.
> No further variants are tried."*

**It is dropped.** Not another cohort, not another covariate set, not another
estimator, not a subtype stratum.

---

## 1. The gate ruled out the one explanation that would have been ours

The SCAN-B covariate set is four terms short of the Block C spine, and without
the calibration a reader could reasonably say the reduction caused the failure.
It did not, and this was established **before SCAN-B was scored**.

Refitting the three TCGA limbs on the same 938 patients, same frozen z-scaling,
with covariates reduced to `PROLIF_DISJOINT + PAM50` - the SCAN-B mapping
exactly:

```
                    full spine                  reduced (SCAN-B mapping)
BCL2L11 GSVA     +0.0509 (0.0154, 0.0864)     +0.0461 (0.0107, 0.0815) p 0.011
BCL2L11 mitoPPS  +0.0681 (0.0313, 0.1050)     +0.0599 (0.0235, 0.0962) p 0.0013
BCL2L1  GSVA     -0.0342 (-0.0651, -0.0033)   -0.0362 (-0.0673, -0.0051)
BCL2L1  mitoPPS  -0.0434 (-0.0760, -0.0108)   -0.0417 (-0.0742, -0.0091)
```

Dropping purity, leukocyte fraction, TP53 status and plate costs about **10% of
the BCL2L11 estimate and none of its sign or significance**, on both
instruments, and the full-spine point estimate sits inside the reduced-set 95%
CI in every case. The gate passed on the condition fixed in advance.

**So the sign reversal in SCAN-B is not the covariate mapping.** That is the
single most useful thing the gate bought, and it could only be bought by fixing
the pass condition before either number existed.

The rebuilt TCGA frame reproduced script 09's saved coefficients to **exactly
zero difference** on all six limb x instrument fits, so the calibration measured
the covariate set and not the rebuild.

## 2. The result

`log2(BCL2L11) ~ M_a * OXPHOS_subunits + PROLIF_DISJOINT + PAM50`, SCAN-B
GSE202203, n = 3,143 of 3,207 (PAM50-complete).

```
instrument  spec     n    estimate    95% CI               p
GSVA        S1     3143   -0.0362   (-0.0665, -0.0059)   0.019   <- primary
GSVA        S2     3143   -0.0365   (-0.0669, -0.0062)   0.018
GSVA        S3     3143   -0.0362   (-0.0666, -0.0059)   0.019
GSVA        C-alt  3090   -0.0367   (-0.0673, -0.0060)   0.019
mitoPPS     S1     3143   -0.0334   (-0.0648, -0.0021)   0.036   <- primary
mitoPPS     S2     3143   -0.0345   (-0.0659, -0.0032)   0.031
mitoPPS     S3     3143   -0.0335   (-0.0648, -0.0022)   0.036
mitoPPS     C-alt  3090   -0.0316   (-0.0633, +0.0001)   0.051
```

The four specifications agree to the third decimal. Whatever this is, it is not
a covariate artefact either: dropping the proliferation term (S2), swapping in
the standard one (S3), and adding age and grade (C-alt) all leave it where it
is.

## 3. The null says the reversal is specific, not a scale artefact

The declaration required the arm to beat its expression-matched null on this
endpoint. It does - **in the wrong direction**, which is exactly why condition 4
fails.

```
instrument   b_obs      null mean   null sd   percentile   draws below   p_emp
GSVA        -0.0362     +0.0034     0.0159      0.50        10 / 2000    0.010
mitoPPS     -0.0334     +0.0090     0.0144      0.05         1 / 2000    0.001
```

Only 10 of 2,000 expression-matched random gene sets of the same size and
abundance profile produce a coefficient this negative on GSVA, and 1 of 2,000 on
mitoPPS. **The negative association is more extreme than 99.5% and 99.95% of
matched sets.** It is a real, specific, inverted signal - not noise and not an
artefact of set size or expression level.

That matters for how the failure is described. This is not "we could not detect
the effect in a second cohort". It is "the effect is there, at n = 3,143, with
the opposite sign."

## 4. The companions

Reported with it, as Block C reported them.

```
endpoint     predicted   GSVA                      mitoPPS
BCL2L11      positive    -0.036  p 0.019           -0.033  p 0.037
BCL2L1       negative    -0.011  p 0.14            -0.026  p 6.4e-4
BBC3         null        -0.011  p 0.52            -0.068  p 7.5e-5
```

**`BCL2L1` does not corroborate.** The direction is right on both instruments,
but it is significant on mitoPPS only, and the co-primary rule makes a
one-instrument result not a result. It is reported as such.

**`BBC3`, the control, fires on mitoPPS** - and it is the LARGEST of the three
effects on that instrument. That is a specificity failure. On mitoPPS in SCAN-B
all three BCL2-family genes fall together with `MYC x OXPHOS` composition
(-0.068, -0.033, -0.026), which is the signature of a general shift across the
family rather than a priming-specific effect. On GSVA the pattern is different:
`BCL2L11` moves and the other two do not.

The two instruments therefore fail the replication for partly different reasons,
and neither rescues the other. **They agree that the sign is negative**, which
is the only thing condition 1 asks.

For completeness, since the co-primary rule depends on the instruments being
different measurements: `rho(GSVA, mitoPPS)` on `OXPHOS subunits` is **0.825 in
SCAN-B against 0.88 in TCGA**. They are no more redundant here than there, so
"both instruments" is doing the same work in both cohorts.

## 5. The four conditions

```
1. interaction POSITIVE on both instruments               NOT MET
2. p < 0.05 on both instruments                           met
3. both co-primary instruments fitted at S1 (STRUCTURAL)  met
4. beats its matched null on this endpoint                NOT MET
```

**Condition 2 is direction-blind and its "met" is misleading. Recorded rather
than quietly fixed.** It tests a two-sided p, so a result significant in the
*wrong* direction passes it. Condition 1 catches that and the verdict is
correct, but a reader scanning this table sees three of four met and could take
that for a near-miss. It is not: the two evidential conditions that carry
direction both fail.

This is the same error class as **F1-b** (F1 note section 6) - a criterion that
tests less than its label implies. Two instances now, both mine, both found
after the fact. The general lesson for any future declaration in this project:
**a directional prediction needs a directional criterion**, and writing
`p < 0.05` beside a predicted sign does not make one.

Condition 3 is structural - it asserts the two fits exist and can fail only
through a coding error. Script 17 stops rather than returning FALSE for it, and
it is flagged as non-evidential in the output.

## 6. What this result is NOT

Four fences, because a reversal is more tempting than a null.

**It is not a fifth hypothesis.** "`MYC x OXPHOS` lowers BIM in SCAN-B" is a new
claim built backwards from a failed replication, and plan section 2 forbids
exactly this move. The finding being reported is that **the TCGA result did not
replicate**. The negative coefficient is the evidence for that, not a result in
its own right.

**It does not close falsification criterion 4.** Criterion 4 names H4 and
chemoresistance: *"H4 fails in the informative direction: OXPHOS-high predicts
chemoresistance regardless of buffering status."* That is about pCR, and H4
failed as a null (F1 note). `BCL2L11` is a Block C limb, a different endpoint,
and an inversion here is not the inversion criterion 4 asks for. **Three of four
criteria remain met; the fourth remains unmet.** The tally is unchanged.

**It is not explained, and the explanation is not being chased.** SCAN-B is a
population-based screening-era cohort, roughly 75% ER-positive and largely
early-stage; TCGA-BRCA is a referral series with more advanced disease. Those
are real differences and any of them could plausibly matter. **None of them is
tested here.** Testing them is precisely the "further variants" the declaration
forbids, and a subgroup that restored the sign would be the least trustworthy
number in the paper. A sign reversal between two large cohorts of the same
disease on the same platform type is a limitation of the original finding, and
it is reported as one.

**It does not reopen SCAN-B.** GSE202203 carries survival, relapse-free
interval and treatment flags; none reached `results/scanb_pheno.rds` and script
17 re-asserted that before fitting. F2, F3 and forkscale in SCAN-B remain
separate decisions with separate notes.

## 7. Where the arm stands now

Every hypothesis is unsupported and the one surviving finding has been tested
independently and dropped.

- **H1 falsified on both clauses**, H2 and H3 not supported, Block D null,
  Block G failed at adequate power, **H4 not supported**.
- **The BIM result is now TCGA-specific and dropped.** It was the single finding
  that survived on both instruments, it carried a pre-specified mouse prior, and
  it had the only pre-registered replication in the arm. That replication ran,
  in the cohort the declaration required, on both instruments, and it failed.
- **Three of four falsification criteria met**, unchanged.
- **BCL-XL** remains the one direction that is at least consistent across
  cohorts - down in TCGA on both instruments, down in SCAN-B on both, though
  significant on only one there. It is not claimed. It was never a hypothesis
  and it does not become one now.

**This is a stronger negative than the arm had yesterday, not a weaker one.**
Yesterday the write-up would have had to say that its one positive was
unreplicated. It can now say that the replication was pre-registered, run in an
independent cohort of 3,143 patients on both required instruments, and reported
when it failed. The dated record shows the declaration was written on
2026-08-29, before any of these cohorts had been opened, and the cohort choice
and every parameter were fixed on 2026-08-31 before a single expression value
was read.

## 8. What is left

**No analysis.** The arm is complete.

Outstanding but not gating, and each a separate decision: F3 proper (METABRIC
identifier file), F2, the CCLE expression-matched null, Block D2, D6.

Next is the write-up, on plan section 15's framing, which the evidence now
supports completely.
