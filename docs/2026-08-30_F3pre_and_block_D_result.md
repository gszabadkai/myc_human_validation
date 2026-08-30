---
date: 2026-08-30
status: RECORDED - result note, read against rules fixed before the fit
relates-to:
  - 2026-08-29_G3_result_forkscale_availability.md (section 4b, the F3-pre read rules)
  - 2026-08-29_block_c_result_H1_not_supported.md (criterion 1, scored on an incomplete list)
  - 2026-08-29_block_b_result_escape_not_supported.md (criterion 2)
  - 2026-08-27_human_validation_plan.md (sections 2, 9, 10, 12)
decides:
  - F3-pre verdict is INTERMEDIATE. The Block F stop gate does NOT close.
  - Block D is null on both instruments. D-a FALSE; the matched-null gate stayed shut.
  - Falsification criterion 1 is now discharged on its COMPLETE stratum list,
    forkscale included. Still three of four. Only H4 untested.
next-action: script 11, to freeze STATE before any H4 outcome data is touched.
---

# F3-pre and Block D - the fork axis is neither a duplicate nor a modifier

Run 2026-08-30 via `scripts/15_forkscale_replication.R`. Results in
`results/forkscale_models.rds`. n = 885 of Block C's 938 carry forkscale
(53 lost); the Block D fits are n = 844, because `er_call` has 46 NAs inside
the 938 and D8 is complete-case. Both numbers were pre-stated.

The analysis frame was asserted to reproduce Block C's spine coefficient
exactly before anything was fitted, and the forkscale definition
(`pc1 / index`) was verified against the companion paper's own assembled
frame - n = 849, max abs diff 0 - before it was used.

---

## 1. F3-pre - INTERMEDIATE

Read mechanically against the table in the G3 note section 4b, on the declared
decision variable (plain `MB1_forkscale`, Spearman, `OXPHOS subunits`, both
instruments):

```
                              rho      threshold
GSVA                        0.529      REDUNDANT >= 0.70 on both
mitoPPS                     0.418      INDEPENDENT <= 0.30 on both
                                       -> INTERMEDIATE
```

Above the independence line on both instruments, below the redundancy line on
both. **The stop gate does not close.** Block F is not killed on redundancy
grounds, and the plan's incremental-value requirement (section 12: nested Cox,
likelihood ratio, C-index delta) stands exactly as written - neither relaxed
nor tightened by this number.

Per the declared rules this verdict is reported and **no action is taken on
it here**.

## 2. Block D - null, and the gate correctly stayed shut

The pre-registered `MYC * OXPHOS * forkscale` three-way on PRIME, MB2 primary:

```
                                   estimate     CI                p
MB2 three-way, subunits, GSVA        0.0181  -0.140 .. 0.176    0.822
MB2 three-way, subunits, mitoPPS     0.0145  -0.118 .. 0.147    0.830
```

**D-a is FALSE.** The matched null is gated on the primary three-way being
significant on both instruments; the gate stayed shut and the null was
correctly not run.

Across the 16 pre-registered coefficients (8 three-way, 8 within-ER),
**not one clears the two-instrument rule.** Two rows are nominally notable and
neither survives:

```
MB2 within ER-Negative, mitoPPS   +1.87   -0.038 .. 3.79   p 0.055   n=195
MB3 within ER-Positive, mitoPPS   -0.430  -0.840 .. -0.020 p 0.040   n=649
```

Both are single-instrument, in strata where GSVA is flat (0.726 p 0.26 and
-0.287 p 0.35 respectively), and the second is on the **control** axis. The
ER-negative CI spans nearly four units. Two nominal hits in 16 coefficients is
the expected rate. Neither is a finding.

### D-b, D-c and D-e pass VACUOUSLY

They are scored TRUE because nothing fired anywhere, exactly as E2b was in
Block B. **They are not evidence of specificity and must not be reported as
though they were.**

One qualification in D-c's favour: MB3 is a genuine control axis rather than a
relabelled MB1. The three fork axes correlate

```
              MB1     MB2     MB3
MB1         1.000   0.712   0.308
MB2         0.712   1.000   0.160
MB3         0.308   0.160   1.000
```

so the check prescribed for exactly this - if MB1/MB2/MB3 were near-identical,
"MB3 as a control" would be no control at all - passes. D-c was capable of
being informative. It simply had nothing to discriminate.

## 3. What this does to the falsification criteria

Plan section 2 states criterion 1 as null "in the full cohort **and** in every
pre-specified stratum (TP53, PIK3CA, PAM50, **forkscale**)".

Forkscale was the one entry never fitted. Block C ran TP53, PIK3CA, PI3K and
PAM50; forkscale was unavailable at the time, because G3 was not discharged
until the following day. Both the Block C note (section 10) and the Block B
note (section 6) therefore marked criterion 1 **MET on an incomplete stratum
list**. This note closes that.

Plan sections 9-10 fix the pre-registered operationalisation of the forkscale
stratum as the continuous three-way - not a median split - and D5 established
that the continuous form is the more powerful one. That is what ran, and it is
null on both instruments.

| Criterion | Status |
|---|---|
| `MYC:OXPHOS` on PRIME null in the full cohort and every stratum | **MET** - now on the complete list, forkscale included |
| No `BUFFER` enrichment in the quadrant at CNV or expression level | **MET** (Block B, E1 and E3) |
| FOXO3 regulon shows no relationship to OXPHOS in any stratum | **MET** (Block C) |
| H4 fails in the informative direction | **UNTESTED** |

**Still three of four. Only H4 remains.** Nothing here moved the tally; what
changed is that criterion 1 is now discharged on evidence rather than on an
unexamined list.

**The pre-stated qualifier stands and is not weakened by the result.** n = 844,
continuous modifier, underpowered by construction. The claim is "the
pre-registered test was run and returned null", not "there is no fork-dependent
coupling". No direction was pre-specified for Block D and none is asserted now.

## 4. What the fork axis actually is - three descriptive observations

These are about the companion paper's axis, not about this arm's finding. None
is a hypothesis and none is pursued.

**It is not an OXPHOS-subunit axis.** On GSVA it tracks OXPHOS *assembly
factors* (0.652) harder than the subunits (0.529), and PROLIF_DISJOINT at 0.503
- nearly as strongly as the subunits. The instruments disagree about which arm
it prefers; mitoPPS reverses the order (subunits 0.418, assembly 0.341).
Whatever MB1 forkscale is, it is a broad mitochondrial-plus-proliferation axis.

**Its MYC correlation is estimator-dependent.** Forkscale vs M-a (Felsher
signature, GSVA) is 0.429; vs M-b (CollecTRI regulon, ULM) it is **-0.031,
p 0.35** - despite M-a and M-b agreeing with each other at 0.613
(`2026-08-28_myc_estimator_validation.md`). Forkscale tracks the expression-
*signature* flavour of MYC and is invisible to the regulon-*activity*
estimator. Read as a shared expression-program correlation, not a MYC-activity
one.

**It is flat against the primary endpoint.** Forkscale vs PRIME is -0.012
(p 0.72).

## 5. One trap in the output - do not lean on the incremental-variance table

The script reports `delta_r2` of 0.0075 (GSVA) and 0.0055 (mitoPPS) for
`OXPHOS ~ covariates + fork`, and 0.0085 the other way. That invites the
conclusion that the two axes are nearly conditionally independent once PAM50,
proliferation and purity are in. **Do not draw it from this table.**

The table fits a **linear** model on the **raw** forkscale, and raw forkscale
is the skew-dominated form. Its Pearson r against OXPHOS-GSVA is 0.146 where
Spearman is 0.529; the log form's Pearson is 0.457. The table is therefore
measuring the linear-on-severely-skewed component and understates by a large
factor.

This changes nothing about the verdict - F3-pre is decided on the marginal
Spearman, and INTERMEDIATE is correct. But **if that `delta_r2` is ever cited
in text, recompute it on the log form first.** One line in the loop.

## 6. Not pursued

`D MB2 DESCRIPTIVE log2_BCL2L1`, GSVA, is -0.087 (p 0.056). It is one of eight
descriptive rows, carries no prediction, is excluded from every pass/fail by
construction, and fails the two-instrument rule anyway (mitoPPS p 0.257). **It
is not a BCL-XL fork-dependence finding**, and promoting it would need its own
declaration first, written before the number was seen. That is no longer
possible for this quantity.

The BIM commitment from the Block C note section 9 is **untouched by this
script** and remains live. Block D's `log2_BCL2L11` rows are descriptive and do
not discharge it.

## 7. A limitation to register now, before the H4 data arrives

MCbiclust forkscale was fitted in TCGA and METABRIC. **It does not exist for
I-SPY2-990 (GSE194040), BrighTNess (GSE164458) or GSE25066**, and constructing
it there would be new methodology, not a lookup.

So **forkscale cannot be adjusted for in the H4 cohorts.** If H4 returns a
positive, the question "is this just the fork axis under another name?" has to
be argued from the TCGA rho of 0.53 recorded above - in text, as a stated
limitation - because it cannot be tested in the cohort where the claim would be
made.

Registered here, before any H4 outcome data has been touched, so that it is a
declared limitation rather than a discovered one.

## 8. Deferred, not dropped

- **The forkscale-vs-STATE panel of F3-pre.** Reported as unavailable because
  script 11 has not run. Re-source script 15 after 11 to populate it. It gates
  nothing.
- **Block D2.** Still needs the Pommier developmental sets snapshotted with a
  provenance README and a G1-style overlap audit against the MitoCarta arms and
  PROLIF. Importing a gene set is a documented decision in this repo, not a side
  effect of writing a model.
- **F3 proper.** A survival test, and it belongs in METABRIC. Plan section 3
  forbids TCGA survival and `data/tcga_cdr/README.md` quantifies why (151 OS /
  145 PFI events at ~2.3 years median follow-up). Still blocked on the METABRIC
  sample-identifier file - G3 note section 5 lists the three acceptable sources,
  best first.
