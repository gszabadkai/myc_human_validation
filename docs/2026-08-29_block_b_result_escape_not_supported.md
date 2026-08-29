---
date: 2026-08-29
status: RECORDED - written before Block F and before the F3 stop gate
relates-to:
  - 2026-08-29_escape_reading_declaration.md (the directions, fixed pre-fit)
  - 2026-08-29_block_c_result_H1_not_supported.md
  - 2026-08-28_G2_result.md
  - 2026-08-27_human_validation_plan.md (sections 2, 10, 15)
decides:
  - H1 has now failed on BOTH clauses.
  - The escape reading is not supported, and fails in the REVERSED direction.
  - Three of the four falsification criteria are now met. Only H4 remains.
next-action: the F3 incremental-value stop gate (plan section 12), before any
  Block F code is written.
---

# Block B - the escape reading is not supported

Run 2026-08-29 via `scripts/10_buffering_models.R`, built to the declaration
signed off the same day with every direction fixed before any of these models
existed. n = 938, both co-primary instruments. The analysis frame was asserted to
reproduce Block C's spine coefficient exactly before anything was fitted.

---

## 1. Every declared prediction fails

| Test | Predicted | Result |
|---|---|---|
| **E1** quadrant enriched for `BUFFER` (H1 clause 2) | positive | **FAIL** |
| **E2a** `BUFFER` predicts OXPHOS subunits, beats matched null | positive | **FAIL** |
| **E2b** and not the assembly factors (specificity) | null / weaker | pass, vacuously |
| **E3** MCL1 expression interaction | positive | **FAIL** |
| **E4** BCL-XL effect attenuated where amplified | positive | **FAIL** (underpowered, as pre-stated) |

E2b "passes" only because nothing fires consistently. It is not evidence of
specificity and must not be reported as though it were.

## 2. E1 - H1 clause 2 fails, and the instruments disagree in sign

```
                                   GSVA                    mitoPPS
categorical (as written)   -0.321  OR 0.73  p 0.308   +0.225  OR 1.25  p 0.463
continuous companion       -0.022  OR 0.98  p 0.782   +0.010  OR 1.01  p 0.906
quadrant vs rest           -0.259  OR 0.77  p 0.138   +0.042  OR 1.04  p 0.814
                           quadrant n = 300           quadrant n = 282
```

No enrichment on either instrument, on any of the three forms, and the two
instruments point in opposite directions. The continuous companions - which D5
showed are the *more* powerful form - are indistinguishable from zero.

**H1 is now falsified on both of its clauses.** Clause 1 in Block C, clause 2
here.

**A precision that matters for the text.** G2's finding stands: MYC amplification
and MCL1 / BCL2L1 amplification genuinely co-occur (adj OR 1.69 and 1.77). E1
shows that co-occurrence is **not conditional on OXPHOS**. Buffering co-occurrence
is a MYC phenomenon, not a MYC-times-OXPHOS phenomenon. Those are different
claims and only the first survives.

## 3. E2 - the reversal, and why it is not a finding

The prediction was positive. **34 of 36 observed coefficients are negative**
(18 of 18 on GSVA, 16 of 18 on mitoPPS). Buffer amplification is associated with
*lower* mitochondrial scores, not higher - a falsification in the informative
direction rather than a null.

But the matched null shows the reversal is not mitochondrial:

```
                          b_obs    null median   percentile   p_emp
OXPHOS subunits, GSVA    -0.138      -0.081        24.7       0.494
OXPHOS subunits, mitoPPS -0.137      +0.041         1.4       0.028
GSVA null medians across all arms: median -0.043
```

**On GSVA a random expression-matched gene set also goes negative.** The negative
shift is a property of any comparably expressed set of genes in these samples -
a global consequence of the copy-number event, not a mitochondrial constraint.
OXPHOS subunits sits at **percentile 24.7**, indistinguishable from random. It
beats its null on mitoPPS but not on GSVA, so the co-primary rule - fixed before
any model was fitted - rejects it.

**This is the clearest demonstration so far of why the 2,000-set null was worth
building.** Without it, "buffer amplification lowers OXPHOS, p = 0.043" is a
publishable-looking sentence. With it, the effect is smaller than what random
genes give.

### 3.1 One arm survives the co-primary rule, and it is not ours

Of 34 arm-by-instrument null tests, 8 reach `p_emp < 0.05` against 1.7 expected -
so there is more than chance here. But only **one arm holds on both
instruments**:

```
Folate and 1-C    GSVA    -0.269   percentile 0.10   p_emp 0.002
                  mitoPPS -0.191   percentile 0.05   p_emp 0.001
```

Buffer amplification is associated with lower one-carbon metabolism, on both
instruments, beating the matched null, marginally surviving Bonferroni over the
34 tests.

**It is reported and it is not pursued.** It is a comparator arm, it bears on no
hypothesis, and - the distinction that matters - **it has no prior**. The BIM
result in Block C earned its pre-registered replication because the mouse arm had
named `Bcl2l11` in advance from independent westerns. This has nothing behind it.
Declaring a second replication off a comparator arm with no prior would be
fishing with extra steps, and the plan's prohibition on a fifth post-hoc
hypothesis exists for exactly this moment.

## 4. E3 - the expression-level buffer prediction fails consistently

```
log2(MCL1)    +0.001 (p 0.943)   +0.017 (p 0.329)     null
log2(BCL2)    -0.073 (p 0.052)   -0.070 (p 0.072)     negative, marginal
log2(BCL2L1)  -0.034 (p 0.030)   -0.043 (p 0.009)     negative  [Block C, restated]
```

All three anti-apoptotic guardians go down or nowhere with `MYC x OXPHOS`. There
is no buffering at the expression level in any direction that helps the escape
reading. Combined with E1, **there is no buffering at the DNA level conditional
on OXPHOS either.**

`log2(BCL2L1)` is restated from Block C's saved object rather than refitted,
because it had already been run as `limb: log2_BCL2L1` and no prediction was
declared for it.

## 5. E4 - null, and uninformative by prior agreement

`MYC:OXPHOS:BUFFER` gives -0.017 (p 0.570) and -0.011 (p 0.713). The declaration
fixed in advance that a three-way term on a 512/426 split is underpowered by
construction and that a null here is uninformative rather than evidence against.
That still holds and is not revisited now that the number is known.

## 6. Where the falsification criteria stand

Plan section 2 requires all four:

| Criterion | Status |
|---|---|
| `MYC:OXPHOS` on PRIME null in the full cohort and every stratum | **MET** (Block C) |
| No `BUFFER` enrichment in the MYC-high / OXPHOS-high quadrant at CNV **or** expression level | **NOW MET** (E1, E3) |
| FOXO3 regulon shows no relationship to OXPHOS in any stratum | **MET** (Block C) |
| H4 fails in the informative direction | **UNTESTED** |

**Three of four. Only H4 remains.**

## 7. What survives the human arm

One finding, and it is the same one Block C reported:

> `MYC x OXPHOS` is associated with **lower BCL-XL and higher BIM**, on both
> co-primary instruments, with PUMA unmoved. `BCL2` moves the same way as
> `BCL2L1`, marginally. `MCL1` does not move.

So in human tumours the MYC-high / OXPHOS-high state carries **more** BH3-only
signal and **less** guardian transcript, with **no compensating amplification**.
Nothing else in Block B or Block C survives both instruments and its own null.

## 8. What this implies for H4, stated before Block F is built

H4 predicts that MYC-high / OXPHOS-high tumours are chemo-**sensitive** unless
buffered. Two consequences follow from the results above, and they pull in
opposite directions:

- **The core prediction gains a mechanism.** If the state carries more BIM and
  less BCL-XL, it should be more primed and therefore more chemosensitive. That
  is H4's direction, now supported by something observed in this cohort rather
  than assumed from the mouse.
- **The conditionality loses its rationale.** `STATE` levels 3 and 4 differ only
  by `BUFFER`, and E1 found no excess buffering in that quadrant to condition on.

Those are not contradictory. E1 asked a **selection** question - is the quadrant
enriched for buffer - and H4 asks a **consequence** question - does buffer change
treatment response. A quadrant can be unenriched for something that still matters
when present. So the level-3-vs-4 contrast remains testable; what it has lost is
its motivation, and that must be said in the text rather than glossed.

The **continuous** three-way primary fixed by D5 is unaffected.

## 9. What happens next

1. **The F3 incremental-value stop gate**, per plan section 12, before any Block F
   code is written. If `STATE` adds nothing over `MB1_forkscale`, there is no
   Panel c and the arm is finished.
2. If it passes, Block F with the reading in section 8 stated up front.
3. The BIM replication declared in the Block C note section 9 is untouched by any
   of this and remains the one thing this arm has to test in independent cohorts.
