---
date: 2026-08-30
status: Block G PRIMARY RESULT recorded; BCL2L11 check RUN (see section 5)
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 2, 10 Block G, 12)
  - 2026-08-29_block_c_result_H1_not_supported.md (section 0, section 9)
  - data/depmap/README.md
decides:
  - G-a and G-b FAIL, and not for want of power. The plan's "a negative result
    here weakens H1 considerably" is cashed in.
  - The BCL2L11 pan-cancer check is declared below, with directions and read
    rules fixed BEFORE it is fitted.
  - How a universal result will be read is fixed here too, not afterwards.
next-action: add section 4's models to script 14, re-source, then read section 5
---

# Block G result, and the declared BCL2L11 check

Script 14 run 2026-08-30 against DepMap Public 26Q1 and PRISM Repurposing 24Q2.
Results in `results/depmap_dependency.rds`.

## 1. The declared predictions fail

| Test | Predicted | Result |
|---|---|---|
| G-a MCL1 dependency, `MYC:OX` negative, both instruments | negative | **FAIL** |
| G-b BCL2L1 dependency, `MYC:OX` negative, both instruments | negative | **FAIL** |
| G-c and not the assembly factors | null / weaker | pass |
| G-d and not BCL2 | null / weaker | pass |

Breast, n = 51: MCL1 +0.057 (p 0.59) GSVA and +0.143 (p 0.35) mitoPPS; BCL2L1
+0.097 (p 0.51) and +0.042 (p 0.85). If anything the estimates are positive,
which is the wrong sign.

## 2. This is not a power excuse

The lineage-adjusted pan-cancer fit is n = 1,130 with confidence intervals of
about +/- 0.04:

| Gene | Instrument | Estimate | 95% CI | p |
|---|---|---|---|---|
| MCL1 | GSVA | -0.009 | -0.040, 0.022 | 0.57 |
| BCL2L1 | GSVA | -0.029 | -0.069, 0.012 | 0.17 |
| MCL1 | mitoPPS | +0.009 | -0.023, 0.040 | 0.59 |
| BCL2L1 | mitoPPS | -0.021 | -0.062, 0.020 | 0.31 |

At that precision an effect of any consequence would have been visible. Plan
section 10 says a negative result here "weakens H1 considerably". That is now
cashed in, and it cannot be attributed to n.

The common-essential floor behaves (RPL3 null, -0.035 / -0.093). `POLR2A` is not
screened in 26Q1, so the floor rests on RPL3 alone. G3 drug sensitivity is
entirely null at n = 26-33, which is too small to mean anything either way; note
also that the two selective BCL-XL inhibitors are absent from PRISM
(`data/depmap/README.md`).

## 3. The one coefficient that is not null

`BCL2L11` (BIM), breast, n = 51:

- mitoPPS **+0.107, CI 0.054 to 0.159, p 0.00018**
- GSVA **+0.040, CI -0.001 to 0.080, p 0.055** - same direction, marginal

Chronos gene effect is more negative for more essential, so a POSITIVE
interaction means that in MYC-high / OXPHOS-high lines, losing BIM helps the
cell more - BIM is more of a liability there. That is the same direction as the
TCGA finding that `MYC x OXPHOS` associates with HIGHER BIM transcript.

**Four reasons this is held loosely.** It is a CONTROL gene, not a declared
endpoint - no prediction was registered for it. n = 51, where one leverage point
can carry a coefficient. There is no expression-matched null in CCLE, so the arm
panel is a rank ordering rather than a calibrated p-value. And one instrument is
significant while the other is only marginal.

**One reason it is not a fishing result.** BIM is pre-specified in the mouse arm
(`myc_mouse/scripts/44`, `PRESPEC <- c("Bbc3", "Bcl2l11")`), it is the finding
that survived Block C on both instruments, and a replication commitment for it
is already live in the Block C note section 9.

**What it is not.** DepMap gene effect is a DEPENDENCY readout; Block C's is a
TRANSCRIPT readout. Consistency of direction across the two is suggestive, not a
replication of one by the other.

---

# 4. The declared BCL2L11 check

Declared 2026-08-30, **before any of the models below has been fitted**. Nothing
in this section has been computed at the time of writing.

## 4.1 Why it is being run

Not to establish universality. **To find out whether the breast result exists at
all.** At n = 51 a coefficient and a pair of leverage points are not
distinguishable; at n = 1,130 they are. This is the same test that has just
disposed of G-a and G-b, applied to the one coefficient that survived.

## 4.2 The three models

A single pooled fit would conflate "is it real?" with "is it universal?", so
they are separated. All are on Chronos gene effect for `BCL2L11`, with
`MYC * OXPHOS subunits + PROLIF` and both instruments, exactly as the existing
Block G models.

```
P1  pooled pan-cancer, lineage-adjusted     Y ~ MYC * OX + PROLIF + lineage
      -> is the effect real?

P2  pan-cancer EXCLUDING breast             Y ~ MYC * OX + PROLIF + lineage
      -> is it universal? and this is genuinely OUT OF SAMPLE: it never sees
         the 51 lines that produced the observation.

P3  breast vs rest                          Y ~ MYC * OX * is_breast + PROLIF + lineage
      -> is it enriched in breast? UNDERPOWERED by construction (51 against
         ~1,079), so a null here is uninformative and is reported as such.
```

`OXPHOS assembly factors` is fitted alongside as the specificity control in P1
and P2, as everywhere else in this repo.

## 4.3 Direction, fixed now

**POSITIVE `MYC:OX`.** That is the direction observed in breast on both
instruments, and the direction implied by the TCGA transcript finding. A
negative result is not a pass and will not be reported as one.

Both instruments must agree in sign, per the standing two-instrument rule.

## 4.4 How each outcome will be read - fixed before the fit

This table exists because the "universal versus breast-specific" question is
exactly the kind that gets rationalised after the answer is known.

| P1 | P2 | Reading | Where it goes |
|---|---|---|---|
| positive | positive | Real, and **not breast-specific**. A property of the MYC/OXPHOS state, not of mammary context. | One sentence in the main text plus Extended Data. **Not a main panel.** |
| positive | null | Consistent with breast enrichment, but P2 null is weak evidence for that on its own - check P3 before saying so. | Extended Data, stated with the caveat |
| null | - | The breast coefficient was noise at n = 51. **Say so and drop it.** | Reported as a negative |
| positive, P3 also positive | - | Enriched in breast. The most interesting outcome and the only one that speaks to the mammary axis directly. | Candidate for Extended Data with the BIM story |

**On "universal is off topic".** Raised by GS before the fit, and worth
answering rather than deflecting. Breast-specificity was never this project's
claim: the manuscript's claim is that mitochondria integrate oncogenic and
metabolic transcriptional programmes - a cell-biological mechanism, tested in a
mammary model. A mechanism that holds across lineages is *on* topic for that
framing and is arguably stronger evidence for it. What the objection does
establish is a **budget** decision, and it is taken here: a universal result
earns a sentence and an Extended Data panel, **never a main display item**,
because the three main panels are about the breast axis.

## 4.5 What a positive still does not buy

- **It does not rescue H1.** H1 is a conjunction whose first clause failed in
  Block C and whose second failed in Block B. Nothing here changes that.
- **It is not reportable without a CCLE expression-matched null.** Script 07's
  nulls are TCGA-specific and do not transfer. Building one in CCLE is a named
  follow-up, not an optional extra, and it is required before any positive
  OXPHOS-arm claim from Block G is written down.
- **It is not the declared BIM replication.** That commitment (Block C note
  section 9) is about transcript in an independent bulk cohort. This is
  dependency in cell lines. Honour that separately or delete it; do not let this
  stand in for it.

## 5. RESULT of the declared check, run 2026-08-30

| Model | Instrument | Estimate | 95% CI | p |
|---|---|---|---|---|
| breast (n=51) | mitoPPS | **+0.107** | 0.054, 0.159 | 0.00018 |
| breast (n=51) | GSVA | +0.040 | -0.001, 0.080 | 0.055 |
| **P1** pooled (n=1,127) | GSVA | +0.00037 | -0.0078, 0.0085 | 0.93 |
| **P1** pooled | mitoPPS | -0.00003 | -0.0083, 0.0082 | 0.99 |
| **P2** excl breast (n=1,076) | GSVA | +0.00055 | -0.0078, 0.0089 | 0.90 |
| **P2** excl breast | mitoPPS | -0.00086 | -0.0093, 0.0076 | 0.84 |
| **P3** breast vs rest (n=1,127) | GSVA | +0.021 | -0.024, 0.067 | 0.36 |
| **P3** breast vs rest | mitoPPS | **+0.089** | 0.035, 0.143 | 0.0012 |

By the rules in section 4.4 the machine verdict is **"P1 null - the breast
BCL2L11 coefficient was noise at n = 51. Report as a negative and drop it."**
That is the registered outcome and it stands as registered.

### 5.1 CORRECTION - section 4.4's P1 row was a bad inference, and it was mine

Written 2026-08-30, after the fit. **This does not change the registered verdict
above.** It records that one line of the declaration was logically invalid, which
is a different thing from disliking its answer.

Section 4.4 says "P1 null -> the breast coefficient was noise at n = 51". **P1
cannot show that.** P1 pools 51 breast lines into 1,127. If the effect were real
and confined to breast, the pooled estimate would be diluted to roughly
51/1127 x 0.107 = 0.005 - comfortably inside the observed CI of +/- 0.008. So
**P1 was guaranteed to be null under both hypotheses** and had no discriminating
power between them. Section 4.1's claim that this is "the same test that
disposed of G-a and G-b" was wrong: G-a and G-b were predicted to hold
pan-cancer, so a pooled null falsified them. This one was not.

### 5.2 What the check DID establish

**The effect is not universal.** P2 fits 1,076 non-breast lines and returns zero
on both instruments with CIs of about +/- 0.008 - roughly a tenth of the breast
estimate. That is a genuine, well-powered negative, and it directly answers the
concern GS raised before the fit: there is no cross-lineage BIM effect here to be
off topic about.

### 5.3 What remains unresolved, and why P3 does not settle it

Whether the breast +0.107 is real is **still open**, and nothing run here
resolves it.

P3 is not independent evidence. It uses the same 51 lines that produced the
observation, and given breast at +0.107 and the rest at ~0, a significant
breast-vs-rest contrast is close to arithmetically implied. It restates the
input rather than testing it. P3 also **fails the two-instrument rule** - GSVA
+0.021, p 0.36 - and so does the breast result it is built on, where GSVA is
only p 0.055.

Two further reasons for restraint. The 14-row table above contains a nominally
significant hit on the **negative control arm** (P2, assembly factors, mitoPPS:
-0.0087, p 0.046). That is the multiplicity floor made visible. And mitoPPS is a
composition measure, so a lineage with distinctive mitochondrial composition can
move it for reasons that have nothing to do with the MYC/OXPHOS state - breast
being one lineage among ~30 here.

**No claim of breast-specificity is made or licensed.** Concluding "not
universal, therefore breast-specific" from a null in the pooled fit would be the
post-hoc rescue this declaration exists to prevent.

### 5.4 What would actually settle it

1. **The CCLE expression-matched null, in breast.** Already a standing
   requirement (section 4.5). It asks the only question that matters at n = 51:
   does +0.107 beat a random expression-matched gene set in these same lines?
2. **A leverage check.** The sandbox plots BCL2L11 against mitoPPS OXPHOS for
   the 51 lines; at that n a couple of points can carry the coefficient.
3. **An independent breast dataset**, which is what the Block C note section 9
   commitment is for - and that is transcript, not dependency, so it remains a
   separate obligation.

## 6. Status

Block G's declared predictions **fail**, at adequate power. The BCL2L11 check
returns a registered verdict of "noise, drop it"; that verdict's stated
reasoning is corrected in 5.1 without changing the verdict. The effect is
**established as not universal** and remains **unresolved in breast**.
