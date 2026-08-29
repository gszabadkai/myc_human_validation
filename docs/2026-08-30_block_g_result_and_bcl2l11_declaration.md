---
date: 2026-08-30
status: Block G PRIMARY RESULT recorded; BCL2L11 check DECLARED, not yet run
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

## 5. Status

Block G's declared predictions **fail**, at adequate power. The BCL2L11 check is
declared above and not yet run.
