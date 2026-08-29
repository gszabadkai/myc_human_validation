---
date: 2026-08-29
status: RECORDED - primary result. Written before anything downstream was looked at.
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 2, 10, 15)
  - 2026-08-29_script09_build_spec.md (the model space, signed off pre-fit)
  - 2026-08-28_D7_proliferation_covariate.md, 2026-08-28_D8_D9_missing_data_and_plate.md
  - results/block_c_models.rds
decides:
  - H1 is NOT supported in TCGA-BRCA. H2 and H3 are not supported either.
  - A pre-registered replication is DECLARED here (section 9) before any
    replication cohort is touched.
next-action: Block B (script 10), which settles the remaining TCGA clause of the
  falsification criteria. Not a fifth hypothesis.
---

# Block C - H1 is not supported

Run 2026-08-29 via `scripts/09_interaction_models.R`, built to the spec signed
off the same day, before any coefficient existed. 116 coefficients over 58
labels, plus 68,000 expression-matched null fits.

This note is written **before** Block B, before script 17, and before any
replication cohort has been opened, so the ordering is checkable.

---

## 0. AMENDMENT, later the same day

Two things were found after this note was first written. **The original text
below is unaltered**; this section says what it got wrong and what it does not
change.

### (a) The plan mis-read the mouse. `BCL2L11` should never have been a negative control

Plan section 2 justified the endpoint-negative panel with *"The mouse says only
PUMA/BCL-XL reverses. Human should show the same."* Checked against the mouse
source, read-only, on 2026-08-29:

```
myc_mouse/scripts/44_collapse_module_and_ownership.R

    PRESPEC <- c("Bbc3", "Bcl2l11")

    "Bbc3 (PUMA) and Bcl2l11 (BIM) are named in advance because the PGC1a
     perturbation induces both at protein level (author's westerns). Those
     experiments are INDEPENDENT of this RNA-seq, so using them to fix which
     genes and sets are tested is pre-registration, not circularity -- and at
     n=24 it is the strongest available position."
```

`myc_mouse/scripts/42` tests `Bcl2l11:Bcl2l1` in the same priming pair panel as
`Bbc3:Bcl2l1`. `myc_mouse/scripts/23` carries a switch, `foxo_puma_select`,
testing whether the FOXO programme is "PUMA-SELECTIVE" or "broad" - so the mouse
was **actively uncertain** about BIM, not settled against it.

What the mouse does say, and it is the sentence the plan compressed:

```
"Bmf and Bcl2l11 have MORE extreme retentions than PUMA but NON-SIGNIFICANT
 6W effects (p6 0.41 and 0.46), and a retention is a ratio of two noisy
 quantities -- you cannot lose an effect you never had. Among the pairs with a
 real 6W effect, PUMA is the ONLY one that reverses sign."
```

**PUMA was chosen because it was resolvable at n = 24, not because BIM was
excluded.** The plan turned "we could not resolve BIM" into "the mouse excluded
BIM", and that is what demoted `BCL2L11` to a control.

**Consequence.** `BCL2L11` is not a post-hoc finding. It is a **co-candidate
pre-specified by the research programme on independent protein evidence, which
this arm's plan mistakenly relabelled as a negative control.** That is a
materially stronger position than section 8 below granted it - and a materially
weaker one than "H1 confirmed".

### (b) The endpoint panel is showing TWO effects, not four

See the block added at the end of section 6. Every ratio is its own numerator
minus the shared `BCL2L1` denominator, exactly; once decomposed, the panel is
BCL-XL down and BIM up, not four scattered hits.

### What does NOT change

- **H1 as pre-registered names PRIME, and PRIME is null.** That stands, and the
  four decision clauses still all fail.
- **The replication declared in section 9 is unaltered.** Its predictions were
  fixed before this amendment and must not move because the rationale improved.
  A prediction that is edited after the fact is not a prediction.
- The BIM result is still S1-only, still untested across the D7 specifications,
  and still unreplicated.

---

## 1. The headline

```
PRIME ~ M_a * OXPHOS subunits + covariates,  n = 938

  GSVA      0.0204   95% CI  -0.0332 to 0.0740   p = 0.455
  mitoPPS   0.0099   95% CI  -0.0452 to 0.0649   p = 0.726
```

Both co-primary instruments, both null. The `MYC:OXPHOS` interaction on the
pre-specified endpoint does not exist in TCGA-BRCA at n = 938.

## 2. The decision rule, clause by clause

Fixed in the build spec section 10 before any fit, evaluated mechanically by the
script rather than in the writing.

| # | Clause | Result |
|---|---|---|
| 1 | interaction positive, CI excludes 0, BOTH instruments | **FAIL** |
| 2 | `OXPHOS subunits` beats its matched null, BOTH | **FAIL** |
| 3 | paired contrast vs `OXPHOS assembly factors` beats its null, BOTH | **FAIL** |
| 4 | endpoint negatives null, and the `BBC3` limb carries the effect | **FAIL** |

## 3. Why the null makes this stronger than a p-value

```
                          b_obs    null median   percentile   p_emp
OXPHOS subunits, GSVA     0.0204     0.0312         31.6      0.632
OXPHOS subunits, mitoPPS  0.0099    -0.0077         73.0      0.541
```

**A random expression-matched gene set gives a slightly LARGER interaction
coefficient than the respiratory chain does.** The observed effect is not merely
non-significant; its magnitude is unremarkable against matched noise. Without
the null this would be "p = 0.455"; with it, the claim is that OXPHOS subunits
behave like an arbitrary set of comparably expressed genes. That is what the
68,000 fits bought, and it is the difference between a null result and an
uninformative one.

The paired contrast **inverts**:

```
subunits - assembly factors   GSVA    -0.0141   percentile 19.4   p_emp 0.387
                              mitoPPS -0.0163   percentile 23.4   p_emp 0.468
```

Assembly factors carry the larger coefficient. The mouse's sharpest specificity
result - subunits at percentile 0.0, their own assembly factors at 50.2 - runs
the other way in human, though not significantly. The other two paired contrasts
(vs mitoribosome, vs nucleotide metabolism) are equally flat.

Across the whole 17-arm battery, exactly one arm reaches `p_emp < 0.05` on
either instrument: `Glycine metabolism` on mitoPPS (percentile 99.2, p_emp
0.017). With 17 arms on 2 instruments that is within expectation and is **not**
pursued.

## 4. H2 and H3 are not supported either

Every pre-specified stratum is flat, on both instruments:

```
                       GSVA              mitoPPS         n
TP53 mutant           -0.016            -0.049          325     H3
TP53 wildtype          0.025             0.051          613     H3
PIK3CA altered         0.041             0.015          328     H2
PIK3CA intact          0.005             0.001          610     H2
PI3K intact            0.008             0.012          549     H2
PI3K altered           0.025            -0.011          389     H2
LumA                  -0.003            -0.017          478
LumB                   0.056             0.092          195
Basal                  0.004             0.016          159
Her2                  -0.102            -0.177           77
```

No CI excludes zero anywhere in that table. H3's prediction (the coupling
persists in TP53-mutant tumours) has nothing to persist. H2's prediction (the
coupling is restored in the PI3K-intact stratum) does not hold: PI3K-intact is
the flattest stratum of all.

FOXO3 as an endpoint is also null - 0.018/0.030 in the full cohort, 0.015/0.017
in PI3K-intact. The mandatory purity-high refit (n = 270) gives 0.069 (p = 0.078)
and 0.052 (p = 0.146), which at that n is uninformative in the direction the spec
said in advance it would be.

## 5. The result does not depend on any specification choice

The whole point of the ladder. Every excursion lands in the same place:

```
S1 PROLIF_DISJOINT    0.0204 / 0.0099      the spine
S2 no proliferation   0.0246 / 0.0119
S3 PROLIF_STD         0.0206 / 0.0100
M2  n = 1,079         0.0104 / -0.0089
M3  n = 1,095         0.0126 / -0.0076
plate random          0.0218 / 0.0134
ER-adjusted           0.0304 / 0.0233       n = 892
within ER-positive    0.0170 / 0.0070       n = 691
within ER-negative    0.0495 / 0.0816       n = 201
purity > 0.7          0.0591 / 0.0370       n = 270
M-c (8q24 amp)        0.0539 / 0.0212
```

D7, D8 and D9 turn out not to have mattered to the answer. They were still worth
settling in advance - that is exactly the claim a reader cannot check after the
fact.

**One instrument-dependent result, killed by the rule that was written for it.**
M-b gives `0.0505, CI 0.0001 to 0.1009, p = 0.049` on GSVA and `0.0190,
p = 0.424` on mitoPPS. Under the co-primary rule fixed before any model was
fitted, an effect on one instrument alone is instrument-dependent and is not a
positive result. Recorded because it is exactly the number a less disciplined
version of this analysis would have led with.

## 6. What DID move - and it is not what was predicted

The endpoint negatives were supposed to be silent. Three of four are not:

```
                       GSVA                mitoPPS
PRIME                  0.020 (p 0.46)      0.010 (p 0.73)    <- pre-specified
BID/BCL2L1             0.062 (p 0.012)     0.062 (p 0.015)
BAX/BCL2L1             0.049 (p 0.010)     0.035 (p 0.071)
BCL2L11/BCL2L1         0.085 (p 5e-4)      0.112 (p 1.9e-5)
BAK1/BCL2L1            0.024 (p 0.26)      0.006 (p 0.80)
```

**The only endpoint that does not move is the one named in advance.** The
limb-wise fits - pre-specified in build-spec section 4 precisely to catch this -
say why:

```
log2(BBC3)     -0.014 (p 0.60)    -0.034 (p 0.21)    PUMA:   unmoved
log2(BCL2L1)   -0.034 (p 0.030)   -0.043 (p 0.009)   BCL-XL: DOWN, both
log2(BCL2L11)  +0.051 (p 0.005)   +0.068 (p 3e-4)    BIM:    UP, both
```

So `MYC x OXPHOS` does raise apoptotic priming in human tumours - the
anti-apoptotic guardian falls and a BH3-only sensitiser rises - but **PUMA is not
the sensor**. The direction of the mouse model survives; its named effector does
not.

The limbs also close the obvious escape route. PRIME did not fail because it was
a poor instrument for PUMA: the `BBC3` limb is null on its own. Script 08 had
already shown PRIME is ~84% its numerator (rho 0.838 with `log2(BBC3)`, -0.128
with `log2(BCL2L1)`), so "PRIME is flat" and "PUMA is flat" are close to the same
statement, and both are true.

### 6.1 Decomposed, the panel is two effects (added in the amendment)

Every ratio endpoint is exactly its own numerator coefficient minus the shared
`BCL2L1` denominator coefficient - an algebraic identity of OLS on the same
design and the same 938 patients, verified to machine precision, not an
empirical finding. Reading it that way:

```
                    ratio     numerator    -denominator    numerator's share
PRIME (PUMA)       0.0204      -0.0138         0.0342          negative
BID/BCL-XL         0.0622      +0.0279         0.0342             45%
BAX/BCL-XL         0.0489      +0.0146         0.0342             30%
BCL2L11/BCL-XL     0.0851      +0.0509         0.0342             60%
BAK1/BCL-XL        0.0241      -0.0101         0.0342          negative
                                                        (GSVA; mitoPPS agrees)
```

**BCL-XL falling is the shared signal in every one of them, and BIM is the only
BH3-only protein whose own level rises.** BID, BAX and BAK1 "fired" mostly by
borrowing the falling denominator. PUMA moves slightly the *wrong* way, and
PRIME is positive at all only because BCL-XL falls underneath it.

So the concern this note originally recorded - "three of four negatives fired,
which is a specificity failure" - largely dissolves. The panel is not four
scattered effects. It is **two**: BCL-XL down, BIM up. That is a coherent and
specific mechanism, and it is the one the mouse pre-specified from westerns.

It also means the four ratio endpoints carry **no information beyond the six
limbs**. Any future reporting should lead with the limbs; the ratios are a
derived view of them.

## 7. What the BIM result survives, stated with its context

- **Both instruments**, which is the bar this project set in advance.
- **Multiplicity**: 21 of 114 coefficients reach p < 0.05 against 5.7 expected by
  chance. `BCL2L11` clears Bonferroni over all 114 on mitoPPS (1.9e-5 against a
  4.4e-4 threshold) and is marginal on GSVA (5e-4).
- **Circularity**: `BCL2L11` is in none of `OXPHOS subunits`, `OXPHOS assembly
  factors`, `OXPHOS umbrella`, `PROLIF_DISJOINT`, or the M-a estimator. It is
  MitoCarta-annotated (`Apoptosis`) and it **is** a CollecTRI MYC target, so it
  would be circular with M-b - but every fit above uses M-a.
- **Not tested across the D7 specifications.** The ladder varies one dimension at
  a time and endpoint variants run at S1 only. So the BIM result has not been
  shown robust to the proliferation covariate, and must not be described as if
  it had.

## 8. What this is not

> **CORRECTED by section 0(a).** The premise of this section - that `BCL2L11`
> was a negative control - was this arm's error, not the mouse's. `BCL2L11` was
> pre-specified alongside `Bbc3` in `myc_mouse/scripts/44` from independent
> western-blot evidence. The paragraph below is kept as written, and the
> corrected reading follows it.

`BCL2L11` was a **pre-registered negative control that came out positive**. That
makes its result legitimate to report. It does not make it a hypothesis.

Plan section 2 is explicit: *"Do not add a fifth post-hoc hypothesis if the first
four fail."* Restructuring the paper around BIM would be that move wearing a
result's clothing. The prohibition is the reason the pre-registration is worth
anything, and it binds hardest exactly here - when the unexpected thing is more
interesting than the expected one.

So: reported as a specificity failure of the mouse model, in the informative
direction. Not promoted, not given a panel, not written into the abstract.

**Corrected reading, after section 0(a).** The prohibition on a fifth post-hoc
hypothesis still binds, and BIM still gets no panel until it replicates. But the
accurate description is no longer "a negative control that fired". It is:

> The mechanism H1 encodes - `MYC x OXPHOS` raises BH3-only-to-BCL-XL priming -
> **is supported in human tumours, with BIM rather than PUMA as the effector.**
> The specific endpoint pre-registered in this arm, PUMA/BCL-XL, is not
> supported. BIM was a co-candidate named in advance by the mouse arm on
> independent protein data; this arm's plan mis-classified it.

Three caveats travel with that sentence wherever it goes: the BIM fits are S1
only and untested across the D7 proliferation specifications; the mouse's own
*RNA* evidence for BIM was weak (p = 0.46), so the prior is protein-level; and
the independent-cohort replication in section 9 has not run.

## 9. PRE-REGISTERED REPLICATION, declared 2026-08-29

Declared **before METABRIC, SCAN-B or any neoadjuvant cohort has been opened in
this repo**, so it is a prediction rather than a description. Script 17 will test
it; nothing else may.

**Prediction.** In each independent cohort, scored in its own cohort-relative run
and meta-analysed as effect estimates (never pooled scores):

```
log2(BCL2L11) ~ M_a * OXPHOS_subunits + <the Block C covariate set>
    PREDICTED: interaction POSITIVE

log2(BCL2L1)  ~ M_a * OXPHOS_subunits + <the Block C covariate set>
    PREDICTED: interaction NEGATIVE

log2(BBC3)    ~ M_a * OXPHOS_subunits + <the Block C covariate set>
    PREDICTED: NULL   (this is the control that makes the other two mean something)
```

**Fixed now:**

- Both co-primary instruments; the claim requires both, as in Block C.
- All three D7 specifications (S1, S2, S3) reported side by side. The TCGA result
  has not been shown robust to the proliferation covariate and the replication
  must settle that.
- Direction is pre-stated above, so the tests are one-sided in interpretation
  even though two-sided in computation.
- The expression-matched null applies as in Block C: `OXPHOS subunits` must beat
  its matched null for the `BCL2L11` endpoint, or the result is a scale artefact.
- **Failure condition:** if `BCL2L11` does not replicate in the independent
  cohorts with the direction above, it is reported as a TCGA-specific
  observation and dropped. No further variants are tried.

This declaration is what separates a pre-registered replication from a story
built backwards. If it is not honoured, delete it rather than amend it.

## 10. Where the falsification criteria actually stand

Plan section 2 requires **all four** to hold before the human arm is called a
failure:

| Criterion | Status |
|---|---|
| `MYC:OXPHOS` on PRIME null in the full cohort **and** every stratum | **MET** |
| No `BUFFER` enrichment at CNV or expression level | **NOT MET** - G2 passed (`MCL1` OR 1.69, `BCL2L1` OR 1.77). Expression level is Block B, not yet run |
| FOXO3 regulon shows no relationship to OXPHOS in any stratum | **MET** in these fits |
| H4 fails in the informative direction | **UNTESTED** - Block F |

So the arm is **not** fully falsified. H1's coupling clause fails; H1's buffering
clause is supported at the DNA level and untested at the expression level; and
the conditional-chemosensitivity prediction has not been looked at.

## 11. What happens next

1. **Block B (script 10).** Pre-specified, cheap, and it settles the one
   remaining TCGA clause above. Run before anything else.
2. Then the Block F decision. Plan section 12 puts the F3 incremental-value check
   against MB1 forkscale ahead of any outcome modelling, and it is a stop gate.
3. Plan section 15's survivor-bias framing is now live: established tumours are
   survivors of the bottleneck, so the absence of the coupling is consistent with
   the window having closed, and the Hannon pre-malignant arm is the proper test.
   **With one complication that must be written into the text rather than
   smoothed over:** something *is* coupled to `MYC x OXPHOS` in established human
   tumours - BCL-XL down, BIM up - just not the effector the mouse identified.
   "The window has closed" and "the wiring is different" are different claims and
   this result does not distinguish them.
