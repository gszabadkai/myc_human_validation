---
date: 2026-08-31
status: DRAFT for the author. Manuscript text, not an analysis note.
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 4, 15)
  - every dated result note in docs/
purpose: the human arm as a well-formed negative, ready to cut into the Letter
---

# The human arm, written up

Plan section 15's framing, which the evidence now supports completely. Nothing
below claims anything the record does not carry.

---

## 1. The recommendation first

**The human arm should be one Extended Data figure and a single main-text
paragraph. It should not take a main display item.**

Plan section 4 allocated it three panels and flagged that the paper was already
at or past four display items before the human arm existed. That tension now
resolves itself: an arm in which no hypothesis is supported and the one
surviving finding failed its own replication does not earn a main figure. It
earns a short, exact paragraph and the Extended Data to back it.

This is not modesty. A main panel showing a null interaction invites the reader
to ask what it is doing there, and the honest answer - "we predicted this and it
did not happen" - is a sentence, not a figure.

## 2. Main-text paragraph (draft, ~250 words)

> To ask whether the MYC-OXPHOS coupling operates in human disease, we tested
> it in 938 primary breast tumours (TCGA-BRCA) with a pre-registered analysis
> plan, fixing the hypotheses, the endpoint and the falsification criteria
> before any model was fitted. The endpoint was apoptotic priming,
> log2(BBC3) - log2(BCL2L1), and MYC x OXPHOS was measured on two independent
> instruments - a pathway-level expression score and a composition-based
> mitochondrial pathway prioritisation. The interaction was null on both
> (+0.020, p = 0.46; +0.010, p = 0.73), and OXPHOS subunits ranked at the 32nd
> percentile of 2,000 expression-matched random gene sets, so the null is not a
> power or scale artefact. It remained null in every pre-specified stratum,
> including TP53-mutant and PI3K-intact tumours, where the mouse mechanism
> predicts it should be strongest. We found no evidence that surviving tumours
> had bought that survival through BCL2-family buffering: MYC, MCL1 and BCL2L1
> amplifications co-occur (adjusted OR 1.69-1.77), but MCL1 expression showed no
> interaction with MYC x OXPHOS (p = 0.94, p = 0.33), and the MYC-high /
> OXPHOS-high quadrant was not enriched for buffering. Nor did the axis predict
> chemotherapy response: in 1,978 patients across three neoadjuvant cohorts, the
> pre-specified three-way interaction on pathological complete response was
> +0.195 (95% CI -0.069 to +0.459), and in the largest cohort it was zero.
> Three of four pre-specified falsification criteria are met.

**Notes on the paragraph.** It leads with the pre-registration because that is
what makes a negative publishable, and it gives the matched-null percentile
immediately because the first objection to any null is "you were underpowered."

## 3. The BIM sentence, and where it goes

One finding survived TCGA on both instruments and was the arm's only
pre-registered replication. It failed. That sequence is worth two sentences,
because a replication that was declared and then honoured is evidence about the
method as well as about the biology.

> A single association survived in TCGA on both instruments: MYC x OXPHOS was
> associated with lower BCL-XL and higher BIM, with PUMA unmoved (BIM +0.051,
> p = 0.005 and +0.068, p = 3 x 10^-4). BIM had been named in advance by the
> mouse arm on independent protein evidence, so we pre-registered a replication
> in an independent cohort before opening it. In 3,143 SCAN-B tumours the
> association was present and inverted (-0.036, p = 0.019 and -0.033, p = 0.037;
> I2 = 93% and 94% against the TCGA estimates), and we report it as a
> TCGA-specific observation.

**Do not soften "inverted" to "not replicated".** The two cohort estimates
exclude each other at 95% and the SCAN-B coefficient is more extreme than 99.5%
of expression-matched null sets. Saying only "did not replicate" would be
weaker than the evidence and would invite the reader to assume a power failure,
which it is not.

## 4. The interpretive paragraph (draft, ~120 words)

This is plan section 15's framing, and it must not be over-claimed. It is an
interpretation of a negative, offered as one.

> Established tumours are the survivors of the transformation bottleneck, and a
> coupling that operates during that bottleneck need not persist in the tissue
> that emerges from it. Human primary breast cancers are sampled long after the
> window in which the mouse model detects this axis, and the absence of the
> coupling in two large cohorts is consistent with that window having closed -
> though our data cannot distinguish a closed window from a mechanism that does
> not transfer to human tumours at all. Testing this requires pre-malignant
> human tissue with matched transcriptomes; the Hannon breast pre-malignancy
> cohort is the appropriate setting and was outside the scope of this study.

**The clause after the dash is not optional.** Without it the paragraph is a
rescue narrative. With it, it is a stated limitation and a named next
experiment.

## 5. Extended Data figure (one item, four panels)

| Panel | Content | Source |
|---|---|---|
| **a** | Coefficient forest: `MYC x OXPHOS` on PRIME, both instruments, beside the 17 pathway comparators and the four endpoint negatives, each with its expression-matched null percentile. Carries the result and its specificity in one panel. | Block C, script 09 |
| **b** | The strata: full cohort, TP53-mutant/wild-type, PIK3CA, PI3K-intact, PAM50, forkscale. All flat, with CIs. | Block C strata, Block D |
| **c** | H4: the three-way interaction on pCR per cohort and pooled, three neoadjuvant cohorts, with the random-effects summary and I2. | Block F1, script 13 |
| **d** | The BIM replication: TCGA and SCAN-B side by side, both instruments, with the matched-null distribution behind each. | Block C + script 17 |

Panel **d** is the one a referee will look at hardest, and it should be drawn so
the reversal is unmissable - two estimates with non-overlapping CIs on opposite
sides of zero, not a forest plot that averages them.

## 6. Methods statements that must appear

Short, and each exists because a referee will otherwise ask.

- **Pre-registration.** Hypotheses H1-H4, the endpoint, the specificity battery
  and the falsification criteria were fixed before any model was fitted; the
  dated analysis plan and decision notes are in the accompanying repository.
- **Two instruments.** OXPHOS was measured as a GSVA pathway level on
  variance-stabilised counts and as mitoPPS composition on DESeq2-normalised
  linear counts. Both are reported and **no claim rests on one alone.**
- **Cohort-relativity.** GSVA and mitoPPS are cohort-relative; each cohort was
  scored in a single run and effect estimates, never scores, were compared
  across cohorts.
- **The matched null.** Every OXPHOS result is calibrated against 2,000 random
  gene sets matched on set size and expression ventile, drawn and scored
  identically to the observed arm.
- **The BIM replication.** Declared 2026-08-29 before any replication cohort was
  opened; the cohort, covariate mapping, specifications and failure condition
  were fixed 2026-08-31 before any SCAN-B expression value was read.
- **The covariate calibration.** SCAN-B lacks purity, leukocyte fraction, TP53
  status and plate. Before scoring SCAN-B we refitted the TCGA model on the same
  938 patients under the reduced covariate set; the BIM estimate changed by
  about 10% and retained sign and significance, so the covariate reduction does
  not explain the reversal.
- **Symbol harmonisation.** SCAN-B is annotated against a 2014 gene build; gene
  sets were harmonised to it through MitoCarta's curated synonyms in the forward
  direction only, and 19 of 89 OXPHOS subunit genes required it.

## 7. Limitations, stated rather than buried

1. **The two cohorts genuinely disagree on BIM, and we cannot say why.** SCAN-B
   is population-based and screening-era, ~75% ER-positive and largely
   early-stage; TCGA-BRCA is a referral series with more advanced disease and
   nearly twice the basal fraction (17.0% vs 10.1%). Any of these could matter.
   We did not test them, because a subgroup selected after seeing a failed
   replication is not evidence.
2. **PRIME is close to a PUMA readout.** It correlates 0.838 with its numerator
   and -0.128 with its denominator, so "PRIME is flat" and "PUMA is flat" are
   nearly the same statement. The limb-wise fits are reported for this reason.
3. **Bulk transcriptomes cannot see apoptotic priming directly.** The endpoint
   is a transcript ratio standing in for a protein-level threshold, and the
   mouse evidence for BIM is protein-level.
4. **The functional test was null at adequate power, not underpowered**, but
   DepMap lines are not tumours.
5. **Purity and immune infiltrate are the dominant confounders in breast**, and
   are adjusted for in TCGA but unavailable in SCAN-B (see the calibration).

## 8. Sentences that must NOT be written

Recorded because each is available and each would be wrong.

- *"MYC x OXPHOS suppresses BIM in ER-positive disease."* A new hypothesis built
  backwards from a failed replication. Forbidden by the plan.
- *"The coupling is restored in [any subgroup]."* No subgroup was pre-specified
  for the reversal and none was tested.
- *"BIM is breast-specific."* The DepMap result was UNRESOLVED, and "not
  universal, therefore specific" is a post-hoc rescue.
- *"MCL1 predicts chemoresistance in TNBC."* One gene, one cohort, I2 = 82%,
  Q p = 0.0036 across three.
- *"The axis is confirmed in human tumours."* It is not, on any endpoint, on
  either instrument, in either cohort.
- Any framing in which the human arm *supports* the mouse. It does not. Its
  value is that it was capable of saying so.

## 9. What this arm contributes to the paper

Stated plainly, because it is easy to lose.

The mouse and cell work is the paper. The human arm's contribution is that a
pre-registered, adequately powered, two-instrument test in 4,081 primary
tumours and 1,978 neoadjuvant patients **did not find the axis in established
human breast cancer**, and said so. That bounds the claim the paper makes: the
mechanism is demonstrated where it was demonstrated, and is not asserted beyond
it.

A reviewer who asks "does this hold in humans?" now has an answer with a date
on it, rather than a selected panel.
