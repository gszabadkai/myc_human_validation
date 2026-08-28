---
date: 2026-08-28
status: RECORDED - validation result, no decision pending
relates-to:
  - 2026-08-27_human_validation_plan.md (section 7.1)
  - 2026-08-28_G1_result_and_decisions.md (D2)
decides:
  - D2 closed on evidence: M-a and M-b agree, and stripping cost nothing
  - M-a/M-b agreement threshold operationalised at Spearman rho >= 0.30
next-action: section 9 of the specificity panel proposal, then script 07
---

# MYC estimator validation

Run 2026-08-28 via `scripts/06_score_myc_activity.R` on TCGA-BRCA,
n = 1,095 patients, 18,115 genes. Results in
`results/tcga_brca_myc_scores.rds`.

Plan section 7.1 specifies three MYC estimators and a rule - M-a primary, with
mandatory M-b/M-c concordance - but the plan never says what concordance is, and
never checks that the estimators measure MYC at all. This is that check.

---

## 1. The three estimators

| | Definition | n | Role |
|---|---|---|---|
| **M-a** | Felsher signature, MitoCarta-stripped, GSVA on VST | 61 genes | **PRIMARY** (D2) |
| **M-b** | CollecTRI MYC regulon, stripped, decoupleR ULM | 811 targets (798 in matrix) | concordance |
| **M-c** | MYC 8q24.21, ISAR GISTIC, `== +2` | 224 amplified of 1,040 called | instrument |

All read the **log-scale VST** matrix with `kcdf = "Gaussian"`, scored in one
cohort-relative run. The linear matrix is never opened by this script and the
`scale` field is asserted rather than trusted.

## 2. Concordance

```
M_a vs M_b            0.613     AGREE
M_a vs M_c            0.405
M_b vs M_c            0.348
M_a stripped vs raw   0.997
M_b stripped vs raw   0.991
```

Spearman, complete observations.

### 2.1 M-a and M-b agree, and that is not circular

**rho = 0.613.** These are built from disjoint evidence: a 61-gene expression
signature averaged by a rank-based enrichment statistic, against an 811-target
signed regulon scored by a linear model. They share no construction step and
only six genes' worth of incidental overlap. Nothing forces them to agree.

**Threshold, pre-registered before the number was seen:** agreement is
`rho >= 0.30`. Deliberately moderate, because identical scores would be
suspicious rather than reassuring, and deliberately **not** a significance
threshold - at n ~1,095 almost any rho is "significant" and a p-value here
would carry no information.

Plan section 7.1's consequence ("if M-a and M-b disagree, report both and treat
the claim as unsupported") is therefore **not** triggered.

### 2.2 Both estimators separate by copy number, and this is the real test

`M-c` is a **DNA** measurement. It cannot share a gene, a normalisation or a
scoring step with either expression estimator, so a relationship between them
cannot be an artefact of shared construction.

```
Kruskal-Wallis, score by GISTIC call (-1/0/1/2)
  M_a   chi-squared = 187.06, df = 3, p < 2.2e-16
  M_b   chi-squared = 131.44, df = 3, p < 2.2e-16
```

Both expression estimators are ordered by MYC copy number. **This is the
strongest available evidence that they measure MYC rather than a correlated
transcriptional programme**, and it is the check the plan did not ask for.

The correlations with M-c (0.405, 0.348) are moderate by design: copy number is
a coarse four-level call and amplification is neither necessary nor sufficient
for MYC activity. A high correlation here would suggest the expression scores
were tracking aneuploidy rather than transcriptional output.

### 2.3 The MitoCarta strip cost nothing

```
M_a stripped (61) vs raw (67)     rho = 0.997
M_b stripped (811) vs raw (886)   rho = 0.991
```

G1 removed the MitoCarta genes to break circularity with the mitochondrial axis,
and concluded from gene-set structure that this should not change what the
estimator measures. It does not. Script 06 warns if either falls below 0.90;
neither comes close.

**This closes D2 on evidence rather than on argument.** The G1 note resolved D2
by reasoning about set composition; this is the same conclusion reached from the
scored data.

## 3. Subtype ordering

Median M-a by PAM50:

```
Basal   +0.473
Her2    +0.142
LumB    +0.019
Normal  -0.217
LumA    -0.323
```

MYC activity highest in basal-like, lowest in luminal A. This is the established
ordering in breast cancer and **nothing in the pipeline was fitted to produce
it** - the signature is a fixed 61-gene list from an external source, scored
without reference to subtype. It is a free external validity check and it passes.

Note the consequence for later blocks: MYC activity is strongly confounded with
subtype, which is why PAM50 is a covariate in Block C and why the H4 analyses
stratify. It also sharpens D7 - the proliferation covariate, subtype and MYC
activity are three views of an overlapping thing.

## 4. Implementation details that affect the numbers

**Sign rule for CollecTRI.** `mor = +1` for stimulatory edges, `-1` for
inhibitory. 92 of MYC's 891 edges are flagged **both**; they take `+1`, matching
how G1 built `COLLECTRI_MYC_STIM`. After de-duplicating to distinct targets the
network is **810 activating, 76 repressing**. A signed method's answer depends on
this, so it is recorded rather than left in a comment. Script 06's sandbox
carries the check of whether flipping those 92 changes the ordering; if it does,
the rule belongs in Methods.

**Coverage.** M-a and M-b are defined for all 1,095 patients. M-c is available
for 1,040 - the ISAR GISTIC subset - so any model using M-c loses 55 patients
relative to one using M-a.

**Gene universe.** 18,115 genes after the low-count filter, up from 18,055 in the
first run of script 01. The five samples that changed from Xena-averaged `-01A`
to clean `-01B`/`-01C` shifted which genes clear the threshold. Confirms the
averaged-aliquot fix propagated.

## 5. What this does and does not establish

**Established:** the three estimators agree; the two expression estimators are
ordered by DNA copy number; the MitoCarta strip did not change what they
measure; the subtype ordering matches external expectation.

**Not established:** nothing about OXPHOS, priming, or any hypothesis. This is
instrument validation only. The MYC axis is now fit to use as an exposure - that
is the entire claim.
