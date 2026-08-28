---
date: 2026-08-28
status: D5 RESOLVED
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 2, 3, 7.5, 10, 11, 14)
  - 2026-08-28_D4_scoping_and_G2_design.md
  - 2026-08-28_G2_result.md
decides:
  - D5 - RESOLVED. I-SPY2 (GSE194040) primary, BrighTNess replication,
    GSE25066 third, I-SPY1 dropped
  - H4 primary test is CONTINUOUS; the median-split STATE version becomes the
    pre-specified secondary. Amends plan section 7.5.
  - H4 uses all three cohorts, META-ANALYSED. Pooling scores stays forbidden.
next-action: D7, then script 01
---

# D5 - which neoadjuvant cohort is primary for H4

---

## 1. Why this was reopened rather than confirmed

Plan section 14 asked me to confirm GSE25066 as primary with BrighTNess as
subtype-matched replication. **The candidate set turned out to be undocumented.**
`GSE25066, GSE22226, GSE164458` appears exactly once in the repository - one table
cell in plan section 3, present in the initial commit - with no search record, no
inclusion criteria and no selection log. Confirming a choice inside an unexamined
set is not a decision, and H4 is the arm's main-panel counter-intuitive claim.

Two further facts made reopening necessary:

- **Lee et al. 2017 used none of these cohorts.** Their evidence is TCGA (n=816),
  METABRIC TNBC (n=320), CCLE (n=20), a NanoString cohort (n=34, 17 matched
  pre/post pairs) and an 18-patient targeted DNA panel. They tested what survives
  chemotherapy, never response prediction. So there is no precedent cohort to
  match, and no risk of re-running them.
- **G2 changed the subtype argument.** The plan's stated reason for a TNBC-matched
  replication was matching MB2_UF. But G2 showed `MCL1` is basal-specific while
  `BCL2L1` - the term carrying Panel b - is subtype-BROAD (Breslow-Day by TNBC
  p=0.85). A TNBC-only primary would test the novel axis in the one stratum where
  its breadth is invisible.

## 2. Search method, fixed before looking

Inclusion criteria, derived from what H4 needs and written down before screening:

| | Criterion |
|---|---|
| **In** | human breast cancer; **pre-treatment** tumour material; genome-wide expression (array or RNA-seq); systemic **neoadjuvant chemotherapy**; response endpoint **pCR and/or RCB**; n >= 100; publicly accessible |
| **Out** | post-treatment/residual only; metastatic; endocrine-only neoadjuvant; targeted panels (cannot support OXPHOS GSVA) |

Three channels:

1. **GEO E-utilities** - the dataset universe. 447 human breast neoadjuvant
   expression series retrieved; **34 passed the screen**.
2. **PubMed** - 122 hits on a multi-omic/neoadjuvant/pCR query, top 25 examined,
   to catch cohorts hosted outside GEO.
3. **Consensus** - 20 papers, review-level coverage.

**Consensus returned no cohort GEO had not already shown, and did not surface
I-SPY2, BrighTNess or NeoTRIP at all.** It indexes papers, not data resources.
Recorded so nobody repeats it expecting dataset discovery.

PubMed contributed one genuine addition GEO could not: **TransNEO**
(Sammut et al. 2022, Nature 601:623, doi:10.1038/s41586-021-04278-5), n=168 with
clinical, digital pathology, genomic and transcriptomic profiling plus a 75-patient
validation set. The deepest neoadjuvant cohort in existence, but **EGA controlled
access** and too small to be primary.

**Constraint held throughout: no score-versus-pCR association was computed in any
candidate cohort during selection.** Marginals only - n, platform, endpoint coding,
subtype composition, regimen, gene coverage. Selecting a cohort on an outcome
association would invalidate the pre-registration.

## 3. The plan's own facts were partly wrong

| | Plan says | Verified |
|---|---|---|
| BrighTNess | "smaller and more recent" | **n=482** vs GSE25066's 508 - the same size, and RNA-seq rather than a 2011 array |
| I-SPY1 GSE22226 | a candidate | **n=150**, Agilent two-colour, superseded by I-SPY2 |

## 4. Decision

| Rank | Cohort | n | Design | Rationale |
|---|---|---|---|---|
| **PRIMARY** | **GSE194040** I-SPY2-990 | **988** | randomised, 14 arms, 179-patient paclitaxel control | largest; only candidate able to test `score x treatment`; RPPA companion in the same patients; 4 `BBC3` probes |
| **Replication** | GSE164458 BrighTNess | 482 | randomised, TNBC | RNA-seq; tests whether BCL2L1's subtype-breadth holds in the matched stratum |
| **Third** | GSE25066 | 508 | single-arm taxane-anthracycline | genuinely independent - different platform generation and era |
| Deferred | TransNEO | 168 | observational | deepest multi-omics; revisit only if EGA access is sought |
| **Dropped** | GSE22226 I-SPY1 | 150 | - | superseded, and much smaller than the plan implied |

### Why I-SPY2 displaces GSE25066

1. **GSE25066 is single-arm.** Everyone receives taxane-anthracycline, so it can
   never test the `score x treatment` interaction that is half of what H4 asks for.
   I-SPY2 has 14 arms and a 179-patient paclitaxel-only control.
2. **n=988 with 319 pCR events** against 508.
3. **RPPA companion GSE196093** - 139 signalling proteins/phosphoproteins in 736 of
   the same patients. Plan Block E / script 16 currently routes RPPA through TCGA,
   i.e. different tumours. This puts it in the same ones.
4. **`BBC3` coverage.** 4 probes on GPL30493 against GSE25066's single probe
   `211692_s_at`, which multi-maps to `MIR3190`/`MIR3191`. `PRIME`'s numerator is
   better measured.

### Verified I-SPY2 metadata

All 988 samples are `Breast cancer biopsy (pre-treatment)`. pCR is complete for
every sample: **319 pCR (32.3%)**, 669 no-pCR. Fields: `patient id`, `tissue`,
`hr`, `her2`, `mp`, `pcr`, `arm`.

| Subtype | n | pCR |
|---|---|---|
| HR+ HER2- | 379 | 16.9% |
| TNBC | 364 | 39.0% |
| HR+ HER2+ | 156 | 36.5% |
| HR- HER2+ | 89 | 62.9% |

Arms range from paclitaxel alone (n=179, 17.3% pCR) to paclitaxel + pertuzumab +
trastuzumab (n=44, 59.1%). **Arm must be adjusted for or stratified**; several arms
are HER2-targeted and therefore confounded with subtype.

Two cautions carried forward:

- The deposited matrix is labelled `adjustment method: ComBat Adjust`, i.e. it is
  **already batch-corrected**. That constrains upstream processing and must be
  stated in Methods.
- The `mp` field is roughly balanced (484 / 504), so it is **not** a high-risk
  enrolment flag as first assumed; it most likely encodes MP1 vs MP2 within
  I-SPY2's high-risk criterion. **Confirm before use.** I-SPY2 still enrols
  high-risk disease and is not a representative breast cancer population.

## 5. Power - the part that changed the design

Simulation, three-way `MYC x OXPHOS x BUFFER` on pCR, baseline 32.3% (I-SPY2
observed), median splits, alpha 0.05. `rr` = pCR risk ratio in H4's single predicted
cell versus the other seven.

| Cohort | n | rr=1.5 | rr=2.0 | rr=2.5 | rr=3.0 |
|---|---|---|---|---|---|
| GSE25066 | 508 | 15% | 36% | 59% | 78% |
| BrighTNess | 482 | 14% | 32% | 58% | 77% |
| **I-SPY2** | **988** | 22% | 58% | **87%** | 98% |
| all three | 1978 | 38% | **86%** | 99% | 100% |

**I-SPY2 alone reaches 80% only for a large effect** (rr >= 2.5: ~62% pCR in the
target cell against ~25% elsewhere). Large, but subtype alone already spans
16.9% to 62.9% in this cohort.

### 5.1 The two-way is NOT a fallback

The plan implies that an underpowered three-way can retreat to `MYC x OXPHOS`.
Tested:

```
I-SPY2, rr = 2.0:    three-way 57.7%    two-way 62.5%
```

**Five percentage points.** The binding constraint is the single-cell contrast, not
the number of terms. Dropping BUFFER sacrifices the entire distinction from Lee et
al. and buys almost nothing. Do not write it into the plan as a fallback.

### 5.2 Median splits are the largest avoidable loss

Same truth (continuous `M`, `O`, `B` with a genuine three-way term), two analysis
strategies, n=988:

| true 3-way coefficient | continuous | median-split | loss |
|---|---|---|---|
| 0.15 | 53.8% | 19.4% | **34 pp** |
| 0.25 | 91.9% | 39.3% | **53 pp** |
| 0.35 | 99.7% | 68.7% | 31 pp |
| 0.50 | 100.0% | 92.9% | 7 pp |

Up to **53 points of power**, recoverable for nothing but a specification choice.

## 6. Consequent decisions

### 6.1 AMENDMENT to plan section 7.5 - continuous becomes primary

Section 7.5 currently reads: *"Continuous version as sensitivity, not as an
alternative to be swapped in if the categorical fails."* **That is inverted, as
follows:**

> **Primary H4 test is the continuous three-way interaction.** The four-level
> `STATE` factor and its level-3-vs-level-4 contrast become the **pre-specified
> secondary**, reported alongside.

**Why this is not the drift section 7.5 prohibits.** The prohibition is against
swapping to a continuous analysis *after the categorical one fails* - a post-hoc
rescue. This amendment is made **before script 11 exists, before any H4 outcome
data has been touched, and on the basis of a simulation containing no project
data.** The reason is power, not result. Both specifications are fixed now and both
will be reported. Recorded here with a date so the ordering is auditable.

`STATE` is still frozen in script 11 and still must not be revised after outcome
data is seen. Its role changes from primary to secondary; its definition does not.

### 6.2 H4 uses all three cohorts, meta-analysed

Pooling reaches 86% power at rr=2.0. But CLAUDE.md forbids pooling GSVA scores
across separately-scored cohorts - they are cohort-relative and not comparable.
**Meta-analysing the three interaction estimates is both the powered route and the
only permitted one**, and Block F already specifies that machinery.

So H4 is not "one cohort with two afterthoughts". It is three cohorts scored
independently, three interaction estimates, meta-analysed. Report the per-cohort
estimates and the pooled estimate with heterogeneity.

### 6.3 Honest statement of what is detectable

Pre-registered now, before data: **H4 is powered to detect a large conditional
effect and is not powered to detect a modest one.** At rr=2.0 even the
meta-analysis sits at 86% while any single cohort is below 60%. Report confidence
intervals rather than leaning on p-values, per plan section 3's power caveat. A
null H4 will be reported as "not powered to exclude a modest effect", not as
falsification of the model.

## 7. Still open

- **D7** - proliferation covariate overlaps the MYC exposure. PROPOSED. Confirm
  before script 09.
- **G1 correlation criterion** - needs script 01. G1 remains not fully discharged.
- **Script 07 specificity panel source** - undecided, plus the FAO interpretive note
  from D4.
- **`TCGAbiolinks`** - not installed. Needed for scripts 01/02/03.
- **`mp` field semantics** in I-SPY2 - confirm before use (section 4).
- **Conditional enrichment** for G2 - proposed in the G2 note, not run.
