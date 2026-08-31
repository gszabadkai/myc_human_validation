---
date: 2026-08-31
status: DECLARATION - written before any SCAN-B expression value has been read
relates-to:
  - 2026-08-29_block_c_result_H1_not_supported.md (section 9, the declaration being honoured)
  - 2026-08-31_block_F1_result_H4_not_supported.md (section 9, the three options)
  - 2026-08-31_handoff.md (section 3)
  - 2026-08-28_D7_proliferation_covariate.md (S1/S2/S3)
  - 2026-08-29_script09_build_spec.md (the spine, the null, the percentile)
decides:
  - Option 1 taken. The BIM replication runs in SCAN-B, cohort GSE202203.
  - The Block C covariate set is mapped, not reproduced. The mapping is fixed here.
  - A mandatory TCGA calibration gates the SCAN-B fits.
next-action: script 16 (fetch SCAN-B), then script 17 (the replication)
---

# The BIM replication, declared for SCAN-B

The Block C note section 9 registered a replication and then bound it:

> *"If it is not honoured, delete it rather than amend it."*

Amendment A1 established that mitoPPS does not exist in the three neoadjuvant
cohorts, which blocked the declaration's co-primary clause there. This note takes
option 1 from the F1 note section 9 - **run it where both instruments exist** -
and fixes every remaining free parameter before any expression value is read.

Everything below is settled in advance. Nothing here has been fitted.

---

## 1. The cohort, and a correction to the plan's assumption

The plan (section 3) names SCAN-B as "RNA-seq, large n" and the handoff recorded
it as "RNA-seq with counts, so both instruments are possible". **That was an
assumption and for the deposit the plan implies it is false.**

Checked 2026-08-31 against the GEO FTP listings and the series-matrix processing
fields:

| Accession | n | Deposited | Counts? |
|---|---|---|---|
| GSE96058 (Brueffer 2018) | 3,273 | `gene_expression_..._transformed.csv.gz`, cufflinks FPKM, log2(FPKM+0.1) | **no** |
| GSE81538 (Brueffer 2016) | 405 | `gene_expression_405_transformed.csv.gz`, FPKM | **no** |
| **GSE202203** (Staaf group, ESR2) | **3,207** | `RawCounts_gene_3207.tsv.gz` **and** `TPM_Raw_gene_3207.tsv.gz` | **yes** |

Under A1's own logic, GSE96058 would have been blocked for exactly the reason the
neoadjuvant cohorts are: a deposit that arrives already logged cannot carry
mitoPPS. **GSE202203 is the cohort, and it is the only SCAN-B deposit on which
this declaration can be honoured.**

GSE202203 and GSE96058 are overlapping draws from the same SCAN-B study. Only
GSE202203 is used here; nothing in this arm reads GSE96058, so there is no
double-counting to guard against.

### 1.1 Independence, stated honestly

SCAN-B is independent of TCGA - different country, different accrual, different
platform generation, no patient overlap. It is **not naive to the companion
paper**: Menegollo, Bentham et al. used SCAN-B in Fig 7H-L. That was the
forkscale/bicluster analysis. The `MYC x OXPHOS` interaction and the `BCL2L11`
endpoint are new to this cohort. This is a note for Methods, not a circularity.

## 2. What has been looked at, and what has not

For the audit trail, because this note claims to precede the data.

**Read:** the GEO FTP supplementary listings; the series-matrix headers
(`!Series_*`, `!Sample_data_processing`, `!Sample_characteristics_ch1` field
*names* and their first few values); and the first ~3 MB of
`GSE202203_RawCounts_gene_3207.tsv.gz`, to establish the file layout, the column
count (3,208 = 3,207 samples + a symbol column) and that the values are
non-integer.

**Not read:** any expression value for any gene named in this declaration, any
outcome, any model. No file has been downloaded into `data/`.

## 3. What is being replicated

From Block C, TCGA, n = 938, spine specification, S1:

```
                 GSVA               mitoPPS
log2(BBC3)      -0.014 (p 0.60)    -0.034 (p 0.21)    PUMA:   unmoved
log2(BCL2L1)    -0.034 (p 0.030)   -0.043 (p 0.009)   BCL-XL: DOWN, both
log2(BCL2L11)   +0.051 (p 0.005)   +0.068 (p 3e-4)    BIM:    UP, both
```

The registered prediction, restated verbatim from the Block C note section 9:

```
log2(BCL2L11) ~ M_a * OXPHOS_subunits + <the Block C covariate set>
    PREDICTED: interaction POSITIVE

log2(BCL2L1)  ~ M_a * OXPHOS_subunits + <the Block C covariate set>
    PREDICTED: interaction NEGATIVE

log2(BBC3)    ~ M_a * OXPHOS_subunits + <the Block C covariate set>
    PREDICTED: NULL   (this is the control that makes the other two mean something)
```

Direction is pre-stated, so the tests are **one-sided in interpretation and
two-sided in computation**, as declared.

## 4. Only M-a. M-b is forbidden for this endpoint

Block C section 7 records that `BCL2L11` **is** a CollecTRI MYC target and would
therefore be circular with M-b. Every Block C fit behind the numbers above used
M-a, and every fit here does too.

`M_a` = GSVA over the 61-gene MitoCarta-stripped Felsher set (`g1$estimators_stripped$FELSHER`),
scored in SCAN-B's own run. M-c (`MYC_amp`, GISTIC) does not exist in SCAN-B - no
CNV is deposited - and is simply unavailable. **The MYC estimator is not a free
parameter in this replication.**

## 5. The covariate mapping - the one substantive decision here

The TCGA spine is:

```
purity + leukocyte_fraction + PROLIF_DISJOINT + PAM50 + TP53_status + plate_pooled
```

Four of six do not exist in SCAN-B. The mapping is fixed **now**, before fitting:

| Spine term | SCAN-B | Decision |
|---|---|---|
| `purity` | no ABSOLUTE call | **drop** |
| `leukocyte_fraction` | TCGA methylation-derived | **drop** |
| `PROLIF_DISJOINT` | computable, same 318 genes | **keep** |
| `PAM50` | `pam50 subtype` (Parker method, fixed reference) | **keep** |
| `TP53_status` | no sequencing in this deposit | **drop** |
| `plate_pooled` | no plate or batch deposited | **drop** |

```
SCAN-B PRIMARY:  log2(GENE) ~ M_a * OXPHOS_subunits + PROLIF_DISJOINT + PAM50
```

**Purity is dropped, not estimated.** ESTIMATE or any deconvolution would put a
different instrument into the model than the one Block C used, and a replication
whose covariate is measured differently is not a replication of the same
specification. The arm already has a reduced-covariate precedent: script 09's
D8 tiers `M2` and `M3` drop purity and PAM50 respectively.

**SCAN-B's PAM50 is its own call**, not the PanCanAtlas call TCGA used. Same
construct, different implementation. Adjusting for subtype is the point; which
centroid assigned it is not.

### 5.1 Specifications reported, never selected between

The declaration requires all three D7 specifications side by side. A fourth is
declared here and fenced.

| Label | Covariates | Status |
|---|---|---|
| **S1** | `PROLIF_DISJOINT + PAM50` | **primary** |
| **S2** | `PAM50` only (proliferation dropped) | declared sensitivity |
| **S3** | `PROLIF_STD + PAM50` | declared sensitivity |
| **C-alt** | S1 `+ age + NHG` | **beyond the declared set** |

`C-alt` exists because SCAN-B carries age and grade and a reviewer will ask why
they were unused. It adds covariates the Block C spine never had, so it **cannot
be the primary and cannot be selected into one**. It is reported alongside, and
if it disagrees with S1 that disagreement is reported as such.

`4 specifications x 2 instruments x 3 endpoints = 24 fits.` That is the whole
model space. There are no others.

## 6. Step 0 - the TCGA calibration, and it is a hard gate

Dropping four covariates is not free, and without this step a SCAN-B null would
be uninterpretable: it could mean BIM does not replicate, or it could mean the
covariate reduction did it. **Those must be distinguishable before SCAN-B is
fitted, and they can be, in TCGA, at no cost.**

Refit the three Block C limbs on the **same 938 TCGA patients, same frozen
z-scaling constants, both instruments**, with the covariate set reduced to
`PROLIF_DISJOINT + PAM50` - the SCAN-B mapping exactly.

**Pass condition, fixed here:** for `BCL2L11`, on **both** instruments, the
reduced-set estimate must (a) keep the positive sign, and (b) have the
full-spine point estimate inside its 95% CI.

If it fails, the reduced covariate set is a confounder of the replication, the
replication is reported as **not runnable as declared**, and the arm falls back
to the F1 note's option 2 (delete the commitment) with this note as the reason.
Script 17 stops at that point and fits nothing in SCAN-B.

The pass condition is a statement about the two estimates being
indistinguishable, not an arbitrary tolerance, and it is stated before either is
computed.

## 7. Scale discipline in this cohort

The two instruments have opposite requirements and must not share an object.
The pattern is script 01's and script 07's, unchanged:

```
RawCounts (non-integer)
  -> round()                      -> integer counts
  -> DESeq2 size factors
       |-- counts(dds, normalized = TRUE)   LINEAR   -> mitoPPS
       |-- vst(dds)                          LOG      -> GSVA (kcdf = "Gaussian")
```

**The counts are not integers.** They are StringTie `-e` estimated counts
(`A1BG 412.35820192`), the prepDE-style output. DESeq2 requires integers, so they
are rounded. Gene-level counts are all that is deposited, so tximport with
transcript lengths is not an option. This is a documented deviation from TCGA,
where GDC counts are already integers; it is recorded here rather than discovered
in a script.

**Cohort-relativity.** Every score - `M_a`, `OXPHOS subunits`, `PROLIF_DISJOINT`,
`PROLIF_STD` - is computed in one SCAN-B run over the SCAN-B matrix. No score
crosses cohorts. Exposures are z-scored **within SCAN-B**, so coefficients are
per SD of the SCAN-B distribution: comparable to TCGA in sign and in rough
magnitude, not identical in units. This is the correct handling under the GSVA
cohort-relativity rule and it is why the comparison is of directions and CIs, not
of decimal places.

**GSVA pin sets.** `.PIN_A` / `.PIN_B` (the two complementary half-matrix sets)
are carried over unchanged. Without them a batch of 20 sets and a batch of 2,000
walk different gene universes and the null is silently invalid.

**mitoPPS universe.** All MitoCarta pathways with >= 3 genes present, mtDNA split
out, the two synthetic mouse Apoptosis pathways still excluded for the reasons in
script 07 section 5.1. mitoPPS values are **never** compared numerically between
SCAN-B and TCGA - only the sign and significance of the interaction transfer.

## 8. The expression-matched null

The declaration requires it: *"`OXPHOS subunits` must beat its matched null for
the `BCL2L11` endpoint, or the result is a scale artefact."*

Scope is exactly that, and no wider: **one arm** (`OXPHOS subunits`), **one
endpoint** (`BCL2L11`), **both instruments**, **S1 only**. 2,000 draws, 20
ventiles of mean linear expression, drawn without replacement within a ventile -
`G7_NULL_NSETS`, `G7_NULL_NBIN` and the draw function from script 07 unchanged.
Ventiles are computed on SCAN-B's own expression, because a null matched on
TCGA's abundance distribution is not matched here.

Percentile and empirical p as in build-spec section 3:

```
percentile(A) = 100 * mean(b_null < b_obs)
p_emp(A)      = 2 * min( mean(b_null <= b_obs), mean(b_null >= b_obs) )
```

The `.OBS_CHECK` assertion - re-score the observed arm inside its own null batch
and stop if the value moves - is carried over. `mtDNA-encoded OXPHOS` remains
outside the null (`G7_NULL_SKIP`) for the abundance reason already recorded.

## 9. What counts as replication, fixed before the fit

**`BCL2L11` replicates if and only if:**

1. the `M_a:OXPHOS_subunits` interaction is **positive**, at p < 0.05,
2. on **both** co-primary instruments,
3. at **S1**,
4. and the arm beats its expression-matched null on the `BCL2L11` endpoint.

**One instrument alone is not a replication.** That is the rule the whole arm has
been run under and it does not soften here.

`BCL2L1` (predicted negative) and `BBC3` (predicted null) are reported with it
and interpreted as they were in Block C: `BCL2L1` corroborates direction,
`BBC3` moving would undercut the specificity of the other two.

### 9.1 There is no meta-analysis

The declaration says "meta-analysed as effect estimates". With SCAN-B as the only
cohort, **k = 1 and there is nothing to pool.** Script 13's DerSimonian-Laird
machinery is not invoked. This is a single-cohort replication and is described as
one.

### 9.2 The failure condition, restated and binding

> *"if `BCL2L11` does not replicate in the independent cohorts with the direction
> above, it is reported as a TCGA-specific observation and dropped. No further
> variants are tried."*

**No further variants are tried.** Not another cohort, not another covariate set,
not another estimator, not a subtype stratum. If it fails, it is written up as
having failed, and that is a better outcome for the manuscript than deleting the
commitment would have been.

## 10. What this replication does NOT do

The fences matter more here than usual, because GSE202203 is a rich deposit and
the arm is otherwise finished.

`GSE202203` carries `overall survival days/event`, `relapse free interval
days/event`, `endocrine treated` and `chemo treated`. **None of them is read by
scripts 16 or 17.** They are not parsed into the analysis frame.

Specifically out of scope, and each would be a separate decision:

- **F2** (treatment-stratified survival in SCAN-B). Not opened by this.
- **F3 / F3-pre** (forkscale). Not computed in SCAN-B.
- **F4**, subtype-specific claims, `STATE`, `PRIME`, `BUFFER_c`, the specificity
  battery beyond the single null arm, and the Block C excursion ladder.
- Any fifth hypothesis. Plan section 2 still binds and this is not one - it is
  the replication of a result registered on 2026-08-29.

SCAN-B is fetched for the BIM replication and for nothing else. If any of the
above is wanted later, it is a new decision with a new note, and the survival
columns will still be in the file.

## 11. Sample handling and the guards that must stop the run

Complete cases on what the model needs: `M_a`, `OXPHOS subunits`,
`PROLIF_DISJOINT`, `PAM50`, and the endpoint gene. `n` is reported, not
pre-set. No clinical or subtype restriction is applied.

Guards, each of which stops rather than warns:

- **Endpoint genes present.** `BCL2L11`, `BCL2L1`, `BBC3` must all be in the
  matrix. If `BBC3` is absent the control is gone and the other two cannot be
  interpreted - stop.
- **Identifier join.** Columns are `S000001`; the phenotype key is
  `scanb external id: Q009012.C009079.S000008.l.r.m.c.lib.g.k2.a`. Script 12's
  `.match_ids` is reused: score every declared transformation, stop below 90%.
  The transformation set is declared in script 16 and not tuned after seeing the
  match rate.
- **Duplicate samples.** GSE96058 shipped 136 replicates; GSE202203's file name
  claims none. If duplicate subjects appear, keep the first by external-ID order
  and **report how many were dropped**.
- **Non-numeric columns**, ragged headers, and the vroom connection-buffer limit:
  script 12's `.read_matrix_tsv`, `.numeric_matrix` and `.vroom_big` guards are
  reused wholesale. All three caught real failures in the neoadjuvant cohorts and
  none was visible in the output.
- **mitoPPS positivity.** Every MitoPathway score must be strictly positive in
  every sample or the pairwise ratio is undefined - script 07's stop.
- **mtDNA genes present.** GENCODE 27 protein-coding should carry the 13 `MT-`
  genes; confirm rather than assume, since the deposit is symbol-collapsed.

## 12. Build consequences

- **Script 16** - fetch GSE202203 (counts + series matrix), parse, join, round,
  DESeq2, VST, linear. Provenance README with URLs and SHA-256 under
  `data/scanb/`. No models.
- **Script 17** - section 2 is the step-0 TCGA calibration and a hard stop;
  sections after it score SCAN-B, fit the 24 models, run the null, and evaluate
  section 9's four conditions.
- `data/raw/` grows by ~250 MB and remains gitignored and unbacked-up. The README
  carries the re-download commands.

Both are written by Claude Code and run by the author in Positron, per Option A.

---

**Written 2026-08-31, before any SCAN-B expression value was read.** The ordering
is what makes this a replication rather than a description, and it is checkable
in the git log.

---

# Addendum A - symbol vocabulary, added 2026-08-31, still pre-fit

Added after the three files were downloaded and their symbol column read, and
**before any expression value entered any model**. It is recorded as an
addendum rather than folded into section 7 so the ordering stays checkable.

## A.1 What was found

SCAN-B is annotated against **UCSC knownGenes downloaded 22 September 2014**
(series `!Sample_data_processing`). Human MitoCarta 3.0 (2020) carries current
HGNC symbols. The ATP synthase subunits were renamed wholesale in 2018, so the
two vocabularies disagree on **19 of the 89 genes in `OXPHOS subunits`** - which
is the exposure of the declared model.

```
ATP5F1A <- ATP5A1    ATP5MC1 <- ATP5G1    ATP5PB  <- ATP5F1
ATP5F1B <- ATP5B     ATP5MC2 <- ATP5G2    ATP5PD  <- ATP5H
ATP5F1C <- ATP5C1    ATP5MC3 <- ATP5G3    ATP5PF  <- ATP5J
ATP5F1D <- ATP5D     ATP5MD  <- USMG5     ATP5PO  <- ATP5O
ATP5F1E <- ATP5E     ATP5ME  <- ATP5I     ATP5IF1 <- ATPIF1
ATP5MPL <- C14orf2   ATP5MF  <- ATP5J2    ATP5MG  <- ATP5L
DMAC2L  <- ATP5S
```

Unharmonised, `OXPHOS subunits` covers **0.787** here against 1.000 in TCGA, and
section 11's coverage floor would have stopped the replication - for the wrong
reason. The genes are all present under their 2014 names.

The failure mode of *not* noticing is worse than the stop: 70 of 89 genes,
still labelled `OXPHOS subunits`, is a Complex V with no F1 head and no c-ring.
Nothing downstream would have looked wrong.

## A.2 What is done about it, and why it is not an amendment

Script 16 section 7.3 harmonises the declared sets to the SCAN-B vocabulary
using **script 07 section 2's existing map**, applied to this matrix instead of
TCGA's. This is the arm's own symbol-harmonisation step running in a second
cohort, not a new decision - TCGA's sets went through exactly the same map. Not
running it would be the deviation.

Its properties are what make it safe:

- the source is **MitoCarta's own curated `Synonyms` column** - not a guess, and
  not a new annotation package introduced mid-arm;
- it runs **forward only**, current symbol -> its listed alias. Script 07
  documents why the reverse is dangerous: `COX1`/`COX2`/`COX3` resolve to the
  prostaglandin synthases, injecting two abundant inflammatory genes into an
  OXPHOS set, invisibly;
- a candidate that is itself a current MitoCarta symbol for a different gene is
  **refused**, and so is any symbol with more than one surviving candidate.
  Ambiguity is left unresolved and reported by name.

Verified on this matrix: **all 19 resolve, none ambiguous, no two map to the
same row, coverage 0.787 -> 1.000.** Script 16 stops if two inputs ever collide
onto one row, because that would double-weight a gene inside a pathway mean with
nothing visible downstream.

## A.3 What is deliberately NOT rescued

> **CORRECTED 2026-08-31, after script 16 ran. The premise below is wrong.**
> I wrote that MitoCarta's synonyms "cannot resolve" the Felsher and Hallmark
> renames because those genes are not MitoCarta genes. They can, and they did.
> **MitoCarta 3.0's sheet 3 is the whole 19,247-gene human background**, not the
> 1,136 mitochondrial genes, and it carries a `Synonyms` column for every gene
> in it. The map is therefore a general HGNC alias source already, and the
> section below argues against introducing something the arm has had since
> 2026-08-28.
>
> What actually happened: **27 symbols resolved, not 19**, all of them genuine
> historical renames.
>
> ```
> EEF1AKNMT -> METTL13    H2AZ2  -> H2AFV      TENT4A -> PAPD7
> H2AX      -> H2AFX      H2BC12 -> HIST1H2BK  VARS1  -> VARS
> H2AZ1     -> H2AFZ      POLR1G -> CD3EAP
> ```
>
> Final coverage, after harmonisation and the low-count filter:
>
> ```
>                            raw     harmonised   still missing
> OXPHOS subunits           0.775      0.989      COX8C
> PROLIF_DISJOINT           0.978      0.994      PRP4K, PTTG3P
> PROLIF_STD                0.979      0.994      PRP4K, PTTG3P
> Felsher M-a (61)          0.951      1.000      -
> ```
>
> `COX8C` is testis-specific and `PTTG3P` is a pseudogene; neither is expected
> in breast at the filter's threshold. Three unresolved, none ambiguous.
>
> **The discipline argument survives intact, and is stronger than the version
> below.** The map was fixed on 2026-08-28 in script 07 section 2, built over
> the union of all four sets at once, and applied uniformly. Nothing was chosen
> after seeing which genes were missing - which is what the paragraph below was
> worried about. The asymmetry it anticipated (exposure harmonised, covariates
> not) **did not occur**: all four sets went through the same map.
>
> **One concrete instance of why the guard matters**, worth recording because it
> is the failure this map is built to prevent. `POLR1G`'s MitoCarta synonyms are
> `ASE-1 | ASE1 | CAST | CD3EAP | PAF49 | RPA34`. Two of those are in the SCAN-B
> matrix: `CD3EAP`, the genuine former symbol, and **`CAST`, which is
> calpastatin - a real and entirely different gene**. The rule that refuses a
> candidate which is itself a current symbol removed `CAST` and left `CD3EAP`
> unique. Without it the MYC estimator would have silently gained calpastatin
> and lost `POLR1G`.
>
> The section below is kept as written, because the reasoning in its last
> paragraph is still the right reasoning for any gene the map genuinely cannot
> reach - there were three, and they were left alone.



The Felsher estimator and the Hallmark proliferation sets are not MitoCarta
genes, so MitoCarta's synonyms cannot resolve their renames:

```
Felsher M-a (61)     0.951    missing POLR1G, VARS1, EEF1AKNMT
PROLIF_DISJOINT      0.978    missing H2AX, H2AZ1, H2AZ2, H2BC12, PRP4K,
PROLIF_STD           0.979            PTTG3P, TENT4A
```

All are above the 0.80 floor. They are **reported by name and left unresolved**.
A general alias source (org.Hs.eg.db, limma's alias tables) would resolve them,
and introducing one for this cohort only - chosen after seeing which genes were
missing - is precisely the move the coverage floor exists to make unnecessary.

This asymmetry is stated plainly because it is real: the exposure is harmonised
and the covariates are not. The justification is that the harmonisation source
was fixed by the arm in advance for MitoCarta sets and does not exist for the
others, not that one mattered more than the other.

## A.4 What this changes in the model

Nothing. The gene sets are the same genes; only the strings differ. The
declared model, the covariate mapping, the specifications, the null, the
replication criteria and the failure branch are all unchanged.

`symbol_map` and `symbol_report` are saved in `results/scanb_pheno.rds`, and
**script 17 must apply `symbol_map` before scoring**. Coverage is reported both
before and after harmonisation so the raw number is never lost.
