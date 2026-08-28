---
date: 2026-08-28
status: RESOLVED 2026-08-28 - all 8 items agreed, applied to plan and CLAUDE.md
relates-to:
  - 2026-08-27_human_validation_plan.md (sections 1, 2, 7.2)
  - 2026-08-28_D4_scoping_and_G2_design.md (the FAO caveat)
  - CLAUDE.md (gene-set rules, OXPHOS umbrella vs subunits)
decides:
  - script 07 specificity panel source and membership
  - RESOLVED: BOTH instruments co-primary; claim only what both support
supersedes: the "build it from MitoCarta pathways, or take the library GMT"
  framing. Both were wrong about what the mouse actually did.
---

# Specificity panel for script 07 - proposal

Written after reading the mouse figure panels and `scripts/43_substrate_
specificity_and_tradeoff.R` (read-only, via `git show`). Everything below is
sourced from what the mouse arm actually ran, not from the plan's paraphrase of
it.

---

## 1. The source question is settled: Human MitoCarta 3.0

The mouse specificity arms are **MitoCarta pathways**, named directly in script
43:

```
OXPHOS subunits          MITOCARTA_OXPHOS_SUBUNITS
OXPHOS (all)             MITOCARTA_OXPHOS
OXPHOS assembly          MITOCARTA_OXPHOS_ASSEMBLY_FACTORS
nucleotide metabolism    MITOCARTA_NUCLEOTIDE_METABOLISM
mitoribosome             MITOCARTA_MITOCHONDRIAL_RIBOSOME
TCA cycle                MITOCARTA_TCA_CYCLE
amino-acid metabolism    MITOCARTA_AMINO_ACID_METABOLISM
lipid metabolism         MITOCARTA_LIPID_METABOLISM
TEB vs ductal (HS)       MG_TEB_VS_DUCTAL_HS_GRAY_UP        (mouse-only)
PROLIF_* pooled          (union of the library's PROLIF_ sets)
```

So the human panel takes **the same MitoPathway names against Human MitoCarta
3.0**, which is already snapshotted in `data/mitocarta_human/`. No library GMT,
no rebuild, no ortholog projection - the cross-species link is the *pathway
name*, and each species uses its own native MitoCarta. That is the cleanest
possible route and it was available all along.

All arms resolve, verified 2026-08-28:

| MitoPathway | human genes | mouse role |
|---|---|---|
| `OXPHOS subunits` | **102** | the claim |
| `OXPHOS` (umbrella) | 169 | reported, not primary |
| `OXPHOS assembly factors` | **68** | **the internal control** |
| `Mitochondrial ribosome` | 83 | growth alternative |
| `Nucleotide metabolism` | 40 | growth alternative |
| `TCA cycle` | 20 | metabolic alternative |
| `Amino acid metabolism` | 90 | metabolic alternative |
| `Lipid metabolism` | 123 | metabolic alternative |
| `ROS and glutathione metabolism` | 27 | **the mouse's own control axis** |
| `Fatty acid oxidation` | 44 | plan negative |
| `Carnitine shuttle` | 5 | plan negative - too small, see 4 |
| `Folate and 1-C metabolism` | 21 | plan negative |
| `Glycine metabolism` | 15 | plan negative |

## 2. The plan is missing the mouse's tightest control

Plan section 2 lists the pathway negatives as: FAO/carnitine shuttle,
one-carbon/glycine cleavage, mitoribosome, TCA, ROS defence, mtDNA-encoded.

**It does not include OXPHOS assembly factors.** That is the mouse's sharpest
specificity result, stated in `fig2_wt_mito_contraction.R`:

> *THE ASSEMBLY FACTORS OF THE SAME COMPLEXES SIT AT THE ORIGIN, a few
> millimetres from their own subunits. That is the internal control, and the
> distance between the two points is the whole specificity claim.*

> *the OXPHOS subunits sit at percentile 0.0 of 2000 expression-matched sets
> while their assembly factors sit at 50.2*

Subunits and assembly factors are the **same complexes** - same organelle, same
pathway umbrella, same transcriptional neighbourhood, different function. A
distant pathway like one-carbon metabolism cannot exclude what that pair
excludes. CLAUDE.md already flags the umbrella-vs-subunits distinction as a
deliberate choice; this is the reason it matters.

**Proposal: `OXPHOS assembly factors` becomes the PRIMARY negative**, with the
plan's list retained as the secondary battery.

## 3. Comparators should represent the alternative, not just be other pathways

Script 43 states its selection logic explicitly:

> *The comparator arms are chosen so that the alternative hypothesis has its best
> possible representatives: NUCLEOTIDE metabolism and the MITORIBOSOME are the
> two mitochondrial arms most tightly tied to growth... If the wild-type
> respiratory withdrawal were a growth withdrawal, those three should move with
> OXPHOS.*

That is a better design principle than "hold other axes separately". Each arm
should answer a named alternative. Mapping the human battery that way:

| Alternative to H1 | Best representative | Arm |
|---|---|---|
| it is mitochondrial biogenesis, not respiration | same complexes, non-catalytic | `OXPHOS assembly factors` |
| it is growth / translation capacity | mito translation apparatus | `Mitochondrial ribosome` |
| it is proliferation | nucleotide supply, plus the D7 proliferation score | `Nucleotide metabolism` + Hallmark E2F/G2M |
| it is redox, not respiration | the mouse's own control axis | `ROS and glutathione metabolism` |
| it is generic mitochondrial metabolism | the three big catabolic arms | `TCA`, `Amino acid`, `Lipid` |
| it is FAO (Lee et al. 2017) | see 4 | `Fatty acid oxidation` |
| it is a scale artefact | expression-matched random sets | the null, see 5 |

## 4. Two adjustments the human data forces

**`Carnitine shuttle` is 5 genes.** Below any sensible floor for GSVA. Merge it
into `Fatty acid oxidation` (44) and report the union as the FAO arm, or drop it
and say so. Do not score a 5-gene set and present it as a negative.

**FAO carries the D4 caveat and needs its own sentence.** Lee et al. 2017 report
that TNBC cancer stem cells are preferentially FAO-dependent (etomoxir suppresses
OCR more than BPTES or UK5099). So if the FAO arm fires, the honest reading is
"consistent with Lee et al." rather than "the OXPHOS coupling is non-specific".
That must be pre-stated, not decided after seeing it.

The same applies more weakly to one-carbon via `MTHFD2`, already recorded in the
G1 note section 7.

## 5. The null design: expression-matched random sets

The mouse does not compare arms to each other and stop. It compares the observed
arm against **2000 expression-matched random gene sets**, and that is what makes
the result interpretable:

- OXPHOS subunits at percentile **0.0** of 2000, assembly factors at **50.2**
- the within-OXPHOS expression gradient sits at the **100th** percentile against
  a null whose median has the **opposite sign** (-0.17, -0.22)

Without an expression-matched null, "OXPHOS moved and TCA did not" is confounded
by set size and expression level. **Proposal: adopt the same 2000
expression-matched null for every arm in the human battery.** It is cheap and it
is the difference between a ranking and a test.

## 6. RESOLVED - BOTH instruments are co-primary

This is the one I cannot settle, and it is not cosmetic. `figS2_oxphos_subunit_
heatmap.R` shows the same 87 mouse subunits on the same contrast giving:

```
unweighted per-gene mean log2FC   +0.061    "returns to baseline"
expression-WEIGHTED mean          +0.201
log2(sum of normalised counts)    +0.226    the compartment-level ruler
```

> *Both are correct. They weight differently.*

The sign of the headline number depends on the weighting, because the
respiratory chain's response is **expression-graded** - abundant subunits move
up, low-expressed ones move down, and that gradient sits at the 100th percentile
against expression-matched nulls.

Three candidate human instruments, and they are not interchangeable:

| Instrument | What it measures | Mouse analogue |
|---|---|---|
| **GSVA on VST** (plan's primary) | rank-based, roughly unweighted | closest to `+0.061` |
| **mitoPPS** (linear, composition) | budget share, content-blind | **the axis the mouse interaction was fitted on** |
| summed/level metric | compartment content | closest to `+0.226` |

**The mouse interaction (p = 0.0052) was fitted on mitoPPS**, not on a GSVA
level. The plan's primary human OXPHOS measure is a GSVA level. So as written,
the human arm would test a different quantity from the one the mouse result
describes, and a null could mean nothing more than an instrument mismatch.

**DECISION 2026-08-28: option 3. Both co-primary, reported side by side.**

> **Rule for disagreement, fixed before any model is fitted:** report both;
> claim only what **both** instruments support. An effect appearing on one
> instrument alone is reported as instrument-dependent and is NOT a positive
> result.

Consequence for Block F: GSVA is the instrument that travels, because mitoPPS is
composition-dependent and cannot be compared across cohorts. So the
meta-analysis runs on GSVA and mitoPPS is TCGA-internal. That asymmetry is a
property of the instruments, not a choice, and it must be stated rather than
quietly worked around.

Options as considered, in the order I ranked them:

1. **mitoPPS primary, GSVA secondary** - matches the mouse instrument, so the
   cross-species claim is like-for-like. Costs: mitoPPS is composition-dependent
   and cohort-relative (CLAUDE.md), so it cannot be compared across cohorts at
   all - only the pattern transfers. That constrains Block F.
2. **GSVA primary, mitoPPS secondary** - the plan as written. Portable across
   cohorts, but not the quantity the mouse measured.
3. **Both as co-primary, reported side by side** - honest, costs a doubled
   battery, and needs a rule for what to conclude when they disagree.

Taken: **3**. The whole arm exists to test a specific mouse claim; running an
instrument the mouse did not use, and calling a null a failure to replicate,
would not be defensible.

## 7. Where `GS_metabolic` fits

`data/raw/user_curated/GS_metabolic_genes_list.xlsx` in the geneset library is
tracked and tag-pinned (blob `21f79d9a...` at `v1.0`) and has a **native `human`
sheet** - 2,348 genes, 74 classifications, human symbols with Entrez ids. The
library maps it *to* mouse, so the human version is the original.

It is **not** needed for the specificity battery, which is MitoCarta. It earns
its place for one thing MitoCarta does not give: **complex-level resolution**.

```
Complex I 49 | Complex II 4 | Complex III 11 | Complex IV 27 | ATPase 5
Proton Transport 49 | Ubiquinone 13
```

The mouse claim is that OXPHOS subunit LFC and mitoPPS drop **across all
complexes**. Testing that in human needs per-complex sets, and this sheet has
them where MitoCarta's `OXPHOS subunits` is one flat list of 102. Proposal:
snapshot the human sheet for that purpose only, with its own provenance README,
and keep it out of the primary battery.

**This also means CLAUDE.md's blanket rejection of library assets is too broad.**
The rejection is correct for `outputs/gmt/human/`, which is mouse-native sets run
through `mouse_to_human()`, gitignored and unpinned. It should not extend to a
tracked raw input with a native human sheet. Those are different objects.
Proposed wording change in section 8 below.

## 8. Proposed CLAUDE.md amendment

Current text rejects the library's human GMT tree. Add, immediately after it:

> **The rejection is of `outputs/gmt/human/`, not of the library.** That tree is
> mouse-native sets pushed through `mouse_to_human()`, gitignored and unpinned by
> the tag. Tracked *raw inputs* that carry their own native human data are a
> different matter: `data/raw/user_curated/GS_metabolic_genes_list.xlsx` has a
> `human` sheet (2,348 genes, 74 classifications) which the library maps *to*
> mouse. Taking the human sheet is not a round trip. Check the sheet, not the
> repository.

## 9. All eight agreed 2026-08-28

| # | Item | Applied to |
|---|---|---|
| 1 | MitoCarta human as the battery source | CLAUDE.md, plan section 2 |
| 2 | `OXPHOS assembly factors` -> PRIMARY negative | plan section 2, CLAUDE.md |
| 3 | Comparators framed as named alternatives | plan section 2 |
| 4 | `Carnitine shuttle` merged into FAO; FAO/Lee caveat pre-stated | plan section 2 |
| 5 | 2,000 expression-matched random sets as the null | plan section 2 |
| 6 | **Both instruments co-primary**; claim only what both support | plan section 7.2, CLAUDE.md |
| 7 | `GS_metabolic` human sheet snapshotted, complex-level only | `data/genesets_metabolic_human/` |
| 8 | CLAUDE.md amendment - the rejection is of `outputs/gmt/human/`, not the library | CLAUDE.md |

Nothing here is open. Script 07 can be written.
