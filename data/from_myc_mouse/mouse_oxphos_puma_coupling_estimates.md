# Mouse MYC x OXPHOS coupling to the PUMA/Bcl-xL ratio - transcribed estimates

Companion note to `mouse_oxphos_puma_coupling_estimates.csv`. Read this before quoting
any number from that file.

## Why this is a transcription and not a copied table

The producing script, `scripts/43_substrate_specificity_and_tradeoff.R` in `myc_mouse`,
writes **no CSV**. Its only output is

    saveRDS(out, results/substrate_specificity_tradeoff.rds)   # 5.3 KB

holding `$tradeoff`, `$tradeoff_perm`, `$ambient`. The ruler-robustness rerun,
`scripts/46_axis_ruler_test.R`, likewise writes only `results/axis_ruler_test.rds`.

Neither object was copied into this repo: they are result objects, not summary tables,
and extracting from them would have required running R. The values below were instead
transcribed by hand from the two mouse-repo markdown tables that already report them.
Every row in the CSV carries its `source_doc` and `source_line`.

## The headline number, with the caveats that must travel with it

The manuscript sentence quotes **p = 0.0052** for the MYC x OXPHOS interaction on the
PUMA:Bcl-xL ratio. Three qualifications are recorded in the mouse repo alongside it, and
a bare "p = 0.0052" in the human arm would misrepresent the mouse result:

1. **The permutation null is not significant.** Script 43's own within-timepoint
   permutation null puts the effect at the 91.7th percentile of 5,000 draws,
   **empirical p = 0.083** (`PANELS.md:1226`). Script 46's rerun gives **0.080**
   (`instrument_choice:250`). The two differ only by permutation draw.
2. **It does not survive a timepoint term.** Adding one moves p from **0.0052 to 0.088**
   (`PANELS.md:1226`).
3. **The specificity is ruler-dependent.** The direction survives swapping mitoPPS for
   absolute OXPHOS levels (beta +0.535 against +0.779, the two rulers correlate
   r = +0.928), but the significance does not (p 0.156). Only **7.1%** of the absolute
   level's variance is separable from the global expression factor (r = +0.964) against
   **55.0%** for mitoPPS. On the absolute ruler the axis is "overall expression", not
   "the respiratory chain" (`instrument_choice:250-263`).

The mouse repo's own conclusion: **keep mitoPPS**, and report the effect as a lead
rather than as an established interaction.

## Two further cautions

- **The two betas are on different scalings.** Script 43's interaction is reported as
  **+6.09** and script 46's reference beta as **+0.779**, both at p = 0.0052. Script 46
  states it reproduces script 43's model exactly (max |delta beta| = 0 over 10 fits), so
  these are the same fit expressed in different units - script 43's panel draws the raw
  log2 ratio, the ruler test a z-scored outcome. The exact relationship was **not
  verified** in the snapshot session that wrote this note. Do not treat +6.09 and +0.779
  as independent estimates, and do not convert between them without checking the source.
- **The control battery is two axes, not many.** Script 43 carries **exactly two**,
  OXPHOS subunits and redox. The draft sentence "no such link exists for other redox and
  metabolic axes" overstates the panel that was run (`PANELS.md:1230`). Note that redox
  is not a null axis - both genotypes couple to it strongly (wild-type slope -7.70);
  what it does not do is couple *differently*. That is what makes it the right control.

## Scope

These are **mouse** estimates, frozen for cross-species comparison. They are a reference
point for what the human arm is testing, never an input to a human model.
