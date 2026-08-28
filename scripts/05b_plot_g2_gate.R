# 05b_plot_g2_gate.R
# =============================================================================
# Diagnostic figure for gate G2.  Reads results/g2_cnv_cooccurrence.rds and
# writes one four-panel PDF/PNG.  No computation happens here - if a number
# looks wrong, fix script 05, not this file.
#
# Separate from 05 so it can be re-run in a second without re-reading the
# 536 MB ISAR file.  Numbered 05b rather than folded into 18_figures.R because
# this is a gate diagnostic, not a manuscript panel; plan section 11 reserves 18
# for display items.
#
# SCALE DISCIPLINE: not applicable, no expression data.  See script 05 header.
# =============================================================================

source(here::here("scripts", "00_setup_packages.R"))

g2 <- readRDS(file.path(DIR_RESULTS, "g2_cnv_cooccurrence.rds"))

DIR_FIGS <- here::here("outputs", "figures")
.ensure_dir(DIR_FIGS)

# -----------------------------------------------------------------------------
# 1. Drop the degenerate fits
# -----------------------------------------------------------------------------
# The "le_neg2" rows (homozygous deletion) have a joint count of ZERO with MYC
# amplification, so glm separates: the odds ratio comes back around 1e-8 with an
# upper CI around 1e88.  Those are not estimates and they would destroy a log
# scale.  They are dropped here and the count is printed, so the omission is
# visible rather than silent.  The same happens in any stratum where n_both = 0
# (BRCA_Normal, mostly).
.finite_or <- function(d) {
  keep <- is.finite(d$adj_or) & is.finite(d$adj_lo) & is.finite(d$adj_hi) &
          d$adj_lo > 1e-3 & d$adj_hi < 1e3
  dropped <- sum(!keep)
  if (dropped > 0) {
    message("   dropped ", dropped, " degenerate fit(s) (separation, n_both = 0)")
  }
  d[keep, , drop = FALSE]
}

.lab_rule <- c(eq2 = "amp (+2)", ge1 = "gain+ (>=+1)",
               le_neg1 = "loss (<=-1)", le_neg2 = "homdel (-2)")

# -----------------------------------------------------------------------------
# 2. Assemble the four panels
# -----------------------------------------------------------------------------
message("1. assembling panels")

p1 <- g2$cooccurrence %>%
  dplyr::filter(tier %in% c("primary", "regional_control")) %>%
  dplyr::transmute(
    panel = "a. Primary tests (ISAR, aneuploidy-adjusted)",
    label = paste0(partner, "  [", .lab_rule[partner_rule], "]"),
    adj_or, adj_lo, adj_hi, n_both,
    flag  = ifelse(tier == "regional_control", "19q13 regional control",
                   ifelse(adj_lo > 1, "passes", "fails"))
  ) %>%
  .finite_or()

p2 <- g2$cooccurrence %>%
  dplyr::filter(tier == "secondary_grid", exposure_rule == "eq2") %>%
  dplyr::transmute(
    panel = "b. Threshold sensitivity (MYC at +2 throughout)",
    label = paste0(partner, "  [", .lab_rule[partner_rule], "]"),
    adj_or, adj_lo, adj_hi, n_both,
    flag  = ifelse(adj_lo > 1, "passes", ifelse(adj_hi < 1, "inverts", "fails"))
  ) %>%
  .finite_or()

p3 <- g2$stratified %>%
  dplyr::filter(strat_var == "pam50", !skipped) %>%
  dplyr::transmute(
    panel = "c. PAM50 strata (primary rule per gene)",
    label = paste0(partner, "  ", sub("^BRCA_", "", stratum)),
    adj_or, adj_lo, adj_hi, n_both,
    flag  = ifelse(adj_lo > 1, "passes", ifelse(adj_hi < 1, "inverts", "fails"))
  ) %>%
  .finite_or()

p4 <- g2$cooccurrence %>%
  dplyr::filter(tier %in% c("primary", "sensitivity_source", "sensitivity_fga")) %>%
  dplyr::transmute(
    panel = "d. Source and covariate sensitivity",
    label = paste0(
      partner, "  ",
      dplyr::case_when(
        tier == "primary"            ~ "A ISAR / aneu",
        tier == "sensitivity_fga"    ~ "A ISAR / FGA",
        source == "B_Firehose"       ~ "B Firehose",
        TRUE                         ~ "B intersection"
      )
    ),
    adj_or, adj_lo, adj_hi, n_both,
    flag = ifelse(adj_lo > 1, "passes", ifelse(adj_hi < 1, "inverts", "fails"))
  ) %>%
  .finite_or()

plot_df <- dplyr::bind_rows(p1, p2, p3, p4) %>%
  dplyr::mutate(
    panel = factor(panel, levels = unique(c(p1$panel, p2$panel, p3$panel, p4$panel))),
    label = factor(label, levels = rev(unique(label))),
    flag  = factor(flag, levels = c("passes", "fails", "inverts",
                                    "19q13 regional control"))
  )

# -----------------------------------------------------------------------------
# 3. Forest plot
# -----------------------------------------------------------------------------
message("2. drawing")

pal <- c("passes" = "#1b7837", "fails" = "#999999", "inverts" = "#762a83",
         "19q13 regional control" = "#d95f02")

gg <- ggplot2::ggplot(plot_df,
                      ggplot2::aes(x = adj_or, y = label, colour = flag)) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = adj_lo, xmax = adj_hi),
                          height = 0, linewidth = 0.6) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::geom_text(ggplot2::aes(label = paste0("n=", n_both)),
                     hjust = -0.25, vjust = -0.9, size = 2.6,
                     colour = "grey30", show.legend = FALSE) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_colour_manual(values = pal, drop = FALSE, name = NULL) +
  ggplot2::facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  ggplot2::labs(
    title    = "Gate G2 - copy-number co-occurrence with MYC amplification in TCGA-BRCA",
    subtitle = paste0(
      "Aneuploidy-adjusted odds ratios, 95% CI, log scale. ",
      "Pass = CI excludes 1 (design note 2.8).\n",
      "Thresholds are gene- and direction-specific; ORs are NOT comparable ",
      "across genes at different thresholds."
    ),
    x = "adjusted odds ratio (log scale)", y = NULL,
    caption = paste0("source: results/g2_cnv_cooccurrence.rds, generated ",
                     format(g2$generated, "%Y-%m-%d"),
                     " | degenerate (separated) fits omitted")
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    legend.position = "bottom",
    strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
    strip.text = ggplot2::element_text(face = "bold", hjust = 0),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(file.path(DIR_FIGS, "g2_gate_forest.pdf"), gg,
                width = 10, height = 8.5)
ggplot2::ggsave(file.path(DIR_FIGS, "g2_gate_forest.png"), gg,
                width = 10, height = 8.5, dpi = 200)

message("05b: done. outputs/figures/g2_gate_forest.{pdf,png}")

# =============================================================================
# Sandbox - skipped by source(), run line by line in Positron
# =============================================================================
if (FALSE) {

  gg   # draw to the Positron plot pane

  # What got dropped as degenerate, and why. Every one of these should have
  # n_both == 0; if any does not, something else is wrong.
  dplyr::bind_rows(
    dplyr::filter(g2$cooccurrence, tier == "secondary_grid"),
    dplyr::filter(g2$stratified, !skipped)
  ) %>%
    dplyr::filter(!is.finite(adj_lo) | adj_hi > 1e3) %>%
    dplyr::select(partner, partner_rule, dplyr::any_of("stratum"),
                  n_both, adj_or, adj_lo, adj_hi) %>%
    print(n = 30)

  # The single most important comparison in the figure: is BBC3 separable from
  # its 19q13 regional control?
  g2$cooccurrence %>%
    dplyr::filter(partner %in% c("BBC3", "BAX"), tier != "secondary_grid") %>%
    dplyr::select(source, tier, partner, adj_or, adj_lo, adj_hi, adj_p) %>%
    print(n = 20)

}
