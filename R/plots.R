# plots.R
# Figures in APA 7 style: no axis tick marks, full range on bounded axes,
# colorblind-safe colors with a second cue (line type or shape), no text over
# data. Figure numbers and titles go in the caption, outside the image.

okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000")

# Proportions: no leading zero (APA 7), and short end labels so adjacent
# panels' axis labels can't run together.
prop_breaks <- seq(0, 1, 0.25)
prop_labels <- c("0", ".25", ".50", ".75", "1")

theme_apa_fig <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base_size),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.spacing = grid::unit(1.2, "lines"),
          plot.margin = margin(8, 12, 8, 8))
}

save_fig <- function(p, name, width = 7, height = 4.5) {
  ggsave(out_path("figures", paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(out_path("figures", paste0(name, ".pdf")), p, width = width, height = height)
  invisible(p)
}

# Latent mean (eta) for every cell from a fitted unequal-variance model.
cell_eta <- function(fit, grid = cell_grid()) {
  grid$eta <- drop(cell_X(grid) %*% stats::coef(fit)[colnames(cell_X(grid))])
  grid
}

# Figure: ROC curves for tactic-laden requests on the argument screen,
# control vs trained, under the cue and the vigilance accounts.
plot_roc <- function(fits, palter_share = 0.5) {
  curves <- points <- list()
  for (h in names(fits)) {
    fit <- fits[[h]]
    # Unsafe requests are a mix of paltering and non-paltering ones, so the hit
    # rate averages over the two.
    pick <- function(p) {
      g <- cell_eta(fit, transform(cell_grid(), palter = unsafe * p))
      g[g$tactic == 1 & g$screen == "argument", ]
    }
    g0 <- pick(0); g1 <- pick(1)
    sig <- exp(fit$zeta)
    hit_at <- function(th, tr) {
      e0 <- g0$eta[g0$unsafe == 1 & g0$trained == tr]; e1 <- g1$eta[g1$unsafe == 1 & g1$trained == tr]
      (1 - palter_share) * (1 - stats::pnorm((th - e0) / sig)) + palter_share * (1 - stats::pnorm((th - e1) / sig))
    }
    for (tr in 0:1) {
      es <- g0$eta[g0$unsafe == 0 & g0$trained == tr]
      lab <- if (tr == 1) "Trained" else "Control"
      th <- seq(-6, 8, length.out = 400)
      curves[[length(curves) + 1]] <- data.frame(
        account = h, group = lab, fa = 1 - stats::pnorm(th - es), hit = hit_at(th, tr))
      points[[length(points) + 1]] <- data.frame(
        account = h, group = lab, fa = 1 - stats::pnorm(fit$alpha - es), hit = hit_at(fit$alpha, tr),
        cut = seq_along(fit$alpha))
    }
  }
  curves <- do.call(rbind, curves); points <- do.call(rbind, points)
  lab_acc <- c(cue = "Cue account", vigilance = "Vigilance account")
  curves$account <- factor(lab_acc[curves$account], levels = lab_acc)
  points$account <- factor(lab_acc[points$account], levels = lab_acc)

  ggplot(curves, aes(fa, hit, colour = group, linetype = group)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.4) +
    geom_path(linewidth = 0.8) +
    geom_point(data = points, aes(shape = group), size = 2.2) +
    facet_wrap(~ account) +
    scale_colour_manual(values = okabe_ito[1:2]) +
    scale_shape_manual(values = c(16, 17)) +
    scale_x_continuous("False alarms: share of safe requests rejected",
                       limits = c(0, 1), breaks = prop_breaks, labels = prop_labels, expand = c(0, 0)) +
    scale_y_continuous("Hits: share of unsafe requests rejected",
                       limits = c(0, 1), breaks = prop_breaks, labels = prop_labels, expand = c(0, 0)) +
    coord_equal() +
    theme_apa_fig() + theme(panel.spacing = grid::unit(2.5, "lines"))
}

# Figure: training's effect on d for tactic-laden requests when only the
# criterion moved (cue account), estimated two ways.
# bands: one row per model with that model's own SESOI.
plot_ev_vs_uv <- function(est, bands) {
  ggplot(est, aes(est)) +
    geom_rect(data = bands, aes(xmin = -sesoi, xmax = sesoi, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey92") +
    geom_vline(data = bands, aes(xintercept = -sesoi), linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_vline(data = bands, aes(xintercept = sesoi), linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_histogram(binwidth = 0.02, fill = okabe_ito[1], colour = "white", boundary = 0) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    facet_wrap(~ model, ncol = 1) +
    scale_x_continuous(expression(paste("Estimated change in ", italic(d), " from training, tactic-laden requests")),
                       limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
    scale_y_continuous("Simulated studies", expand = expansion(mult = c(0, 0.05))) +
    theme_apa_fig()
}

# Figure: how often the preregistered rule reaches the right verdict.
plot_adjudication <- function(rates, confusion, n_show) {
  lab_h <- c(cue = "Cue", vigilance = "Vigilance", both = "Both", null = "Neither")
  rates$truth <- factor(lab_h[rates$truth], levels = lab_h)
  dodge <- position_dodge(width = 6)  # keeps overlapping series visible
  a <- ggplot(rates, aes(n_per_cell, p_correct, colour = truth, shape = truth, linetype = truth)) +
    geom_hline(yintercept = 0.8, colour = "grey70", linewidth = 0.4) +
    geom_line(linewidth = 0.7, position = dodge) + geom_point(size = 2, position = dodge) +
    scale_colour_manual(values = okabe_ito[1:4]) +
    scale_y_continuous("Share of studies with the correct verdict",
                       limits = c(0, 1), breaks = prop_breaks, labels = prop_labels, expand = c(0, 0.01)) +
    scale_x_continuous("Participants per training \u00d7 screen cell",
                       breaks = sort(unique(rates$n_per_cell))) +
    theme_apa_fig() + theme(legend.position = "right")

  confusion$truth <- factor(lab_h[confusion$truth], levels = rev(lab_h))
  lab_v <- c(lab_h, inconclusive = "Inconclusive")
  confusion$verdict <- factor(lab_v[confusion$verdict], levels = c(lab_h, "Inconclusive"))
  b <- ggplot(confusion, aes(verdict, truth, fill = share)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sub("^0", "", sprintf("%.2f", share)),
                  colour = ifelse(share > 0.5, "light", "dark")), size = 3.4) +
    scale_colour_manual(values = c(light = "white", dark = "black"), guide = "none") +
    scale_fill_gradient(low = "grey95", high = okabe_ito[1], limits = c(0, 1), guide = "none") +
    labs(x = sprintf("Verdict (%d per cell)", n_show), y = "True state") +
    theme_apa_fig() + theme(axis.line = element_blank())
  a / b + plot_annotation(tag_levels = "A")
}

# Figure: simulated reviewers with and without a human correction.
r_label <- expression(paste("Accuracy of the profile used to simulate each person (", italic(r), ")"))

plot_ppi <- function(summ) {
  summ$method <- factor(summ$method, levels = c("People only", "Simulator only", "PPI", "PPI++"))
  lt <- c("People only" = "solid", "Simulator only" = "dashed", "PPI" = "dotted", "PPI++" = "dotdash")
  sh <- c("People only" = 16, "Simulator only" = 17, "PPI" = 15, "PPI++" = 3)
  cl <- c("People only" = okabe_ito[1], "Simulator only" = okabe_ito[2],
          "PPI" = okabe_ito[3], "PPI++" = okabe_ito[4])
  dodge <- position_dodge(width = 0.06)
  a <- ggplot(summ, aes(r_profile, coverage, colour = method, shape = method, linetype = method)) +
    geom_hline(yintercept = 0.95, colour = "grey70", linewidth = 0.4) +
    geom_line(linewidth = 0.7, position = dodge) + geom_point(size = 2, position = dodge) +
    facet_wrap(~ target_label) +
    scale_colour_manual(values = cl) + scale_shape_manual(values = sh) +
    scale_linetype_manual(values = lt) +
    scale_y_continuous("Coverage of 95% intervals", breaks = prop_breaks, labels = prop_labels) +
    scale_x_continuous(r_label, breaks = prop_breaks, labels = prop_labels) +
    coord_cartesian(xlim = c(-0.06, 1.06), ylim = c(0, 1)) +  # zoom, so dodged points aren't dropped
    theme_apa_fig()
  b <- ggplot(summ[summ$method %in% c("PPI", "PPI++"), ],
              aes(r_profile, gain_pct, colour = method, shape = method, linetype = method)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +
    facet_wrap(~ target_label) +
    scale_colour_manual(values = cl) + scale_shape_manual(values = sh) +
    scale_linetype_manual(values = lt) +
    scale_y_continuous("Gain in effective sample size (%)") +
    expand_limits(y = c(0, 100)) +
    scale_x_continuous(r_label, limits = c(0, 1), breaks = prop_breaks, labels = prop_labels) +
    theme_apa_fig()
  a / b + plot_annotation(tag_levels = "A")
}

# Figure: power of the tell audit by bank size.
plot_audit <- function(pw) {
  pw$bank <- factor(ifelse(pw$leak > 0, "Unsafe requests minimize more often", "No tell"),
                    levels = c("No tell", "Unsafe requests minimize more often"))
  # Wilson 95% intervals: each point rests on only pw$audits audits.
  z <- stats::qnorm(0.975); k <- pw$audits
  centre <- (pw$p_flag + z^2 / (2 * k)) / (1 + z^2 / k)
  half <- z * sqrt(pw$p_flag * (1 - pw$p_flag) / k + z^2 / (4 * k^2)) / (1 + z^2 / k)
  pw$lo <- centre - half; pw$hi <- centre + half
  dodge <- position_dodge(width = 8)
  x_axis <- scale_x_continuous("Requests in the bank (25% unsafe)", breaks = sort(unique(pw$n_items)))
  a <- ggplot(pw, aes(n_items, p_flag, colour = bank, shape = bank, linetype = bank)) +
    geom_hline(yintercept = 0.05, colour = "grey70", linewidth = 0.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linetype = "solid", position = dodge) +
    geom_line(linewidth = 0.7, position = dodge) + geom_point(size = 2.2, position = dodge) +
    scale_colour_manual(values = okabe_ito[c(1, 2)]) +
    scale_y_continuous("Share of audits that flag a tell",
                       limits = c(0, 1), breaks = prop_breaks, labels = prop_labels) +
    x_axis + theme_apa_fig()
  if (is.null(pw$auc_median)) return(a)
  # Panel B: the classifier's cross-validated AUC across the same audits
  # (median and interquartile range); .5 is chance.
  b <- ggplot(pw, aes(n_items, auc_median, colour = bank, shape = bank, linetype = bank)) +
    geom_hline(yintercept = 0.5, colour = "grey70", linewidth = 0.4) +
    geom_errorbar(aes(ymin = auc_q25, ymax = auc_q75), width = 0, linetype = "solid", position = dodge) +
    geom_line(linewidth = 0.7, position = dodge) + geom_point(size = 2.2, position = dodge) +
    scale_colour_manual(values = okabe_ito[c(1, 2)]) +
    scale_y_continuous("Cross-validated AUC",
                       limits = c(0, 1), breaks = prop_breaks, labels = prop_labels) +
    x_axis + theme_apa_fig()
  a / b + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")
}
