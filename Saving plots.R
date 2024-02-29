# load in the data and get plotting functions
source("Tools and data.R")

library("qpcR") # for AIC comparison
library("lmtest") # for Granger causality
library("ggh4x") # plotting stuff
library("ggpubr") # combining cilia and cyano plots

# from here https://ggplot2.tidyverse.org/reference/labeller.html
capitalise <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

tag_facet <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf,
    y = Inf, hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
    gb <- ggplot_build(p)
    lay <- gb$layout$layout
    tags <- cbind(
        lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
    p + geom_text(
        data = tags, aes_string(x = "x", y = "y", label = "label"), ...,
        hjust = hjust, vjust = vjust, fontface = fontface, family = family,
        inherit.aes = FALSE)
}

pch_map <- c(
    "Spiro_C" = 0, "Spiro_D" = 5,
    "Tetra_1" = 2, "Tetra_2" = 6,
    "Para_4" = 1,
    "Loxo_1" = 3, "Loxo_2" = 4,
    "VIII_2383" = 7, "VIII_2434" = 9,
    "V_2375" = 10, "V_2524" = 13
)

pch_merged <- scale_shape_manual(values = pch_map)
pch_cilia <- scale_shape_manual(values = pch_map[1:7])
pch_cyano <- scale_shape_manual(values = pch_map[8:11])

# model_system <- "cilia"
# model_system <- "cyano"

# some plots separate
for (model_system in c("cilia", "cyano")) {

    # path to save figures in ----
    # outpath <- paste0("figures/plot_data/", model_system, "_")
    outpath <- paste0("figures/", model_system, "_")

    # load the raw data ----
    data <- get(paste("data_", model_system, sep = ""))
    long_data <- get(paste("long_data_", model_system, sep = ""))

    print(paste("Loaded", model_system, "data"))

    # plot pcgr vs. pop. ----
    p <- ggplot(long_data) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        aes(x = day, y = value, col = treat) +
        geom_point(size = 1) +
        facet_grid2(name ~ strain, independent = "y",
            scales = "free", switch = "y", labeller = label_parsed) +
        labs(x = "Time (days)",
            y = NULL,
            col = "Treatment") +
        theme(
            legend.position = "bottom",
            strip.placement.y = "outside",
            strip.background.y = element_rect(fill = NA, color = NA))

    p <- tag_facet(p)

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "pop-traits-time.png"),
            width = 5, height = 6.5, device = "png", dpi = 450)
    } else {
        ggsave(paste0(outpath, "pop-traits-time.png"),
            width = 8, height = 6.5, device = "png", dpi = 450)
    }

    print(paste("Saved", model_system, "basic plots"))

    # plot pcgr vs. pop. ----
    p <- ggplot(data) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        aes(x = density, y = pcgr, col = treat) +
        geom_hline(yintercept = 0) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5, formula = "y ~ x") +
        facet_wrap(vars(strain), ncol = 2, scales = "free") +
        labs(x = "Density",
            y = "PCGR",
            col = "Treatment")

    p <- tag_facet(p)

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "pcgr.pdf"),
            width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "pcgr.pdf"),
            width = 5, height = 8, device = "pdf")
    }

    print(paste("Saved", model_system, "PCGR plots"))

    # plot dT vs. trait ----
    p <- ggplot(data) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        aes(x = trait, y = dT, col = treat) +
        geom_hline(yintercept = 0) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5, formula = "y ~ x") +
        facet_wrap(vars(strain), ncol = 2, scales = "free") +
        labs(x = "Trait value",
            y = "Trait change",
            col = "Treatment")

    p <- tag_facet(p)

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "dT.pdf"),
            width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "dT.pdf"),
            width = 5, height = 8, device = "pdf")
    }

    print(paste("Saved", model_system, "dT plots"))

    # Check correlations between traits and abundance ----
    p <- ggplot(data) +
        scale_shape_manual(values = 0:10) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        aes(x = log10(density), y = trait, col = treat) +
        geom_point(pch = 1) +
        facet_wrap(vars(strain), ncol = 2) +
        labs(col = "treatment") +
        labs(x = expression(paste("log"[10], " density")),
            y = "Trait value",
            col = "Treatment")

    p <- tag_facet(p)

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "corr.pdf"),
            width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "corr.pdf"),
            width = 5, height = 8, device = "pdf")
    }
    print(paste("Saved", model_system, "corr plots"))
}

# stats of the model outputs ----
stats_result <- modelling(
        data = data_merged,
        var_to_nest_by = c("strain", "treat"),
        formulas = c(
            "dT ~ trait",
            "dT ~ trait + density",
            "pcgr ~ density",
            "pcgr ~ trait + density")) %>%
    rowwise %>%
    mutate(
        model_summary = list(tibble(
            coefficient = names(coef(model)),
            estimate = unname(coef(model)),
            "Std. Error" = summary(model)$coefficients[, "Std. Error"],
            "Conf. Int." = confint(model),
            "Pr(>|t|)" = summary(model)$coefficients[, "Pr(>|t|)"]
            )),
        obs_pred = list(tibble(
            obs = model$model[, 1],
            pred = model$fitted.values,
            error = abs(model$residuals)
            )),
        R2 = abs(summary(model)$r.squared),
        aic = MuMIn::AICc(model),
        response = ifelse(grepl("dT ~ ", form), "trait change", "growth"),
        response = factor(response, levels = c("growth", "trait change")),
        predictor = ifelse(grepl("\\+", form), "augmented", "single"),
        predictor = factor(predictor, levels = c("single", "augmented"))
        ) %>%
    left_join(data_merged %>%
        dplyr::select(system, strain, treat, pca_varexp) %>%
        distinct, relationship = "many-to-many")

print("Created statistics")

# create bars to better separate strains
vline_pos <- (0.5):(length(unique(stats_result$strain)) - 0.5)

# Plot the PCA variance explained
p <- ggplot(stats_result) +
    aes(x = strain, y = pca_varexp, fill = treat) +
    facet_grid2(. ~ system, scale = "free_x", space = "free_x") +
    scale_x_discrete() +
    geom_vline(xintercept = vline_pos, color = "grey92") +
    stat_summary(
        inherit.aes = FALSE,
        fun.y = mean,
        aes(x = 2, y = pca_varexp, yintercept = after_stat(y)),
        geom = "hline", lty = 2) +
    geom_bar(stat = "identity", position = "dodge", width = 0.75) +
    labs(y = "Proportion of variance explained by\nfirst principal component",
        x = "Strain",
        fill = "Treatment") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = cbPalette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
        labels = scales::percent) +
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

ggsave("figures/PC_var_explained.pdf", p,
        width = 6, height = 2.5, device = "pdf")

print(paste("Saved PCA variance explained plots"))

# plot R squared

# plot
p <- ggplot(stats_result) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = cbPalette) +
    aes(x = strain, y = R2, fill = treat) +
    scale_x_discrete() +
    geom_vline(xintercept = vline_pos, color = "grey92") +
    stat_summary(
        inherit.aes = FALSE,
        fun.y = mean,
        aes(x = 2, y = R2, yintercept = after_stat(y)),
        geom = "hline", lty = 2) +
    geom_bar(stat = "identity", position = "dodge", width = 0.75) +
    geom_hline(yintercept = 0, col = "black") +
    geom_hline(yintercept = 1, col = "grey46") +
    facet_grid2(response + predictor ~ system, scales = "free",
        labeller = labeller(.default = capitalise), space = "free_x",
        strip = strip_nested()) +
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    labs(y = expression("R "^2),
        x = "Strain",
        fill = "Treatment")

p <- tag_facet(p)

ggsave("figures/R2.pdf", p, width = 5, height = 6, device = "pdf")

print("Saved R-squared plots")

# join models to data and make predictions ----
data_preds <- stats_result %>%
    unnest(obs_pred) %>%
    dplyr::select(-data, -model, -model_summary)

# form "pcgr ~ density"
p <- ggplot(data_preds) +
    theme_bw() +
    pch_merged +
    scale_colour_manual(values = cbPalette) +
    aes(x = obs, y = pred, col = treat, pch = strain) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_smooth(method = "lm", se = FALSE, lwd = 0.5, formula = "y ~ x") +
    facet_nested(response ~ system + predictor, scales = "free", switch = "y",
        labeller = labeller(.default = capitalise), independent = TRUE) +
    geom_abline(slope = 1, intercept = 0) +
    labs(x = "Observed value",
        y = "Predicted value",
        pch = "Strain",
        col = "Treatment") +
    theme(strip.placement.y = "outside",
	    strip.background.y = element_blank(),
        strip.text.y = element_text(colour = "black")) +
    guides(pch = guide_legend(override.aes = list(size = 2, alpha = 1)),
        col = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

p <- tag_facet(p)

pdf("figures/growth_dtrait.pdf", width = 4.5, height = 9)
grid::pushViewport(grid::viewport(name = "rotate", angle = 270, width = 2, height = 0.5))
print(p, vp = "rotate")
dev.off()

print(paste("Saved", model_system, "obs. vs pred. plots"))
    
# Now error plot
data_synth <- data_preds %>%
    group_by(system, strain, treat, response, predictor) %>%
    summarise(error = sum(error), aic = first(aic), pca_varexp = first(pca_varexp)) %>%
    pivot_wider(names_from = predictor, values_from = error:aic) %>%
    rowwise %>%
    mutate(delta_error = (error_augmented - error_single) / error_single,
        AIC.weights = list(akaike.weights(c(aic_augmented, aic_single))$weights),
        p_better = AIC.weights[1])

# view the delta error of using the augmented model vs single model
p <- ggplot(data_synth) +
    scale_shape_manual(values = 0:10) +
    theme_bw() +
    scale_fill_manual(values = cbPalette) +
    aes(x = strain, y = delta_error, fill = treat) +
    scale_x_discrete() +
    geom_vline(xintercept = vline_pos, color = "grey92") +
    stat_summary(
        inherit.aes = FALSE,
        fun.y = mean,
        aes(x = 2, y = delta_error, yintercept = after_stat(y)),
        geom = "hline", lty = 2) +
    geom_bar(stat = "identity", position = "dodge", width = 0.75) +
    geom_hline(yintercept = 0) +
    facet_grid2(response ~ system, labeller = labeller(.default = capitalise),
        scales = "free_x", space = "free_x") +
    labs(x = "Strain",
        fill = "Treatment",
        y = expression(paste(
            "(Error"[augmented], " - Error"[single], ") /  Error"[single]))) +
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

p <- tag_facet(p)

ggsave("figures/delta_error.pdf", p, width = 4.5, height = 4, device = "pdf")

print(paste("Saved", model_system, "delta error plots"))
    
# view the probability that augmented model is best
p <- ggplot(data_synth) +
    scale_shape_manual(values = 0:10) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = cbPalette) +
    aes(x = strain, y = p_better, fill = treat) +
    scale_x_discrete() +
    geom_vline(xintercept = vline_pos, color = "grey92") +
    stat_summary(
        inherit.aes = FALSE,
        fun.y = mean,
        aes(x = 2, y = p_better, yintercept = after_stat(y)),
        geom = "hline", lty = 2) +
    geom_bar(stat = "identity", position = "dodge", width = 0.75) +
    # stat_summary(geom = "hline", fun. = "mean") +
    geom_hline(yintercept = 0, col = "black") +
    geom_hline(yintercept = 1, col = "grey46") +
    facet_grid2(response ~ system, labeller = labeller(.default = capitalise),
        scales = "free_x", space = "free_x") +
    labs(x = "Strain",
        fill = "Treatment",
        y = "Probability that augmented model is best") +
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )

p <- tag_facet(p)

ggsave("figures/AIC.pdf", width = 6, height = 4, device = "pdf")

print(paste("Saved", model_system, "AIC plots"))

p <- ggplot(data_synth) +
    aes(pch = strain, x = p_better, y = -delta_error, col = treat) +
    theme_bw() +
    # scale_x_log10() +
    scale_color_manual(values = cbPalette) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = cbPalette) +
    geom_hline(yintercept = 0) +
    geom_point(size = 2, stroke = 1.5) +
    facet_grid2(response ~ system, labeller = labeller(.default = capitalise)) +
    labs(pch = "Strain",
        y = "Predictive accuracy difference",
        color = "Treatment",
        x = "Probability that augmented model is best") +
    # theme(legend.position = "bottom", legend.box = "vertical") +
    pch_merged +
    guides(pch = guide_legend(override.aes = list(size = 2, alpha = 1, stroke = 0.5)),
        col = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

p <- tag_facet(p)

ggsave("figures/AIC_deltaerror.pdf", p, width = 6, height = 4.75, device = "pdf")

print(paste("Saved", model_system, "AIC vs delta error plots"))

# results about regression coefficients ----
stats_coef <- stats_result %>%
    dplyr::select(system, strain, treat, form, response,
        predictor, model_summary) %>%
    unnest(model_summary) %>%
    ungroup() %>%
    rowwise %>%
    mutate(estimate_sig = ifelse(`Pr(>|t|)` < 0.05, estimate, NA),
        CI_low = `Conf. Int.`[1],
        CI_high = `Conf. Int.`[2])

# coefficient values and significance
plotCoef <- function(dat) {
    p <- ggplot(dat) +
        aes(y = strain, x = estimate, col = treat) +
        geom_abline(intercept = 0, slope = 1e16) +
        geom_point(shape = 1, size = 3,
            position = position_dodge(width = 0.5)) +
        scale_x_continuous(n.breaks = 3) +
        # scale_x_continuous(trans = scales::pseudo_log_trans()) +
        geom_errorbarh(
            aes(xmin = CI_low, xmax = CI_high),
            position = position_dodge(width = 0.5), height = .2) +
        geom_point(aes(y = strain, x = estimate_sig),
            position = position_dodge(width = 0.5),
            pch = "*", cex = 5, show.legend = FALSE) +
        facet_grid2(response + predictor ~ coefficient, scales = "free_x",
            independent = "x", render_empty = FALSE, strip = strip_nested(),
            labeller = labeller(.default = capitalise)) +
        # theme(axis.text.x = element_text(angle = 90)) +
        labs(x = "Estimate", y = "Strain", col = "Treatment") +
        theme_bw() +
        scale_colour_manual(values = cbPalette) 
    # p <- tag_facet(p)
}

# ciliate plot
p <- plotCoef(subset(stats_coef, system == "Ciliates"))
ggsave("figures/cilia_td_general.pdf", p, width = 6.5, height = 7.5)

p <- plotCoef(subset(stats_coef, system == "Cyanobacteria"))
ggsave("figures/cyano_td_general.pdf", width = 6.5, height = 7.5)

print(paste("Saved", model_system, "model coefficient plots"))

# varexp vs predictive acc
p <- data_synth %>%
    ggplot() +
        aes(x = pca_varexp, y = -delta_error) +
        theme_bw() +
        # scale_x_log10() +
        scale_color_manual(values = cbPalette) +
        scale_y_continuous(labels = scales::percent) +
        geom_hline(yintercept = 0) +
        geom_point(aes(col = treat, pch = strain), size = 2) +
        geom_smooth(aes(col = treat), method = "lm", se = FALSE) +
        # geom_smooth(method = "lm") +
        facet_grid2(response ~ .,
            labeller = labeller(.rows = label_both, .default = capitalise)) +
        # facet_grid2(. ~ system, labeller = labeller(.default = capitalise)) +
        pch_merged +
        guides(pch = guide_legend(override.aes = list(size = 2, alpha = 1)),
            col = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
        labs(
            x = "Proportion of variance explained by\nfirst principal component",
            y = "Predictive accuracy difference",
            col = "Treatment", pch = "Strain")

p <- tag_facet(p)

ggsave("figures/varexp_deltaerror.pdf", width = 4, height = 5)
