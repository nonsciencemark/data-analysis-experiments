# load in the data and get plotting functions
source("Tools and data.R")

library("qpcR") # for AIC comparison
library("lmtest") # for Granger causality
library("ggh4x") # plotting stuff
library("qpdf") # merging pdf at the end

# model_system <- "cilia"
# model_system <- "cyano"

for (model_system in c("cilia", "cyano")) {

    # path to save figures in ----
    outpath <- paste0("figures/", model_system, "/", model_system, "_")

    # load the raw data ----
    data <- get(paste("data_", model_system, sep = ""))

    print(paste("Loaded", model_system, "data"))

    # plot pcgr vs. pop. ----
    ggplot(data) +
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

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "pcgr.pdf"),
                     width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "pcgr.pdf"),
                     width = 5, height = 8, device = "pdf")
    }

    print(paste("Saved", model_system, "PCGR plots"))

    # plot dT vs. trait ----
    ggplot(data) +
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

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "dT.pdf"),
            width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "dT.pdf"),
            width = 5, height = 8, device = "pdf")
    }

    print(paste("Saved", model_system, "dT plots"))

    # Check correlations between traits and abundance ----
    ggplot(data) +
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

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "corr.pdf"),
            width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "corr.pdf"),
            width = 5, height = 8, device = "pdf")
    }
    
    print(paste("Saved", model_system, "corr plots"))

    # stats of the model outputs ----
    stats_result <- modelling(
            data = data,
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
                "Pr(>|t|)" = summary(model)$coefficients[, "Pr(>|t|)"]
                )),
            # intrinsic_growth = model_summary["(Intercept)", "Pr(>|t|)"],
            obs_pred = list(tibble(
                obs = model$model[, 1],
                pred = model$fitted.values,
                error = abs(model$residuals)
                )),
            R2 = abs(summary(model)$r.squared),
            aic = AIC(model),
            response = ifelse(grepl("dT ~ ", form), "trait change", "growth"),
            response = factor(response, levels = c("growth", "trait change")),
            predictor = ifelse(grepl("\\+", form), "both", "single"),
            predictor = factor(predictor, levels = c("single", "both"))
            ) %>%
        left_join(data %>% 
            dplyr::select(strain, treat, pca_varexp) %>%
            distinct, relationship = "many-to-many")

    print(paste("Created", model_system, "statistics"))

    # create bars to better separate strains
    vline_pos <- (1.5):(length(unique(stats_result$strain)) - 0.5)

    # Plot the PCA variance explained
    ggplot(stats_result) +
        aes(x = strain, y = pca_varexp, fill = treat) +
        scale_x_discrete() +
        geom_vline(xintercept = vline_pos, color = "grey92") +
        stat_summary(
            inherit.aes = FALSE,
            fun.y = mean,
            aes(x = 2, y = pca_varexp, yintercept = after_stat(y)),
            geom = "hline", lty = 2) +
        geom_bar(stat = "identity", position = "dodge", width = 0.75) +
        labs(y = "Variance explained by\nfirst principal component",
            x = "Strain",
            fill = "Treatment") +
        theme_bw() +
        scale_fill_manual(values = cbPalette) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
            labels = scales::percent) +
        theme(panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank())

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "PC_var_explained.pdf"),
            width = 5.5, height = 2.5, device = "pdf")
    } else {
        ggsave(paste0(outpath, "PC_var_explained.pdf"),
            width = 5.5, height = 2.5, device = "pdf")
    }

    print(paste("Saved", model_system, "PCA variance explained plots"))

    # plot R squared

    # plot
    ggplot(stats_result) +
        scale_shape_manual(values = 0:10) +
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
        geom_hline(yintercept = c(0, 1)) +
        facet_grid2(response + predictor ~ ., scales = "free",
            strip = strip_nested()) +
        theme(panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank()) +
        labs(y = expression("R "^2),
            x = "Strain",
            fill = "Treatment")

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "R2.pdf"),
            width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "R2.pdf"),
            width = 6, height = 4, device = "pdf")
    }

    print(paste("Saved", model_system, "R-squared plots"))

    # join models to data and make predictions ----
    data_preds <- stats_result %>%
        unnest(obs_pred) %>%
        dplyr::select(-data, -model, -model_summary)

    # form "pcgr ~ density"
    ggplot(data_preds) +
        theme_bw() +
        scale_shape_manual(values = 0:10) +
        scale_colour_manual(values = cbPalette) +
        aes(x = obs, y = pred, col = treat, pch = strain) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5, formula = "y ~ x") +
        facet_grid2(predictor ~ response, scales = "free", independent = TRUE) +
        geom_abline(slope = 1, intercept = 0) +
        labs(x = "Observed value",
            y = "Predicted value",
            pch = "Strain",
            col = "Treatment")

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "growth_dtrait.pdf"),
            width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "growth_dtrait.pdf"),
            width = 5, height = 4, device = "pdf")
    }

    print(paste("Saved", model_system, "obs. vs pred. plots"))
    
    # Now error plot
    data_synth <- data_preds %>%
        group_by(strain, treat, response, predictor) %>%
        summarise(error = sum(error), aic = first(aic)) %>%
        pivot_wider(names_from = predictor, values_from = error:aic) %>%
        rowwise %>%
        mutate(delta_error = (error_both - error_single) / error_single,
            AIC.weights = list(akaike.weights(c(aic_both, aic_single))$weights),
            p_better = AIC.weights[1])

    # view the delta error of using the full model vs single model
    ggplot(data_synth) +
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
        facet_grid(response ~ .) +
        labs(x = "Strain",
            fill = "treatment",
            y = expression(paste(
                "(Error"[full], " - Error"[single], ") /  Error"[single]))) +
        theme(panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank())

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "delta_error.pdf"),
            width = 4.5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "delta_error.pdf"), # ciliate plot bigger
            width = 6, height = 4, device = "pdf")
    }

    print(paste("Saved", model_system, "delta error plots"))
    
    # view the probability that full model is best
    ggplot(data_synth) +
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
        geom_hline(yintercept = c(0, 1)) +
        facet_grid(response ~ .) +
        labs(x = "Strain",
            fill = "Treatment",
            y = "Probability that full model is best") +
        theme(panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank()
        )

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "AIC.pdf"),
            width = 4.5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "AIC.pdf"), # ciliate plot bigger
            width = 6, height = 4, device = "pdf")
    }

    print(paste("Saved", model_system, "AIC plots"))

    ggplot(data_synth) +
        aes(pch = strain, x = p_better, y = -delta_error, col = treat) +
        scale_shape_manual(values = 0:10) +
        theme_bw() +
        scale_x_log10() +
        scale_color_manual(values = cbPalette) +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values = cbPalette) +
        geom_hline(yintercept = 0) +
        geom_point(size = 2, stroke = 1.5) +
        facet_wrap(response ~ .) +
        labs(pch = "Strain",
            y = "Predictive accuracy difference",
            color = "Treatment",
            x = "Probability that full model is best") +
        theme(legend.position = "bottom", legend.box = "vertical")

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "AIC_deltaerror.pdf"),
            width = 6, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "AIC_deltaerror.pdf"),
            width = 6, height = 4.5, device = "pdf")
    }

    print(paste("Saved", model_system, "AIC vs delta error plots"))

    # results about regression coefficients ----
    stats_coef <- stats_result %>%
        dplyr::select(strain, treat, form, response,
            predictor, model_summary) %>%
        unnest(model_summary) %>%
        ungroup() %>%
        mutate(estimate_sig = ifelse(`Pr(>|t|)` < 0.05, estimate, NA))

    # coefficient values and significance
    ggplot(stats_coef) +
        aes(y = strain, x = estimate, col = treat) +
        geom_abline(intercept = 0, slope = 1e16) +
        geom_point(shape = 1, size = 3,
            position = position_dodge(width = 0.5)) +
        geom_errorbarh(
            aes(xmin = estimate - `Std. Error`, xmax = estimate + `Std. Error`),
            position = position_dodge(width = 0.5), height = .2) +
        geom_point(aes(y = strain, x = estimate_sig),
            position = position_dodge(width = 0.5),
            pch = "*", cex = 5, show.legend = FALSE) +
        facet_grid2(response + predictor ~ coefficient, scales = "free_x",
            independent = "x", render_empty = FALSE, strip = strip_nested()) +
        theme(axis.text.x = element_text(angle = 90)) +
        labs(x = "Estimate", y = "Strain", col = "Treatment") +
        theme_bw() +
        scale_colour_manual(values = cbPalette)

    if (model_system == "cyano") {
        ggsave(paste(outpath, "td_general.pdf", sep = ""),
            width = 6.5, height = 7.5)
    } else {
        ggsave(paste(outpath, "td_general.pdf", sep = ""),
            width = 6.5, height = 7.5)
    }

    print(paste("Saved", model_system, "model coefficient plots"))

    # merge pdfs to see them all
    all_pdfs <- list.files(paste0("figures/", model_system), 
        pattern = "\\.pdf$", full.names = TRUE)
    all_pdfs <- all_pdfs[!grepl("merged.pdf", all_pdfs)]
    pdf_combine(input = all_pdfs,
        output = paste0("figures/", model_system, "/merged.pdf"))

    print(paste("Saved", model_system, "model merged plots file"))

}
