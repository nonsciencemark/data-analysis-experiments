# load in the data and get plotting functions
source("Tools and data.R")

library("qpcR") # for AIC comparison
library("lmtest") # for Granger causality
library("ggh4x") # plotting stuff
library("egg")

# model_system <- "cilia"
# model_system <- "cyano"

for (model_system in c("cilia", "cyano")) {

    # path to save figures in
    outpath <- paste0("plots/", model_system, "/", model_system, "_")

    # GET THE DATA
    data <- get(paste("data_", model_system, sep = ""))

    # plot pcgr vs. pop. ----
    ggplot(data) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        aes(x = density, y = pcgr, col = treat) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
        facet_wrap(vars(strain), ncol = 2, scales = "free")

    # save
    if (model_system == "cyano") {
        ggsave(paste0(outpath, "pcgr.pdf"), 
                     width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "pcgr.pdf"),
                     width = 5, height = 8, device = "pdf")
    }

    # plot dT vs. trait ----
    ggplot(data) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        aes(x = trait, y = dT, col = treat) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
        facet_wrap(vars(strain), ncol = 2, scales = "free")

    # save
    if (model_system == "cyano") {
        ggsave(paste0(outpath, "dT.pdf"),
            width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "dT.pdf"),
            width = 5, height = 8, device = "pdf")
    }

    # join models to data and make predictions ----
    data_preds <- data %>%
        ungroup() %>%
        nest_by(strain) %>%
        left_join(stats_result, by = c("strain"), multiple = "all") %>%
        mutate(response = case_when(
            length(grep("dT", form)) > 0 ~ "trait change",
            length(grep("pcgr", form)) > 0 ~ "growth")) %>%
        mutate(predictor = case_when(
            length(grep("trait + density", form, fixed = TRUE)) > 0 ~ "both",
            length(grep("density", form)) > 0 ~ "one",
            length(grep("trait", form)) > 0 ~ "one")) %>%
        mutate(data_and_pred = list(cbind(
            data,
            prediction = predict.lm(object = model, newdata = data)))) %>%
        dplyr::select(c("strain", "response", "form",
            "predictor", "data_and_pred")) %>%
        unnest(c("data_and_pred")) %>%
        mutate(error = case_when(
            response == "trait change" ~ abs(prediction - dT),
            response == "growth" ~ abs(prediction - pcgr)))

    # Plot model fits
    data_preds <- data_preds %>%
        dplyr::select(-c("data", "model")) %>%
        unnest(predictions)

    # form "pcgr ~ density"
    ggplot(data_preds %>% filter(response == "growth")) +
        theme_bw() +
        scale_shape_manual(values = 0:10) +
        scale_colour_manual(values = cbPalette) +
        aes(x = pcgr, y = pred, col = treat, pch = strain) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
        facet_grid(form ~ .) +
        geom_abline(slope = 1, intercept = 0) +
        labs(col = "treatment")

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "growth.pdf"),
            width = 4, height = 4.5, device = "pdf")
    } else {
        ggsave(paste0(outpath, "growth.pdf"),
            width = 4, height = 4.5, device = "pdf")
    }

    ggplot(data_preds %>% filter(response == "trait change")) +
        theme_bw() +
        scale_shape_manual(values = 0:10) +
        scale_colour_manual(values = cbPalette) +
        aes(x = dT, y = pred, col = treat, pch = strain) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
        facet_grid(form ~ .) + #
        geom_abline(slope = 1, intercept = 0)    +
        labs(col = "treatment")

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "trait.pdf"),
            width = 4, height = 4.5, device = "pdf")
    } else {
        ggsave(paste0(outpath, "trait.pdf"),
            width = 4, height = 4.5, device = "pdf")
    }

    #Now plot, basic plot first: growth
    ggplot(data_preds %>% filter(response == "growth")) +
        aes(x = pcgr, y = prediction, col = treat, pch = strain) + #
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        facet_grid(form ~ .) + #
        geom_abline(intercept = 0, slope = 1)

    #Now plot trait change
    ggplot(data_preds %>% filter(response == "trait change")) +
        aes(x = dT, y = prediction, col = treat, pch = strain) + #
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        facet_grid(form ~ .) + #
        geom_abline(intercept = 0, slope = 1)

    #Now error plot
    data_preds_synth <- data_preds %>%
        group_by(strain, response, predictor, treat) %>%
        summarise(error = sum(error)) %>%
        pivot_wider(names_from = predictor, values_from = error) %>%
        mutate(delta_error = (both - one) / one)

    ggplot(data_preds_synth) +
        scale_shape_manual(values = 0:10) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        aes(x = response, y = delta_error, col = treat) +
        geom_jitter(width = 0.25) +
        geom_hline(yintercept = 0) +
        facet_wrap(vars(strain), ncol = 2) +
        labs(y = expression(paste("Error"[full],"-Error"[single])),
            col = "treatment")

# stats of the model outputs
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
            mod_coefs = list(summary(model)$coefficients),
            intrinsic_growth = mod_coefs["(Intercept)", "Pr(>|t|)"],
            model_summary = list(as_tibble(
                rownames_to_column(
                    as.data.frame(summary(model)$coefficients),
                        var = "predictor"))),
            R2 = abs(summary(model)$r.squared),
            aic = AIC(model),
            response = ifelse(grepl("dT ~ ", form), 
                "trait change", "growth")
            # type_of_pred = ifelse(
            # length(grep("density", predictor) > 0) |
            #     length(grep("trait", predictor) > 0),
            # "slope", "intercept")
            )

    ggplot(stats_result) +
        aes(x = intrinsic_growth, y = aic) +
        geom_point()

    stats_result <- modelling(data = data,
        var_to_nest_by = c("strain", "treat"),
        formulas = c("dT ~ trait")) %>%
        mutate(model_summary =
            list(as_tibble(rownames_to_column(
                as.data.frame(summary(model)$coefficients),
                var = "predictor")))) %>%
        unnest(model_summary)    %>%
        rowwise %>%
        mutate(type_of_pred = ifelse(
            length(grep("density", predictor) > 0) |
                length(grep("trait", predictor) > 0),
            "slope", "intercept")) %>%
        ungroup() %>%
        mutate(Estimate_sig = ifelse(`Pr(>|t|)` < 0.05, Estimate, NA))

    # Now join these models to all the data and make the predictions
    data_preds <- data %>%
        ungroup() %>%
        nest_by(species, strain, treat) %>%
        left_join(stats_result, by = c("strain", "treat"), multiple = "all") %>%
        mutate(predictions = list(cbind(data,
            pred = predict.lm(model, newdata = data)))) %>%
        rowwise %>%
        mutate(AIC = AIC(model)) %>%
        ungroup() %>%
        mutate(response = ifelse(grepl("dT ~ ", form),
            "trait change", "growth"))

    data_preds_synth <- data_preds %>%
        mutate(predictors = ifelse(grepl("\\+", form), "both", "focal")) %>%
        dplyr::select(c("species", "strain", "treat",
            "response", "AIC", "predictors")) %>%
        pivot_wider(names_from = "predictors", values_from = "AIC") %>%
        rowwise %>%
        mutate(
            AIC.weights = list(akaike.weights(c(both, focal))$weights),
            p.better = AIC.weights[1])

    plot(
        data_preds_synth %>%
            dplyr::select(-p.better) %>%
            unnest_longer(p.better)
    )

    ggplot(data_preds_synth) +
        scale_shape_manual(values = 0:10) +
        theme_bw() +
        scale_fill_manual(values = cbPalette) +
        aes(x = response, y = p.better, fill = treat) +
        geom_col(position = position_dodge(width = 0.5), width = 0.25) +
        geom_hline(yintercept = c(0, 1)) +
        facet_wrap(vars(strain), ncol = 2, scales = "free") +
        labs(col = "treatment",
                 y = expression("Prob. full model outperforms single model"))

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "AIC.pdf"),
            width = 4.5, height = 3.5, device = "pdf")
    } else {
        ggsave(paste0(outpath, "AIC.pdf"), # ciliate plot bigger
            width = 4.5, height = 7, device = "pdf")
    }

    ## Check correlations between traits and abundance ----

    ggplot(data) +
        scale_shape_manual(values = 0:10) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        aes(x = log10(density), y = trait, col = treat) +
        geom_point(pch = 1) +
        facet_wrap(vars(strain), ncol = 2) +
        labs(col = "treatment")

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "corr.pdf"),
            width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "corr.pdf"),
            width = 5, height = 8, device = "pdf")
    }

    plot <- ggplot(stats_result) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        aes(x = strain, y = Estimate, col = treat) +
        geom_point(shape = 1, size = 3,
            position = position_dodge(width = 0.5)) +
        geom_errorbar(
            aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`),
            position = position_dodge(width = 0.5), width = .2) +
        geom_point(aes(x = strain, y = Estimate_sig),
                             position = position_dodge(width = 0.5),
                             pch = "*", cex = 5, show.legend = F) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        coord_flip() +
        facet_wrap(vars(type_of_pred), scales = "free") +
        theme(axis.text.x = element_text(angle = 90)) +
        labs(col = "treatment")

    if (model_system == "cyano") {
        ggsave(paste(outpath, "td_general.pdf", sep = ""), plot = plot,
            width = 7, height = 7)
    } else {
        ggsave(paste(outpath, "td_general.pdf", sep = ""), plot = plot,
            width = 7, height = 7)
    }

    ggplot(stats_result) +
        scale_shape_manual(values = 0:10) +
        theme_bw() +
        scale_fill_manual(values = cbPalette) +
        aes(x = form, y = R2_max, fill = treat) +
        geom_col(position = position_dodge(width = 0.5), width = 0.25) +
        geom_hline(yintercept = c(0, 1)) +
        facet_wrap(vars(strain), ncol = 2, scales = "free") +
        labs(col = "treatment",
            y = expression("Prob. full model outperforms single model"))

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "R2.pdf"),
            width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "R24.pdf"),
            width = 5, height = 8, device = "pdf")
    }

    ggplot(stats_result) +
        scale_shape_manual(values = 0:10) +
        theme_bw() +
        scale_fill_manual(values = cbPalette) +
        aes(x = form, y = R2_max, fill = treat) +
        geom_col(position = position_dodge(width = 0.5), width = 0.25) +
        geom_hline(yintercept = c(0, 1)) +
        facet_wrap(vars(strain), ncol = 2, scales = "free") +
        labs(col = "treatment",
            y = expression("Prob. full model outperforms single model"))

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "R2.pdf"),
            width = 5, height = 4, device = "pdf")
    } else {
        ggsave(paste0(outpath, "R2.pdf"),
            width = 5, height = 8, device = "pdf")
    }
}
