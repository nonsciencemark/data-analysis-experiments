source("Tools and data.R")

library("qpcR") # for AIC comparison
library("lmtest") # for Granger causality
library("ggh4x") # plotting stuff
library("egg")

# model_system <- "cilia"
# model_system <- "cyano"
model_sytems <- c("cilia", "cyano")

for (model_system in model_systems) {

    outpath <- paste0("plots/", model_system, "/", model_system, "_")

    # GET THE DATA
    data <- get(paste("data_", model_system, sep = ""))

    # BASIC OF DNDT-D AND DTDT-T

    # plot pcgr vs. pop. ----
    ggplot(data) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        aes(x = density, y = pcgr, col = treat) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
        facet_wrap(vars(strain), ncol = 2, scales = "free")

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

stats_result <- modelling(data = data,
    var_to_nest_by = c("strain", "treat"),
    formulas = c("dT ~ trait",
        "dT ~ trait + density",
        "pcgr ~ density",
        "pcgr ~ trait + density")) %>%
        dplyr::select(-data)

    # Now join these models to all the data and make the predictions
    data_preds <- data %>%
        ungroup() %>%
        nest_by(species, strain, treat) %>%
        left_join(stats_result, by = c("strain", "treat"), multiple = "all") %>%
        mutate(predictions = list(cbind(data, 
            pred = predict.lm(model, newdata = data)))) %>%
        rowwise() %>%
        mutate(AIC = AIC(model),
                     ) %>%
        ungroup() %>%
        mutate(response = ifelse(grepl('dT ~ ', form), 
            "trait change", "growth")) %>%
        ungroup()

    data_preds_synth <- data_preds %>%
        mutate(predictors = ifelse(grepl("\\+", form), "both", "focal")) %>%
        dplyr::select(c("species", "strain", "treat", 
            "response", "AIC", "predictors")) %>%
        pivot_wider(names_from = "predictors", values_from = "AIC") %>%
        rowwise %>%
        mutate(
            AIC.weights = list( akaike.weights(c(both, focal))$weights),
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
        aes(x     = response, y = p.better, fill = treat) +
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

    # Plot model fits
    data_preds <- data_preds %>%
        dplyr::select(-c("data", "model")) %>%
        unnest(predictions)

    # form == "pcgr ~ density"
    ggplot(data_preds %>% filter(response == "growth")) +
        theme_bw() +
        scale_shape_manual(values = 0:10) +
        scale_colour_manual(values = cbPalette) +
        aes(x = pcgr, y = pred, col = treat, pch = strain) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5)+
        facet_grid(form~.) + #, scales = "free"
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
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5)+
        facet_grid(form~.) + #
        geom_abline(slope = 1, intercept = 0)    +
        labs(col = "treatment")

    if (model_system == "cyano") {
        ggsave(paste0(outpath, "trait.pdf"),
                     width = 4, height = 4.5, device = "pdf")
    } else {
        ggsave(paste0(outpath, "trait.pdf"),
                     width = 4, height = 4.5, device = "pdf")
    }

    #Now join these models to all the data and make the predictions
    data_preds <- data %>%
        ungroup() %>%
        nest_by(strain) %>%
        left_join(stats_result, by = c("strain"), multiple = "all") %>%
        mutate(response = case_when(
            length(grep("dT", form)) > 0 ~ "trait change",
            length(grep("pcgr", form)) > 0 ~ "growth")) %>%
        mutate(predictor = case_when(
            length(grep("trait + density", form, fixed = T)) > 0 ~ "both",
            length(grep("density", form)) > 0 ~ "one",
            length(grep("trait", form)) > 0 ~ "one")) %>%
        #dplyr::select(-form) %>%
        mutate(data_and_pred = list(cbind(
            data,
            prediction = predict.lm(object = model, newdata = data)))) %>%
        dplyr::select(c("strain", "response", "form", 
            "predictor", "data_and_pred")) %>%
        unnest(c("data_and_pred")) %>%
        mutate(error = case_when(
            response == "trait change" ~ abs(prediction - dT),
            response == "growth" ~ abs(prediction - pcgr)))

    #Now plot, basic plot first: growth
    ggplot(data_preds %>% filter(response == "growth")) +
        aes(x = pcgr, y = prediction, col = treat, pch = strain) + #
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        facet_grid(form~.) + #
        geom_abline(intercept = 0, slope = 1)

    #Now plot trait change
    ggplot(data_preds %>% filter(response == "trait change")) +
        aes(x = dT, y = prediction, col = treat, pch = strain) + #
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        facet_grid(form~.) + #
        geom_abline(intercept = 0, slope = 1)

    #Now error plot
    data_preds_synth <- data_preds %>%
        group_by(strain, response, predictor, treat) %>%
        summarise(error = sum(error)) %>%
        pivot_wider(names_from = predictor, values_from = error) %>%
        mutate(delta_error = (both - one)/one)

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

    # Leftovers ----

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


    stats_result <- modelling(data = data,
        var_to_nest_by = c("strain", "treat"),
        formulas = c("dT ~ trait")) %>%
        mutate(model_summary =
            list(as_tibble(rownames_to_column(
                as.data.frame(summary(model)$coefficients),
                var = "predictor")))) %>%
        unnest(model_summary)    %>%
        rowwise() %>%
        mutate(type_of_pred = ifelse(
            length(grep("density", predictor) > 0) |
                length(grep("trait", predictor) > 0),
            "slope", "intercept")) %>%
        ungroup() %>%
        mutate(Estimate_sig = ifelse(`Pr(>|t|)` < 0.05, Estimate, NA))

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

    # max lag
    max.order <- 2

    if (model_system == "cyano") {

        # get the p values of causality
        stats_result <- data %>%
            mutate(order = list(1:max.order)) %>%
            unnest(order) %>%
            group_by(strain, treat, repl, order) %>%
            # for time series data we need to inject NAs between the different sets
            summarise(dat = list(
                data.frame(
                    density = c(density, rep(NA, max.order)),
                    trait = c(trait, rep(NA, max.order))
                )
            )) %>%
            group_by(strain, treat, order) %>%
            summarise(dat = list(bind_rows(dat))) %>%
            rowwise %>%
            mutate(dens_causes_trait = list(
                grangertest(trait ~ density, order = order, data = dat)),
                trait_causes_dens = list(
                    grangertest(density ~ trait, order = order, data = dat))) %>%
            dplyr::select(-dat) %>%
            rowwise %>%
            mutate(dens_causes_trait = na.omit(dens_causes_trait$`Pr(>F)`),
                trait_causes_dens = na.omit(trait_causes_dens$`Pr(>F)`)) %>%
            pivot_longer(cols = dens_causes_trait:trait_causes_dens,
                names_to = 'model', values_to = 'p value')

    } else {

        # get the p values of causality
        stats_result <- data %>%
            mutate(order = list(1:max.order)) %>%
            unnest(order) %>%
            group_by(strain, treat, order) %>%
            summarise(dat = list(
                data.frame(
                    density = density,
                    trait = trait
                )
            )) %>%
            rowwise %>%
            mutate(dens_causes_trait = list(
                grangertest(trait ~ density, order = order, data = dat)),
                trait_causes_dens = list(
                    grangertest(density ~ trait, order = order, data = dat))) %>%
            dplyr::select(-dat) %>%
            rowwise %>%
            mutate(dens_causes_trait = dens_causes_trait$`Pr(>F)`[2],
                         trait_causes_dens = trait_causes_dens$`Pr(>F)`[2]) %>%
            pivot_longer(cols = dens_causes_trait:trait_causes_dens,
                                     names_to = 'model', values_to = 'p value')

    }

    stats_result <- stats_result %>%
        mutate(model = ifelse(grepl('dens_causes', model),
            'Density~causes~trait', 'Trait~causes~density'))

    plot <- ggplot(stats_result, aes(x = as.factor(order), y = `p value`)) +
        theme_bw() +
        scale_colour_manual(values = cbPalette) +
        geom_boxplot() +
        geom_point(aes(col = treat), shape = 16, size = 3, 
            position = position_dodge(width = 0.5)) +
        # critical value = 0.05 / number of models of different lags fitted
        geom_hline(yintercept = 0.05 / max.order, linetype = "dashed") +
        facet_grid(strain ~ model, scales = "free", labeller = label_parsed) +
        scale_y_log10() +
        labs(col = "treatment", x = "lag steps")

    if (model_system == "cyano") {
        ggsave(paste(outpath, "Granger.pdf", sep = ""), 
            plot = plot, width = 5, height = 5.5)
    } else {
        ggsave(paste(outpath, "Granger.pdf", sep = ""), 
            plot = plot, width = 5, height = 8)
    }

    stats_result <- modelling(data = data,
        var_to_nest_by = c("strain", "treat"),
        formulas = c("dT ~ trait",
            "dT ~ trait + density",
            "pcgr ~ density",
            "pcgr ~ trait + density")) %>%
        dplyr::select(-data) %>%
        rowwise %>%
        dplyr::mutate(R2 = abs(summary(model)$r.squared),
            adjR2 = abs(summary(model)$adj.r.squared)) %>%
        mutate(delta_adj = R2 - adjR2,
            R2_max = max(c(R2, adjR2)))

    ggplot(stats_result) +
        scale_shape_manual(values = 0:10) +
        theme_bw() +
        scale_fill_manual(values = cbPalette) +
        aes(x     = form, y = R2_max, fill = treat) +
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

    stats_result <- modelling(data = data,
        var_to_nest_by = c("strain", "treat"),
        formulas = c("dT ~ trait",
            "dT ~ trait + density",
            "pcgr ~ density",
            "pcgr ~ trait + density")) %>%
        dplyr::select(-data) %>%
        rowwise %>%
        dplyr::mutate(R2 = abs(summary(model)$r.squared),
            adjR2 = abs(summary(model)$adj.r.squared)) %>%
        mutate(delta_adj = R2 - adjR2,
            R2_max = max(c(R2, adjR2)))

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
