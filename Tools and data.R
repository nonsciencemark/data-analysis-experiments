library("tidyverse")
library("lubridate")
library("mgcv")

# colour palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7")

# Functions -----------------
modelling <- function(data = data, var_to_nest_by = "strain", formulas) {
    var_to_nest_by <- c(var_to_nest_by, "form")
    test <- expand_grid(data %>% ungroup, form = formulas) %>%
        nest_by(across({{var_to_nest_by}})) %>%
        mutate(model = list(lm(as.formula(form), data = data)))
    return(test)
}

# Import and make uniform the ciliate data ----------
data_cilia <- read_csv("data/ciliates/DIVERCE_TdB_Ciliates_Traits_FULL.csv") %>%
    rename(density = Parts_permL,
        treat = Treatment) %>%
    mutate(strain = as.factor(ID_spec),
        temperature = as.factor(Temp),
        atrazine = as.factor(Atrazine)) %>%
    group_by(atrazine, temperature, strain) %>%
    mutate(pcgr = log(
        lead(density, 1) / density) / (lead(Time_Days, 1) - Time_Days),
        .after = density) %>%
    # remove very negative pcgr and Spiro_5 (not enough data)
    dplyr::filter(pcgr > -5, strain != "Spiro_5") %>%
    ungroup %>%
    # link to species
    separate(strain, into = c("species", "sp.strain"), remove = FALSE) %>%
    # PCA
    group_by(strain, treat) %>%
    nest %>%
    rowwise %>%
    mutate(
        pca_data = list(dplyr::select(data, contains("mean"))),
        pca = list(prcomp(pca_data, center = TRUE, scale = TRUE)),
        data = list(cbind(data, pca1 = pca$x[, 1])),
        pca_varexp = (pca$sdev^2)[1] / sum(pca$sde^2)) %>%
    dplyr::select(strain, treat, data, pca_varexp) %>%
    unnest(data) %>%
    # epxand pca and compute dT and clean up stuff
    rename(trait = pca1) %>%
    group_by(strain, atrazine, temperature) %>%
    mutate(dT = lead(trait, 1) - trait / (lead(Time_Days, 1) - Time_Days)) %>%
    dplyr::filter(!is.na(dT), Temp > 20, Atrazine != 10) %>%
    mutate(treat = case_when(
            (Atrazine == 0) & (Temp == 22) ~ "C",
            # (Atrazine == 10) & (Temp == 22) ~ "a",
            (Atrazine == 20) & (Temp == 22) ~ "A",
            (Atrazine == 0) & (Temp == 24) ~ "T",
            # (Atrazine == 10) & (Temp == 24) ~ "aT",
            (Atrazine == 20) & (Temp == 24) ~ "AT"),
        species = factor(species, levels = c("Loxo", "Spiro", "Tetra", "Para")),
        treat = factor(treat, levels = c("C", "T", "A", "AT")),
        strain = factor(strain, levels = c(
            "Spiro_C", "Spiro_D", "Tetra_1", "Tetra_2",
                "Loxo_1", "Loxo_2", "Para_4")))

# Import and make uniform the cyano data-------------------
data_cyano <- read_csv("data/cyanobacteria/mono_data.csv") %>%
    rename(density = `Population density`, day = date) %>%
    # fixed species names and link to species
    mutate(species = case_when(strain %in% c(2375, 2524) ~ "V", TRUE ~ "VIII"),
        strain = as_factor(strain),
        atrazine = as_factor(case_when(
            treat %in% c("C", "T") ~ "no", TRUE ~ "yes")),
        temperature = as_factor(case_when(
            treat %in% c("A", "C") ~ "normal", TRUE ~ "hot"))) %>%
    group_by(atrazine, temperature, strain) %>%
    ungroup %>%
    # PCA
    group_by(strain, treat) %>%
    nest %>%
    rowwise %>%
    mutate(
        pca_data = list(dplyr::select(data, Size:Phycocyanin)),
        pca = list(prcomp(pca_data, center = TRUE, scale = TRUE)),
        data = list(cbind(data, pca1 = pca$x[, 1])),
        pca_varexp = (pca$sdev^2)[1] / sum(pca$sde^2)) %>%
    dplyr::select(strain, treat, data, pca_varexp) %>%
    unnest(data) %>%
    rename(trait = pca1) %>%
    group_by(strain, atrazine, temperature, repl) %>%
    mutate(dT = lead(trait, 1) - trait / (lead(day, 1) - day)) %>%
    dplyr::filter(!is.na(dT)) %>%
    mutate(treat = factor(treat, levels = c("C", "T", "A", "AT")),
        strain = paste(species, "_", strain, sep = ""))
