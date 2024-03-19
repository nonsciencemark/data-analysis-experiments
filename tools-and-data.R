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
        data = list(cbind(data, pca1 = pca$x[, 1:2])),
        pca_varexp = sum((pca$sdev^2)[1:2]) / sum(pca$sde^2)) %>%
    dplyr::select(strain, treat, data, pca_varexp) %>%
    unnest(data) %>%
    # epxand pca and compute dT and clean up stuff
    rename(trait_1 = pca1.PC1, trait_2 = pca1.PC2) %>%
    group_by(strain, atrazine, temperature) %>%
    mutate(dT1 = lead(trait_1, 1) - trait_1 / (lead(Time_Days, 1) - Time_Days),
        dT2 = lead(trait_2, 1) - trait_2 / (lead(Time_Days, 1) - Time_Days)) %>%
    dplyr::filter(!is.na(dT1), Temp > 20, Atrazine != 10) %>%
    mutate(treat = case_when(
            (Atrazine == 0) & (Temp == 22) ~ "C",
            (Atrazine == 20) & (Temp == 22) ~ "A",
            (Atrazine == 0) & (Temp == 24) ~ "T",
            (Atrazine == 20) & (Temp == 24) ~ "AT"),
        species = factor(species, levels = c("Loxo", "Spiro", "Tetra", "Para")),
        treat = factor(treat, levels = c("C", "T", "A", "AT")),
        strain = factor(strain, levels = c(
            "Spiro_C", "Spiro_D", "Tetra_1", "Tetra_2",
                "Loxo_1", "Loxo_2", "Para_4")))

cilia_trait_levels <- c("mean_area", "mean_speed", "mean_ar",
    "mean_linearity", "trait_1", "trait_2")
cilia_trait_labels <- c("Size", "Speed", "Aspect~ratio",
    "Linearity", "Trait 1", "Trait 2")
# strain_labels <- c("Spiro[1]", "Spiro[2]", "Tetra[1]", "Tetra[2]", "Loxo[1]", "Loxo[2]", "Para[4]")

long_data_cilia <- data_cilia %>%
    mutate(density = log10(density)) %>%
    pivot_longer(all_of(c("density", cilia_trait_levels))) %>%
    mutate(
        name = factor(name,
            levels = c("density", cilia_trait_levels),
            labels = c("log[10]~density", cilia_trait_labels)),
        # strain = factor(strain, levels = levels(strain), labels = strain_labels),
        day = Time_Days)

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

cyano_trait_levels <- c("Size", "Chlorophyll", "Phycocyanin",
    "Phycoerythrin", "trait")
cyano_trait_labels <- c("Size", "Chlorophyll", "Phycocyanin",
    "Phycoerythrin", "Trait")

long_data_cyano <- data_cyano %>%
    mutate(density = log10(density)) %>%
    pivot_longer(all_of(c("density", cyano_trait_levels))) %>%
    mutate(name = factor(name,
        levels = c("density", cyano_trait_levels),
        labels = c("log[10]~density", cyano_trait_labels)))

# merged data
data_merged <- bind_rows(
    data_cilia %>% mutate(system = "Ciliates"),
    data_cyano %>% mutate(system = "Cyanobacteria")
)

long_data_merged <- bind_rows(
    long_data_cilia %>% mutate(system = "Ciliates"),
    long_data_cyano %>% mutate(system = "Cyanobacteria")
)
