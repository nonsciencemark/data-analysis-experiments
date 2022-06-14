# parallel setup
library('tidyverse')
library('parallel')
library('multidplyr')

# initialise clusters
ncores <- parallel::detectCores() - 1 # i want to be able to use my computer while it processes
cluster <- new_cluster(ncores)
cluster # init cluster

# send them relevant packages
cluster_library(cluster, "tidyverse")
cluster_library(cluster, "flowCore")
cluster_library(cluster, "flowDensity")
cluster_library(cluster, "flowClust")
cluster_library(cluster, "cyanoFilter")
cluster_library(cluster, "gauseR")
cluster_library(cluster, "PeacoQC")

# source useful R files
source('C:/Users/mholmes/Documents UNamur/Personnel/flowClust_cluster.R')
source('C:/Users/mholmes/Documents UNamur/Personnel/flowCore_get_datasets.R')
source('C:/Users/mholmes/Documents UNamur/Personnel/CytoTree_ExprsMerge.R')
cluster_send(cluster, source('C:/Users/mholmes/Documents UNamur/Personnel/flowClust_cluster.R'))
cluster_send(cluster, source('C:/Users/mholmes/Documents UNamur/Personnel/flowCore_get_datasets.R'))

# send cytometer varnames to clusters
varnames <- c("FSC.HLin","SSC.HLin","GRN.B.HLin","YEL.B.HLin","RED.B.HLin",
              "NIR.B.HLin","RED.R.HLin","NIR.R.HLin")
varmins <- c(-1,1,-10,-10,2,-10,-10,-10)
cluster_assign(
  cluster, 
  varnames = c("FSC.HLin","SSC.HLin","GRN.B.HLin","YEL.B.HLin","RED.B.HLin",
               "NIR.B.HLin","RED.R.HLin","NIR.R.HLin"),
  varmins = c(-1,1,-10,-10,2,-10,-10,-10))

# create structure of well plates ====

# plate well names
wellnames <- paste0(rep(LETTERS[1:8], each = 12), 
                    c('01', '02', '03', '04', '05', '06', 
                      '07', '08', '09', '10', '11', '12'))

# monoculture
mono.order <- c('2375', '2383', '2434', '2524') # strain order in plates

mono.id <- data.frame(well = wellnames[1:48], # create 96-well plate structure
                      strain = mono.order, 
                      treat = rep(c('C', 'A', 'T', 'AT'), each = 12),
                      repl = rep(1:3, each = 4))

# NOTE THAT THERE MAY BE LOST DATA ON 2022-01-25, WELLS B10-12, 
# BECAUSE THE MEASURED DENSITY WAS ABOUT 39 CELLS/UL

# get files ====

# file location
# get wdirs ====
wdirs <- list.dirs('C:/Users/mholmes/Documents UNamur/Lab work/03_Experiments/12_2022-01-14_DIVERCE//', 
                   full.names = T, recursive = T)
mono.wdir <- wdirs[grepl('MONO', wdirs)] # get mono subdirs

# monoculture processing ====
mono.data <- list.files(# list all fcs files
  mono.wdir, recursive = TRUE, full.names = TRUE, pattern = '\\.fcs') %>% 
  lapply(function(i) { # get info from files
    info <- getdatasets(i) %>% dplyr::filter(well %in% mono.id$well)}) %>% 
  bind_rows %>% as_tibble %>% # put into tibble
  full_join(mono.id, by = 'well') %>% # add the info from the mono.id table
  na.omit %>% # remove any bad rows
  mutate(date.time = as.POSIXct(date.time, format = '%d-%b-%Y_%H:%M:%S'),
         Date = as.POSIXct(as.Date(date.time))) %>% # reformat time
  group_by(strain, treat, repl, Date) %>% # define grouping variables
  dplyr::filter(date.time == max(date.time)) %>% # only take the last sample measured (i.e. correct one) for each well
  rowwise %>% # rowwise operations
  dplyr::filter(Date < as.POSIXct('2022-01-21 23:00:00'),
                treat == 'T' | treat == 'AT') %>% 
  mutate(treat = recode(treat, 'T' = 'Control', 'AT' = 'Atrazine')) %>%
  partition(cluster) %>% # split across clusters
  mutate(
    fcs = list( 
      read.FCS(filename = filename, # load fcs files into tibble as list-col
               dataset = n, 
               emptyValue = FALSE, 
               alter.names = TRUE,
               truncate_max_range = FALSE) %>% 
        noNeg %>% noNA %>% lnTrans %>% # remove bad data and log-transform
        RemoveMargins(channels = 1:10) %>% # remove margin events and doublets
        RemoveDoublets(channel1 = 'SSC.ALin', channel2 = 'SSC.HLin', nmad = 2) %>%
        exprs), # don't need fcs metadata really
    fcs = list(fcs[,1:8]), # don't need later columns anymore
    clusts = list(flowClust(fcs, varNames = varnames, K = 2, min = varmins)), # cluster to separate debris
    debris.col = which.min(rowMeans(clusts@mu)), # debris has lower values
    cells.col = (1:2)[-debris.col], # cells other col
    cells = colSums(clusts@z, na.rm = T)[cells.col], # count cells
    debris = colSums(clusts@z, na.rm = T)[debris.col], # count debris
    uncertainty = list(clusts@uncertainty), # record uncertainty for each cell
    weights = list(clusts@z[,cells.col] %>% replace_na(0)), # get cell weighting for multiplying with variables later
    channel.means = list(apply(fcs, 2, weighted.mean, w = weights)), # get means weighted by likelihood of cell
    channel.vars = list(apply(fcs, 2, weighted.var, w = weights)), # get vars weighted by likelihood of cell
    pop = cells / vol * dil, # calculate cell concentration
    debris = debris / vol * dil, # calculate debris concentration
    ratio = pop / debris) %>% # pop:debris ratio 
  collect %>% # get all data frames from different clusters
  arrange(strain, treat, repl, Date) %>% # order tibble
  dplyr::select(date.time, strain, treat, repl, Date, # keep only useful columns
                fcs, clusts, weights, cells.col, debris.col,
                uncertainty, channel.means, channel.vars, pop, debris, ratio) %>%
  group_by(Date, strain, treat, repl) # re-group data

varsDf <- mono.data %>%
  #dplyr::select(strain, treat, repl, Date, fcs, clusts, cells.col) %>%
  rowwise %>%
  mutate(fcs = list(data.frame(fcs, 'lab' = clusts@label) %>% 
                       dplyr::filter(lab == cells.col) %>%
                       dplyr::select(-lab)))

fcsVars <- varsDf %>% 
  ungroup %>%
  dplyr::select(strain, treat, repl, Date, fcs) %>%
  unnest(fcs)
  
PCA <- prcomp(fcsVars %>% dplyr::select(FSC.HLin:RED.R.HLin),
              center = TRUE, scale. = TRUE, retx = TRUE)

loadings <- summary(PCA)$importance[2,]
rLoadings <- round(loadings, 3)

varsDf2 <- fcsVars %>%
  ungroup %>%
  mutate(PC1 = PCA$x[,'PC1'],
         PC2 = PCA$x[,'PC2']) %>%
  group_by(strain, treat, Date) %>% 
  summarise(across(FSC.HLin:PC2, list(mean = mean, sd = sd))) %>%
  mutate(treat = factor(treat, levels = c('Control', 'Atrazine')))

# create arrows
arrows <- PCA$rotation[,1:2] %>%
  as.data.frame() %>%
  rownames_to_column(var = 'variable') %>%
  mutate(origin = 0,
         end = 1) %>%
  pivot_longer(origin:end) %>%
  mutate(PC1 = PC1 * value,
         PC2 = PC2 * value)

ggpubr::ggarrange(
  
  ggplot(varsDf2, aes(x = PC1_mean, y = PC2_mean, col = strain, fill = Date)) +
    facet_wrap(.~treat) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_path() +
    geom_errorbar(aes(ymin = PC2_mean - PC2_sd, ymax = PC2_mean + PC2_sd), alpha = .25) + 
    geom_errorbarh(aes(xmin = PC1_mean - PC1_sd, xmax = PC1_mean + PC1_sd), alpha = .25) + 
    geom_point(pch = 21, stroke = 1.1, size = 2) +
    scale_fill_viridis_c(trans = scales::time_trans()) +
    scale_color_brewer(palette = 'Dark2') +
    theme_bw() +
    labs(x = paste0('PC1 (', rLoadings[1], ')'), 
         y = paste0('PC2 (', rLoadings[2], ')'), 
         col = 'Strain', title = 'Trait changes over experiment', 
         subtitle = 'Points indicate mean values and error bars indicate SD.'),
  
  ggplot() +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_path(data = arrows, inherit.aes = FALSE, 
              arrow = arrow(length = unit(0.1, 'inches')),
              mapping = aes(x = PC1 * 2, y = PC2 * 2, group = variable)) +
    ggrepel::geom_text_repel(data = arrows %>% dplyr::filter(name == 'end'), inherit.aes = FALSE, 
                             mapping = aes(x = PC1 * 2, y = PC2 * 2, label = variable)) +
    theme_bw() +
    scale_x_continuous(limits = range(varsDf2$PC1_mean * 1.05)) +
    scale_y_continuous(limits = range(varsDf2$PC2_mean * 1.05)) +
    labs(x = paste0('PC1 (', rLoadings[1], ')'), 
         y = paste0('PC2 (', rLoadings[2], ')'), 
         title = 'PCA loadings of different channels',
         subtitle = paste0('Total explained = ', sum(loadings[1:2]) * 100, '%')),
  
  ncol = 2, widths = c(2,1)
)

# alpha stuff

traitCompare <- function(y1, y2) {
  exp(-(y1 - y2)^2)
}

## only relevant traits
usefulChannels <- c('FSC.HLin', 'SSC.HLin', 'RED.R.HLin', 'RED.B.HLin', 'YEL.B.HLin')

# GET PARAMETERS AND RELEVENAT CHANNELS AND DO THE TRAIT COMPARE FUNCTION AND STUFF

alphaSelf <- mono.data %>% 
  dplyr::select(strain, treat, repl, Date, usefulChannels)

