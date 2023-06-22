# plotting theme ====
theme_bw2 <- theme_bw() + 
  theme(legend.position = 'bottom', 
        legend.box = 'vertical', 
        legend.margin = margin(0, 5.5, 0, 5.5),
        # aspect.ratio = 1, 
        strip.background = element_rect(fill = 'black'), 
        strip.text = element_text(colour = 'white'),
        axis.ticks.length = unit(-2.75, 'pt'))

# per-capita growth rate, growth phase ====
rate_of_change <- function(exp_data, ...) {
  exp_data <- exp_data %>% 
    mutate(repl = factor(repl)) %>%
    group_by(repl) %>%
    mutate(# time difference
      dt = date - lag(date),
      # lagged values
      lag_pop = lag(pop), 
      lag_FSC = lag(FSC), 
      lag_RED.B = lag(RED.B), 
      lag_RED.R = lag(RED.R), 
      lag_YEL.B = lag(YEL.B),
      # change in values
      dN = log(pop / lag_pop), #(pop - lag_pop) / lag_pop,
      dX_FSC = FSC - lag_FSC,
      dX_RED.B = RED.B - lag_RED.B,
      dX_RED.R = RED.R - lag_RED.R,
      dX_YEL.B = YEL.B - lag_YEL.B,
      # rates of change
      dNdt = dN / dt,
      dX_FSC_dt = dX_FSC / dt,
      dX_RED.B_dt = dX_RED.B / dt,
      dX_RED.R_dt = dX_RED.R / dt,
      dX_YEL.B_dt = dX_YEL.B / dt
    ) %>%
    dplyr::select(-(dN:dX_YEL.B)) %>%
    ungroup
  return(exp_data)
}

# removing outliers ====
filter_outliers <- function(exp_data, ...) {
  # pops and traits
  varnames <- c('dNdt', 'dX_FSC_dt', 'dX_RED.B_dt',  'dX_RED.R_dt',  'dX_YEL.B_dt')
  # one by one
  for (i in varnames) {
    var <- pull(exp_data, i)
    Q1 <- quantile(var, .25, na.rm = TRUE)
    Q3 <- quantile(var, .75, na.rm = TRUE)
    iqr <- IQR(var, na.rm = TRUE)
    outliers <- (var < (Q1 - 1.5 * iqr) | var > (Q3 + 1.5 * iqr))
    var[outliers] <- NA
    exp_data[[i]] <- var
  }
  return(exp_data)
}

# to fit the regressions of whatever model ====
fit_regressions <- function(grow, modnames, ...) {
  # model definitions
  LHS <- 'dNdt'
  RHS <- modnames
  formulae <- paste(LHS, ' ~ ', RHS)
  grow <- grow %>% 
    dplyr::select(dNdt, date, repl, lag_pop, lag_FSC, lag_RED.B, lag_RED.R, lag_YEL.B)
  # model
  mod <- lapply(formulae, function(f) {
    lm.mod <- lm(as.formula(f), grow)
    # remove variables that drive relationship too much
    c.dist <- cooks.distance(lm.mod)
    n <- nrow(lm.mod$model)
    newdat <- lm.mod$model[c.dist < 4/n,]
    newmod <- lm(as.formula(f), newdat)
  })
}

# predicting both growth and recovery ====
model_predict <- function(mods, modnames, grow, reco, ...) {
  grow <- grow %>% 
    dplyr::select(dNdt, date, repl, lag_pop, lag_FSC, lag_RED.B, lag_RED.R, lag_YEL.B)
  reco <- reco %>% 
    dplyr::select(dNdt, date, repl, lag_pop, lag_FSC, lag_RED.B, lag_RED.R, lag_YEL.B)
  # predicting growth phase with growth model
  pred_grow <- lapply(mods, function(i) {
    predict(i, grow)
  })
  grow <- cbind(grow, do.call(what = 'cbind', pred_grow))
  colnames(grow) <- c(colnames(grow)[1:8], names(modnames)) 
  # predicting recovery phase with growth model
  pred_reco <- lapply(mods, function(i) {
    predict(i, reco)
  })
  reco <- cbind(reco, do.call(what = 'cbind', pred_reco))
  colnames(reco) <- c(colnames(reco)[1:8], names(modnames))
  return(list('pred_grow' = grow, 'pred_reco' = reco))
}

# get R squared of both growth and recovery ====
resid_compare <- function(grow, pred_grow, reco, pred_reco, mods, modnames, ...) {
  obs <- c(grow$dNdt, reco$dNdt)
  pred <- lapply(1:length(mods), function(i) {
    c(pred_grow[,8+i], pred_reco[,8+i])
  })
  # model
  predobs_mods <- lapply(pred, function(y) {lm(y ~ obs)})
  # residuals
  mods_res <- lapply(predobs_mods, residuals)
  mods_res <- bind_cols(mods_res)
  colnames(mods_res) <- names(modnames)
  # R2
  mods_R2 <- lapply(predobs_mods, function(i) {summary(i)$adj.r.squared})
  return(list(predobs_mods, mods_res, mods_R2))
}


# modelling functions from FDL ====

## lm method ====
modelling <- function(data=data, var_to_nest_by="strain", 
                      response="pcgr", predictor="pop") {
  test <- data %>%
    ungroup() %>%
    group_by(across({{var_to_nest_by}})) %>%
    nest() %>%
    rename(data_to_model = data) %>%
    expand(nesting(data_to_model), predictor) %>%
    #expand(nesting({{var_to_nest_by}}, data_to_model), predictor) %>%
    ungroup() %>%
    mutate(formula = paste0(response, "~`", predictor, '`')) %>%
    rowwise() %>%
    mutate(model = list(lm(as.formula(formula), data = data_to_model))) %>%
    dplyr::select(-c("data_to_model", "formula")) 
  return(test)
}

## pcr method ====
modelling_PC <- function(data=data, var_to_nest_by="strain", 
                         response="pcgr", predictor='Size+Chlorophyll+Phycocyanin+Phycoerythrin') {
  test <- data %>%
    ungroup() %>%
    group_by(across({{var_to_nest_by}})) %>%
    nest() %>%
    rename(data_to_model = data) %>%
    expand(nesting(data_to_model), predictor) %>%
    #expand(nesting({{var_to_nest_by}}, data_to_model), predictor) %>%
    ungroup() %>%
    mutate(formula = paste0(response, "~", predictor)) %>%
    rowwise() %>%
    mutate(model = list(pcr(as.formula(formula), data = data_to_model))) %>%
    dplyr::select(-c("data_to_model", "formula")) 
  return(test)
}
