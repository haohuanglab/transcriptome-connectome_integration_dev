library(nlme)
library(mgcv)
library(knitr)
library(readr)
library(openxlsx)

rm(list=ls())
gc()

folder_thesedata <- "/Data"

filename_agetable_tofit <- "age_table_tofit.txt"
fn_agetable_tofit <- file.path(folder_thesedata, filename_agetable_tofit)

output_folder <- file.path(folder_thesedata, "Connectomic")

network_folder <- file.path(output_folder, "network")  

filename_aic_table <- file.path(output_folder, "aic_table.xlsx")

nls_funcs <- c("exponential", "poisson", "sigmoid")

recode_cohort <- function(x) {
  factor(as.character(x))
}

gof_model <- function(age_vec,
                      tmetric_vec,
                      cohort_vec,
                      flag_model,
                      flag_transfxBeforefit,
                      param_model,
                      age_vec_topred = NULL,
                      df_xtick = NULL,
                      unit_age = NULL,
                      GRETNA_metric = NULL,
                      y_max = NULL,
                      flag_xaxis = "plot",
                      flag_yaxis = "plot",
                      return_plot = TRUE) {

  num_age_vec_topred <- 200

  if (is.null(age_vec_topred)) {
    age_vec_topred <- seq(min(age_vec), max(age_vec), length.out = num_age_vec_topred)
  }

  if (is.null(df_xtick)) {
    df_xtick <- data.frame(
      age = 38 * 7 / 30.5 + 12 * c(seq(0, 4, by = 1), 6, 8, 15, 22),
      tick = c("Birth", "1", "2", "3", "4", "6", "8", "15  ", "  22")
    )
  }

  if (is.null(tmetric_vec)) {
    stop("tmetric_vec is NULL!")
  }

  cohort <- recode_cohort(cohort_vec)

  if (flag_transfxBeforefit == "none") {
    x <- age_vec
    x_topred <- age_vec_topred
  } else if (flag_transfxBeforefit == "log10") {
    x <- log10(age_vec)
    x_topred <- log10(age_vec_topred)
  } else if (flag_transfxBeforefit == "loge") {
    x <- log(age_vec, base = exp(1))
    x_topred <- log(age_vec_topred, base = exp(1))
  } else {
    stop("Invalid flag_transfxBeforefit")
  }

  y <- tmetric_vec
  df <- data.frame(x = x, y = y, cohort = cohort)

  ref_cohort <- levels(df$cohort)[1]
  df_topred <- data.frame(
    x = x_topred,
    cohort = factor(ref_cohort, levels = levels(df$cohort))
  )

  if (flag_model == "gam") {
    m1 <- gam(
      y ~ s(x, k = param_model$k_gam, bs = param_model$bs_gam) +
        s(cohort, bs = "re"),
      data = df,
      method = param_model$method_gam
    )
    pred <- predict.gam(
      m1,
      newdata = df_topred,
      se.fit = TRUE,
      exclude = "s(cohort)"
    )
    y_pred <- as.numeric(pred$fit)
    y_pred_se <- as.numeric(pred$se.fit)
  } else if (flag_model == "linear") {
    m1 <- lme(
      fixed = y ~ x,
      random = ~1 | cohort,
      data = df,
      method = "ML"
    )
    y_pred <- as.numeric(predict(m1, newdata = df_topred, level = 0))
    y_pred_se <- rep(NA_real_, length(y_pred))
  } else if (flag_model == "quadratic") {
    df$x.sq <- df$x^2
    df_topred$x.sq <- df_topred$x^2
    m1 <- lme(
      fixed = y ~ x + x.sq,
      random = ~1 | cohort,
      data = df,
      method = "ML"
    )
    y_pred <- as.numeric(predict(m1, newdata = df_topred, level = 0))
    y_pred_se <- rep(NA_real_, length(y_pred))
  } else if (flag_model == "log") {
    df$logx <- log(df$x, base = exp(1))
    df_topred$logx <- log(df_topred$x, base = exp(1))
    m1 <- lme(
      fixed = y ~ logx,
      random = ~1 | cohort,
      data = df,
      method = "ML"
    )
    y_pred <- as.numeric(predict(m1, newdata = df_topred, level = 0))
    y_pred_se <- rep(NA_real_, length(y_pred))
  } else if (flag_model == "inverse") {
    df$x.inv <- 1 / df$x
    df_topred$x.inv <- 1 / df_topred$x
    m1 <- lme(
      fixed = y ~ x.inv,
      random = ~1 | cohort,
      data = df,
      method = "ML"
    )
    y_pred <- as.numeric(predict(m1, newdata = df_topred, level = 0))
    y_pred_se <- rep(NA_real_, length(y_pred))
  } else if (flag_model == "sqrt") {
    df$sqrtx <- sqrt(df$x)
    df_topred$sqrtx <- sqrt(df_topred$x)
    m1 <- lme(
      fixed = y ~ sqrtx,
      random = ~1 | cohort,
      data = df,
      method = "ML"
    )
    y_pred <- as.numeric(predict(m1, newdata = df_topred, level = 0))
    y_pred_se <- rep(NA_real_, length(y_pred))
  } else if (flag_model == "cubic") {
    df$x.sq <- df$x^2
    df$x.cu <- df$x^3
    df_topred$x.sq <- df_topred$x^2
    df_topred$x.cu <- df_topred$x^3
    m1 <- lme(
      fixed = y ~ x + x.sq + x.cu,
      random = ~1 | cohort,
      data = df,
      method = "ML"
    )
    y_pred <- as.numeric(predict(m1, newdata = df_topred, level = 0))
    y_pred_se <- rep(NA_real_, length(y_pred))
  } else if (flag_model == "arctan") {
    df$atanx <- atan(df$x)
    df_topred$atanx <- atan(df_topred$x)
    m1 <- lme(
      fixed = y ~ atanx,
      random = ~1 | cohort,
      data = df,
      method = "ML"
    )
    y_pred <- as.numeric(predict(m1, newdata = df_topred, level = 0))
    y_pred_se <- rep(NA_real_, length(y_pred))
  } else if (flag_model == "tanh") {
    df$tanhx <- tanh(df$x)
    df_topred$tanhx <- tanh(df_topred$x)
    m1 <- lme(
      fixed = y ~ tanhx,
      random = ~1 | cohort,
      data = df,
      method = "ML"
    )
    y_pred <- as.numeric(predict(m1, newdata = df_topred, level = 0))
    y_pred_se <- rep(NA_real_, length(y_pred))
  } else if (flag_model == "power") {
    pooled_nls <- nls(
      y ~ a * x^b + c,
      data = df,
      start = list(a = 1, b = 1, c = 0),
      control = nls.control(maxiter = 500, warnOnly = TRUE)
    )
    start_vals <- coef(pooled_nls)
    m1 <- nlme(
      y ~ a * x^b + c,
      data = df,
      fixed = a + b + c ~ 1,
      random = c ~ 1 | cohort,
      start = start_vals,
      control = nlmeControl(
        maxIter = 500,
        pnlsMaxIter = 50,
        msMaxIter = 500,
        returnObject = TRUE
      )
    )
    y_pred <- as.numeric(predict(m1, newdata = df_topred, level = 0))
    y_pred_se <- rep(NA_real_, length(y_pred))
  } else if (flag_model == "exponential") {
    pooled_nls <- nls(
      y ~ a * exp(-b * x) + c,
      data = df,
      start = list(a = 100, b = 1, c = 0),
      control = nls.control(
        maxiter = 500,
        warnOnly = TRUE,
        minFactor = 1/4096
      )
    )

    start_vals <- coef(pooled_nls)

    m1 <- nlme(
      y ~ a * exp(-b * x) + c,
      data = df,
      fixed = a + b + c ~ 1,
      random = c ~ 1 | cohort,
      start = start_vals,
      control = nlmeControl(
        maxIter = 500,
        pnlsMaxIter = 50,
        msMaxIter = 500,
        minScale = 1e-8,
        tolerance = 1e-6,
        pnlsTol = 1e-4,
        returnObject = TRUE
      )
    )

    y_pred <- as.numeric(predict(m1, newdata = df_topred, level = 0))
    y_pred_se <- rep(NA_real_, length(y_pred))
  } else if (flag_model == "poisson") {
    pooled_nls <- nls(
      y ~ a * exp(-b * x) + c,
      data = df,
      start = list(a = exp(4.56), b = -0.7528, c = 0),
      control = nls.control(
        maxiter = 500,
        warnOnly = TRUE,
        minFactor = 1/4096
      )
    )
    start_vals <- coef(pooled_nls)
    m1 <- nlme(
      y ~ a * exp(-b * x) + c,
      data = df,
      fixed = a + b + c ~ 1,
      random = c ~ 1 | cohort,
      start = start_vals,
      control = nlmeControl(
        maxIter = 500,
        pnlsMaxIter = 50,
        msMaxIter = 500,
        returnObject = TRUE
      )
    )
    y_pred <- as.numeric(predict(m1, newdata = df_topred, level = 0))
    y_pred_se <- rep(NA_real_, length(y_pred))
  } else if (flag_model == "sigmoid") {
    pooled_nls <- nls(
      y ~ SSlogis(x, Asym, xmid, scal),
      data = df
    )

    start_vals <- coef(pooled_nls)

    m1 <- nlme(
      y ~ SSlogis(x, Asym, xmid, scal),
      data = df,
      fixed = Asym + xmid + scal ~ 1,
      random = Asym ~ 1 | cohort,
      start = start_vals,
      control = nlmeControl(
        maxIter = 500,
        pnlsMaxIter = 50,
        msMaxIter = 500,
        returnObject = TRUE
      )
    )

    y_pred <- as.numeric(predict(m1, newdata = df_topred, level = 0))
    y_pred_se <- rep(NA_real_, length(y_pred))
  } else {
    stop(paste("Invalid flag_model:", flag_model))
  }

  toreturn <- list(
    m1 = m1,
    AIC = AIC(m1),
    BIC = BIC(m1),
    y_pred = y_pred,
    y_pred_se = y_pred_se
  )

  toreturn$r2 <- NA_real_

  return(toreturn)
}

AIC_metric <- function(df) {
  df$delta_AIC <- df$AIC - min(df$AIC)
  df$likelihood_AIC <- exp(-0.5 * df$delta_AIC)
  df$wAIC <- df$likelihood_AIC / sum(df$likelihood_AIC)
  return(df)
}

BIC_metric <- function(df) {
  df$delta_BIC <- df$BIC - min(df$BIC)
  return(df)
}

flag_splitdataset <- "full" 

metrics <- c("Ne", "Lp", "Neloc") #"Bc", Cp", "Dc"

list_11NCX_regions <- c(
  "MFC", "OFC", "DFC", "VFC", "M1C", "S1C",
  "IPC", "A1C", "STC", "ITC", "V1C"
)

models_to_use <- c("gam", "exponential", "gam")

wb <- createWorkbook()

bs_gam <- "cs"
method_gam <- "REML"
k_gam <- 3
param_model <- data.frame(k_gam, bs_gam, method_gam)

df_full <- list()

flag_savefit <- TRUE
flag_savefit_overwrite <- FALSE

flag_calcu_AICmetrics <- 1

for (i in seq_along(metrics)) {

  GRETNA_metric <- metrics[i]
  model_to_plot <- models_to_use[i]

  region <- ifelse(
    GRETNA_metric == "NodalStrength",
    "avgNodalStrength",
    "avgMetric"
  )

  message("Working on ", GRETNA_metric, " ...")

  fn_tmetric <- file.path(
    network_folder,
    paste0(GRETNA_metric, "_ROI_recoded.txt")
  )

  tmetric_everything <- read.table(fn_tmetric, header = TRUE)

  if ("cohort" %in% names(tmetric_everything)) {
    cohort_col <- "cohort"
  } else if ("whichdataset" %in% names(tmetric_everything)) {
    cohort_col <- "whichdataset"
  } else {
    stop(
      "Neither 'cohort' nor 'whichdataset' column found in ",
      basename(fn_tmetric), " -- cannot determine cohort assignment."
    )
  }

  tmetric_everything$cohort_model <- recode_cohort(
    tmetric_everything[[cohort_col]]
  )

  tmetric_everything$agepcy <- tmetric_everything$agepcd / 365
  tmetric_everything$agepcm <- tmetric_everything$agepcd / 30.5

  df_xtick_full <- data.frame(
    agepcm = 38 * 7 / 30.5 +
      12 * c(seq(0, 4, by = 1), 6, 8, 15, 22),
    tick = c("Birth", "1", "2", "3", "4", "6", "8", "15  ", "  22")
  )

  agetable_tofit <- read.table(fn_agetable_tofit, header = TRUE)
  agetable_tofit$pcm <- agetable_tofit$pcd / 30.5

  f_pred_list <- list()

  models_to_try <- c(
    "linear", "quadratic", "log", "sqrt",
    "gam", "power", "inverse", "cubic",
    "arctan", "tanh", nls_funcs
  )

  index_of_model_to_plot <- which(models_to_try == model_to_plot)

  df_list <- list()

  df_pred <- data.frame(
    matrix(
      data = NA,
      nrow = nrow(agetable_tofit),
      ncol = 2 + length(list_11NCX_regions)
    )
  )

  colnames(df_pred) <- c(
    "agepcd",
    as.character(c(region, list_11NCX_regions))
  )

  df_pred$agepcd <- agetable_tofit$pcd

  for (i_model in seq_along(models_to_try)) {
    model_to_try <- models_to_try[i_model]
    is_plotting <- model_to_try == model_to_plot

    df <- data.frame(
      region = c(region, list_11NCX_regions),
      model = rep(NA, length(list_11NCX_regions) + 1),
      r2 = rep(NA, length(list_11NCX_regions) + 1),
      AIC = rep(NA, length(list_11NCX_regions) + 1),
      BIC = rep(NA, length(list_11NCX_regions) + 1)
    )

    for (i_region in seq_len(nrow(df))) {
      ROIstr <- as.character(df$region[i_region])
      temp <- tmetric_everything[[ROIstr]]

      output <- tryCatch(
        gof_model(
          age_vec = tmetric_everything$agepcm,
          age_vec_topred = agetable_tofit$pcm,
          unit_age = "pcm",
          tmetric_vec = temp,
          cohort_vec = tmetric_everything$cohort_model,
          flag_model = model_to_try,
          flag_transfxBeforefit = "log10",
          param_model = param_model,
          df_xtick = data.frame(
            age = df_xtick_full$agepcm,
            tick = df_xtick_full$tick
          ),
          GRETNA_metric = GRETNA_metric,
          return_plot = is_plotting
        ),
        error = function(e) {
          message(
            "Model failed: metric=", GRETNA_metric,
            ", region=", ROIstr,
            ", model=", model_to_try,
            " | ", e$message
          )
          NULL
        }
      )

      if (!is.null(output)) {
        df$AIC[i_region] <- output$AIC
        df$BIC[i_region] <- output$BIC
        df$r2[i_region] <- output$r2
        df$model[i_region] <- model_to_try

        if (is_plotting) {
          df_pred[[ROIstr]] <- output$y_pred
        }
      }
    }

    df_list[[i_model]] <- df
  }

  df_full[[GRETNA_metric]] <- df_list[[index_of_model_to_plot]]

  if (flag_calcu_AICmetrics == 1) {

    df_aic <- data.frame(name = models_to_try)

    df_aic$AIC <- unlist(
      lapply(
        df_list,
        function(x) {
          ifelse(
            sum(is.na(x$AIC)) > 0,
            1.2e100,
            median(x$AIC)
          )
        }
      )
    )

    df_aic$BIC <- unlist(
      lapply(
        df_list,
        function(x) {
          ifelse(
            sum(is.na(x$BIC)) > 0,
            1.2e100,
            median(x$BIC)
          )
        }
      )
    )

    df_aic[is.na(df_aic)] <- 1.2e100

    df_aic <- AIC_metric(df_aic)
    df_aic <- BIC_metric(df_aic)

    df_aic$delta_comb <- df_aic$delta_AIC + df_aic$delta_BIC

    df_aic <- df_aic[
      order(
        df_aic$delta_comb,
        df_aic$delta_BIC,
        df_aic$delta_AIC
      ),
    ]

    print(knitr::kable(df_aic))

    addWorksheet(wb, GRETNA_metric)
    writeData(wb, sheet = GRETNA_metric, x = df_aic)
  }

  if (flag_savefit) {
    df_pred_tosave <- df_pred[, c("agepcd", list_11NCX_regions)]
    names(df_pred_tosave)[names(df_pred_tosave) == "agepcd"] <- "agepcd_todo"

    fit_output_file <- file.path(
      output_folder,
      paste0(
        GRETNA_metric,
        "_ROIbest_fit_agefit.txt"
      )
    )

    if (!file.exists(fit_output_file) || flag_savefit_overwrite) {
      write.table(
        df_pred_tosave,
        fit_output_file,
        row.names = FALSE
      )
    } else {
      message(
        "Did not overwrite existing file: ",
        fit_output_file
      )
    }
  }

}

if (flag_calcu_AICmetrics == 1) {
  saveWorkbook(
    wb,
    filename_aic_table,
    overwrite = TRUE
  )
}

message("Done. Outputs saved to: ", output_folder)
