# these are the functions for 2_Model_Structure_for_Baseline_and_Dynamic_Model

# libraries
library(survival)

# keep results
init_preds_BASE_DYN <- function(){
  
  results <- tibble(preds = double(),
                    y_true_cat = character(),
                    y_true_time = double(),
                    train_set = character(),
                    test_set = character(),
                    model = character(),
                    model_type = character(),
                    horizon = character(),
                    imputation = character(),
                    LM = double(),
                    functioneelDossierNr = double(),
                    CAT_catheter_episode = double())
  
  return(results) 
}

init_results_BASE_DYN <- function(){
  
  results <- tibble(train_set = character(),
                    test_set = character(),
                    metric = character(),
                    value = numeric(),
                    model = character(),
                    model_type = character(),
                    horizon = character(),
                    imputation = character(),
                    LM = double())
  
  return(results)
}

init_coefs <- function(){
  
  results <- tibble(variable = character(),
                    value = numeric(),
                    train_set = character(),
                    model = character(),
                    model_type = character(),
                    horizon = character(),
                    imputation = character(),
                    LM = double())

  
  return(results)
}

init_results_sumstats <- function(){
  
  summary_results <- tibble(median_AUC = character(),
                            mean_AUC = character(),
                            median_Calibration_Slope = character(),
                            mean_Calibration_Slope = character(),
                            median_Calibration_Intercept = character(),
                            mean_Calibration_Intercept = character(),
                            median_OE_ratio = character(),
                            mean_OE_ratio = character(),
                            median_ECE = character(),
                            mean_ECE = character(),
                            median_ECI = character(),
                            mean_ECI = character(),
                            median_BS = character(),
                            mean_BS = character(),
                            median_Scaled_BS = character(),
                            mean_Scaled_BS = character(),
                            model = character(),
                            imputation = character())
  
  return(summary_results)
}

init_results_plot <- function(){
  
  summary_plot <- tibble(metric = character(),
                         mean = numeric(),
                         lower.ci = numeric(),
                         upper.ci = numeric(),
                         median = numeric(),
                         min = numeric(),
                         max = numeric(), 
                         Q1 = numeric(),
                         Q3 = numeric(),
                         model = character(),
                         imputation = character())
  
  return(summary_plot)
}

# for computing performances in terms of various calibration and discrimination measures
# using only y_test and pred_test for calculation, treating all observations/predictions as binary

#' @param dd test dataset
#' @param fm The formula that will be called by the model, of the form \code{y_test ~ pred_test}
#'  
#' @describeIn Area under the ROC curve
#' @importFrom pROC roc auc
#' @export
c_statistic <- function(obs, pred){
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package \"pROC\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  m <- pROC:::roc(obs, pred, direction="<", quiet = TRUE)
  pROC:::auc(m)
}

#' @describeIn calibration_slope Estimate calibration slope
#' @export
calibration_slope <- function(obs, pred){
  dat <- data.frame(e = pred, o = obs)
  dat$e[dat$e == 0] = 0.0000000001
  dat$e[dat$e == 1] = 0.9999999999
  fm <- as.formula(o~e)
  fm <- update(fm, . ~ qlogis(.))
  m <- glm(fm, dat, family = binomial, na.action = na.omit)
  coef(m)[2]
}


#' @describeIn O/E ratio Estimate ratio of observed to expected number of events
#' @export
oe_ratio <- function(obs, pred){
  sum(obs) / sum(pred)
}


#' @describeIn ECI Estimate estimated calibration index (average squared difference)
#' @export
ECI <- function(obs, pred){
  idx <- order(pred)
  obs <- obs[idx]
  pred <- pred[idx]
  argzLoess = alist(degree = 2)
  # pred <- jitter(pred, factor=0.2)
  argzLoess$formula = obs ~ pred
  Sm <- do.call("loess", argzLoess)
  Sm <- data.frame(Sm$x, Sm$fitted)
  cal.smooth <- approx(Sm, xout = pred, ties = "ordered")$y
  coom <- cbind(pred, cal.smooth)
  mean((pred - cal.smooth) ^ 2) * 100
}


#' @describeIn scaled brier score Estimate Scaled Brier Score
#' @export
scaled_brier_score <- function(obs, pred) {
  1 - (brier_score(obs, pred) / (mean(obs, na.rm = TRUE) * (1 - mean(obs, na.rm = TRUE))))
}



# apply the same sequence of thresholds to all samples and take the mean of the sensitivity and specificity per threshold to get the "mean ROC curve"

mean_roc <- function(data, cutoffs = seq(from = 0, to = 1, by = 0.01)) {
  map_df(cutoffs, function(cp) {
    out <- cutpointr(data = data, x = preds, class = y_test,
                     subgroup = test_set, method = oc_manual, cutpoint = cp,
                     pos_class = 1, direction = ">=")
    data.frame(cutoff = cp, 
               sensitivity = mean(out$sensitivity),
               specificity = mean(out$specificity))
  })
}


# for making the 10 deciles of calibration plot

lift_table <- function(table, bin_number = 10) {
  table %>%
    dplyr::mutate(yhat_bin = ggplot2::cut_number(yhat, bin_number)) %>%
    dplyr::group_by(yhat_bin) %>%
    dplyr::summarise(mean_y = mean(y), mean_yhat = mean(yhat)) 
}

# for adding rcs smooth calibration curves with confidence bands
# this comes from https://github.com/BavoDC/CalibrationCurves/blob/master/R/rcspline.plot.noprint.R

#' Internal function
#'
#' Adjusted version of the \code{\link[Hmisc]{rcspline.plot}} function where only the output is returned and no plot is made
#'
#'
#' @param x a numeric predictor
#' @param y a numeric response. For binary logistic regression, \code{y} should be either 0 or 1.
#' @param model \code{"logistic"} or \code{"cox"}. For \code{"cox"}, uses the \code{coxph.fit} function with \code{method="efron"} argument set.
#' @param xrange range for evaluating \code{x}, default is \eqn{f} and \eqn{1 - f} quantiles of \code{x},
#' where \eqn{f = \frac{10}{\max{(n, 200)}}}{f = 10/max(\code{n}, 200)} and \eqn{n} the number of observations
#' @param event event/censoring indicator if \code{model="cox"}. If \code{event} is present, \code{model} is assumed to be \code{"cox"}
#' @param nk number of knots
#' @param knots knot locations, default based on quantiles of \code{x} (by \code{\link[Hmisc]{rcspline.eval}})
#' @param show \code{"xbeta"} or \code{"prob"} - what is plotted on \verb{y}-axis
#' @param adj optional matrix of adjustment variables
#' @param xlab \verb{x}-axis label, default is the \dQuote{label} attribute of \code{x}
#' @param ylab \verb{y}-axis label, default is the \dQuote{label} attribute of \code{y}
#' @param ylim \verb{y}-axis limits for logit or log hazard
#' @param plim \verb{y}-axis limits for probability scale
#' @param plotcl plot confidence limits
#' @param showknots show knot locations with arrows
#' @param add add this plot to an already existing plot
#' @param plot logical to indicate whether a plot has to be made. \code{FALSE} suppresses the plot.
#' @param subset subset of observations to process, e.g. \code{sex == "male"}
#' @param lty line type for plotting estimated spline function
#' @param noprint suppress printing regression coefficients and standard errors
#' @param m for \code{model="logistic"}, plot grouped estimates with triangles. Each group contains \code{m} ordered observations on \code{x}.
#' @param smooth plot nonparametric estimate if \code{model="logistic"} and \code{adj} is not specified
#' @param bass smoothing parameter (see \code{supsmu})
#' @param main main title, default is \code{"Estimated Spline Transformation"}
#' @param statloc location of summary statistics. Default positioning by clicking left mouse button where upper left corner of statistics should appear.
#'  Alternative is \code{"ll"} to place below the graph on the lower left, or the actual \code{x} and \code{y} coordinates. Use \code{"none"} to suppress statistics.
#'
#' @return list with components (\samp{knots}, \samp{x}, \samp{xbeta}, \samp{lower}, \samp{upper}) which are respectively the knot locations, design matrix,
#' linear predictor, and lower and upper confidence limits
#' @seealso   \code{\link[rms]{lrm}}, \code{\link[rms]{cph}}, \code{\link[Hmisc]{rcspline.eval}}, \code{\link[graphics]{plot}}, \code{\link[stats]{supsmu}},
#' \code{\link[survival:survival-internal]{coxph.fit}}, \code{\link[rms]{lrm.fit}}
rcspline_plot <- function(x, y, model=c("logistic","cox","ols"), xrange,
                           event, nk=5, knots=NULL, show=c("xbeta", "prob"),
                           adj=NULL, xlab, ylab, ylim, plim=c(0,1),
                           plotcl=TRUE, showknots=TRUE, add=FALSE, plot = TRUE, subset,
                           lty=1, noprint=FALSE, m, smooth=FALSE, bass=1,
                           main="auto", statloc)
{
  model <- match.arg(model)
  show  <- match.arg(show)
  if(plot) {
    oldpar = par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }
  
  if(! missing(event))
    model<-"cox"
  
  if(model == "cox" & missing(event))
    stop('event must be given for model="cox"')
  
  if(show == "prob" & ! missing(adj))
    stop('show="prob" cannot be used with adj')
  
  if(show == "prob" & model != "logistic")
    stop('show="prob" can only be used with model="logistic"')
  
  if(length(x) != length(y))
    stop('x and y must have the same length')
  
  if(! missing(event) && length(event) != length(y))
    stop('y and event must have the same length')
  
  if(! missing(adj)) {
    if(! is.matrix(adj)) adj <- as.matrix(adj)
    if(dim(adj)[1] != length(x))
      stop('x and adj must have the same length')
  }
  
  if(missing(xlab))
    xlab <- Hmisc::label(x)
  
  if(missing(ylab))
    ylab <- Hmisc::label(y)
  
  isna <- is.na(x) | is.na(y)
  if(! missing(event))
    isna <- isna | is.na(event)
  
  nadj <- 0
  if(! missing(adj)) {
    nadj <- ncol(adj)
    isna <- isna | apply(is.na(adj), 1, sum) > 0
  }
  
  if(! missing(subset))
    isna <- isna | (! subset)
  
  x <- x[! isna]
  y <- y[! isna]
  if(! missing(event))
    event <- event[! isna]
  
  if(! missing(adj))
    adj <- adj[! isna, ]
  
  n <- length(x)
  if(n<6)
    stop('fewer than 6 non-missing observations')
  
  if(missing(xrange)) {
    frac<-10./max(n, 200)
    xrange<-quantile(x, c(frac, 1.-frac))
  }
  
  if(missing(knots))
    xx <- Hmisc::rcspline.eval(x, nk=nk)
  else xx <- Hmisc::rcspline.eval(x, knots)
  
  knots <- attr(xx, "knots")
  nk <- length(knots)
  
  df1 <- nk-2
  if(model == "logistic") {
    b <- rms::lrm.fit(cbind(x, xx, adj),  y)
    beta <- b$coef
    cov <- b$var
    model.lr <- b$stats["Model L.R."]
    offset <- 1 	#to skip over intercept parameter
    ylabl <-
      if(show == "prob")
        "Probability"
    else "log Odds"
    
    sampled <- paste("Logistic Regression Model,  n=", n," d=", sum(y), sep="")
  }
  
  if(model == "cox") {
    if(! existsFunction('coxph.fit'))
      coxph.fit <- getFromNamespace('coxph.fit', 'survival')
    ##11mar04
    
    ## added coxph.control around iter.max, eps  11mar04
    lllin <- coxph.fit(cbind(x, adj), cbind(y, event), strata=NULL,
                       offset=NULL, init=NULL,
                       control=coxph.control(iter.max=10, eps=.0001),
                       method="efron", rownames=NULL)$loglik[2]
    b <- coxph.fit(cbind(x, xx, adj), cbind(y, event), strata=NULL,
                   offset=NULL, init=NULL,
                   control=coxph.control(iter.max=10, eps=.0001),
                   method="efron", rownames=NULL)
    beta <- b$coef
    if(! noprint) {
      print(beta);
      print(b$loglik)
    }
    
    beta <- b$coef
    cov <- b$var
    model.lr<-2*(b$loglik[2]-b$loglik[1])
    offset <- 0
    ylabl <- "log Relative Hazard"
    sampled <- paste("Cox Regression Model, n=",n," events=",sum(event),
                     sep="")
  }
  
  if(model == "logistic"|model == "cox") {
    model.df <- nk - 1 + nadj
    model.aic <- model.lr-2.*model.df
    v <- solve(cov[(1 + offset) : (nk + offset - 1), (1 + offset) : (nk + offset - 1)])
    assoc.chi <- beta[(1 + offset) : (nk + offset - 1)] %*% v %*%
      beta[(1 + offset) : (nk + offset - 1)]
    assoc.df <- nk - 1   #attr(v,"rank")
    assoc.p <- 1.-pchisq(assoc.chi, nk - 1)
    v <- solve(cov[(2 + offset) : (nk + offset - 1), (2 + offset) : (nk + offset - 1)])
    linear.chi <- beta[(2 + offset) : (nk + offset - 1)] %*% v %*%
      beta[(2 + offset) : (nk + offset - 1)]
    linear.df <- nk - 2   #attr(v,"rank")
    linear.p <- 1. - pchisq(linear.chi, linear.df)
    if(nadj > 0) {
      ntot <- offset + nk - 1 + nadj
      v <- solve(cov[(nk + offset) : ntot, (nk + offset) : ntot])
      adj.chi <- beta[(nk + offset) : ntot] %*% v %*%
        beta[(nk + offset) : ntot]
      adj.df <- ncol(v)   #attr(v,"rank")
      adj.p <- 1. - pchisq(adj.chi, adj.df)
    } else {
      adj.chi <- 0
      adj.p <- 0
    }
  }
  
  ## Evaluate xbeta for expanded x at desired range
  xe <- seq(xrange[1], xrange[2], length=600)
  if(model == "cox")
    xx <- Hmisc::rcspline.eval(xe, knots, inclx=TRUE)
  else
    xx<- cbind(rep(1, length(xe)), Hmisc::rcspline.eval(xe, knots, inclx=TRUE))
  
  xbeta <- xx %*% beta[1 : (nk - 1 + offset)]
  var <- drop(((xx %*% cov[1 : (nk - 1 + offset), 1 : (nk - 1 + offset)])*xx) %*%
                rep(1, ncol(xx)))
  lower <- xbeta - 1.96*sqrt(var)
  upper <- xbeta + 1.96*sqrt(var)
  if(show == "prob") {
    xbeta <- 1./(1. + exp(-xbeta))
    lower <- 1./(1. + exp(-lower))
    upper <- 1./(1. + exp(-upper))
  }
  
  xlim <- range(pretty(xe))
  if(missing(ylim))
    ylim <- range(pretty(c(xbeta, if(plotcl) lower, if(plotcl) upper)))
  
  if(main == "auto") {
    if(show == "xbeta")
      main <- "Estimated Spline Transformation"
    else main <- "Spline Estimate of Prob{Y=1}"
  }
  
  if(! interactive() & missing(statloc))
    statloc<-"ll"
  
  if(plot) {
    if(! add) {
      oldmar<-par("mar")
      if(! missing(statloc) && statloc[1] == "ll")
        oldmar[1]<- 11
      
      plot(xe, xbeta, type="n", main=main, xlab=xlab, ylab=ylabl,
           xlim=xlim, ylim=ylim)
      lines(xe, xbeta, lty=lty)
      ltext<-function(z, line, label, cex=.8, adj=0)
      {
        zz<-z
        zz$y<-z$y-(line - 1)*1.2*cex*par("csi")*(par("usr")[4]-par("usr")[3])/
          (par("fin")[2])   #was 1.85
        text(zz, label, cex=cex, adj=adj)
      }
      
      sl<-0
      if(missing(statloc)) {
        message("Click left mouse button at upper left corner for statistics\n")
        z<-locator(1)
        statloc<-"l"
      } else if(statloc[1] != "none") {
        if(statloc[1] == "ll") {
          z<-list(x=par("usr")[1], y=par("usr")[3])
          sl<-3
        } else z<-list(x=statloc[1], y=statloc[2])
      }
      
      if(statloc[1] != "none" & (model == "logistic" | model == "cox"))	{
        rnd <- function(x, r=2) as.single(round(x, r))
        
        ltext(z, 1 + sl, sampled)
        ltext(z, 2 + sl, "    Statistic        X2  df")
        chistats<-format(as.single(round(c(model.lr, model.aic,
                                           assoc.chi, linear.chi, adj.chi), 2)))
        pvals<-format(as.single(round(c(assoc.p, linear.p, adj.p), 4)))
        ltext(z, 3 + sl, paste("Model        L.R. ", chistats[1], model.df,
                               " AIC=", chistats[2]))
        ltext(z, 4 + sl, paste("Association  Wald ", chistats[3], assoc.df,
                               " p= ", pvals[1]))
        ltext(z, 5 + sl, paste("Linearity    Wald ", chistats[4], linear.df,
                               " p= ", pvals[2]))
        if(nadj > 0)ltext(z, 6 + sl, paste("Adjustment   Wald " , chistats[5],
                                           adj.df, " p= ", pvals[3]))}
    } else lines(xe, xbeta, lty=lty)
    
    if(plotcl) {
      #prn(cbind(xe, lower, upper))
      lines(xe, lower, lty=2)
      lines(xe, upper, lty=2)
    }
    
    if(showknots) {
      bot.arrow <- par("usr")[3]
      top.arrow <- bot.arrow + .05 * (par("usr")[4]-par("usr")[3])
      for(i in 1 : nk)
        arrows(knots[i], top.arrow, knots[i], bot.arrow, length=.1)
    }
    
    if(model == "logistic" & nadj == 0) {
      if(smooth) {
        z<-supsmu(x, y, bass=bass)
        if(show == "xbeta") z$y <- logb(z$y/(1.-z$y))
        points(z, cex=.4)
      }
      
      if(! missing(m)) {
        z<-groupn(x, y, m=m)
        if(show == "xbeta") z$y <- logb(z$y/(1.-z$y))
        
        points(z, pch=2, mkh=.05)}
    }
  }
  
  invisible(list(
    knots = knots,
    x = xe,
    xbeta = xbeta,
    lower = lower,
    upper = upper
  ))
}



# riskRegression::predictCox does not know how to handle frailty, manual calculation of risk for CSC is required
# this is adapted from https://github.com/survival-lumc/ValidationCompRisks/blob/main/sharing_CSC_model.R

predictLM_CSC <- function(fit_model, model_info, # List object (see above)
                          newdata, # Data.frame for which to make predictions
                          horizon, # Prediction time horizon (numeric)
                          primary_event) { # Primary event (numeric)
  
  # n_causes <- unique(vapply(model_info, length, FUN.VALUE = integer(1L)))
  n_causes <- length(model_info)
  # -- Absolute risk prediction
  
  causes_ind <- seq_len(n_causes)
  
  # Calculate linear predictors for all causes in new dataset
  linpreds <- lapply(causes_ind, function(cause) {
    mod_matrix <- predict(fit_model$models[[cause]], newdata, type="lp", reference = "zero")
  })
  
  
  # Compute absolute risks for each individual
  preds <- vapply(seq_len(nrow(newdata)), function(id) {
    
    # Calculate individual-specific cause-specific hazards
    time_points <- model_info$baseline_hazards[[primary_event]][["time"]]
    hazards <- vapply(causes_ind, function(cause) {
      cumhaz <- model_info$baseline_hazards[[cause]][["hazard"]] * exp(linpreds[[cause]][[id]])
      diff(c(0, cumhaz))
    }, FUN.VALUE = numeric(length(time_points)))
    
    hazards <- cbind(hazards, time_points)
    lmsi <- newdata$LM[id]
    hazards<- hazards[hazards[,"time_points"] <= lmsi+horizon & hazards[,"time_points"] > lmsi,]
    # Calculate event-free survival
    surv <- cumprod(1 - rowSums(hazards[,1:n_causes]))
    # Calculate cumulative incidence
    cuminc <- cumsum(hazards[, primary_event] * c(1, surv[-length(surv)]))
    cuminc_horizon <- tail(cuminc, n=1)
    
    return(cuminc_horizon)
  }, FUN.VALUE = numeric(1L))
}