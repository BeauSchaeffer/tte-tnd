##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Equi-Confounding Analysis Pooled

# setwd("~/Desktop/Research/Kaiser/KP_analysis")


# Packages ----------------------------------------------------------------

library(tidyverse)
library(data.table)
library(speedglm)
library(splines)
# library(geepack) # called without loading

# Data --------------------------------------------------------------------


data_Y3 <- readr::read_rds("cleaned_data/data_Y3.rds")


# Downsample --------------------------------------------------------------


subclass_ids <- data_Y3 |> dplyr::select(subclass) |> unique()
set.seed(345)
subclass_ids_subset <- dplyr::slice_sample(subclass_ids, n=10000)
dat_downsamp <- data_Y3 |> dplyr::filter(subclass %in% subclass_ids_subset$subclass) |> droplevels()
rm(subclass_ids, subclass_ids_subset)

# IPTW --------------------------------------------------------------------

  ### see standard_analysis_pooled if needed

# ITT long format expansion -----------------------------------------------


  ### calc number of rows needed for each individual
time_unit <- 1
dat_downsamp$max_units <- ceiling(dat_downsamp$Y3_itt_t_trunc/time_unit)+1
  # analysis_data_Y_TWO$max_units <- pmax(1, ceiling(analysis_data_Y_TWO$Y2_itt_t / time_unit))
  ### above line adjusted from causal lab code to ensure that
  ### at least 1 row created for each individual
dat_downsamp.long.itt <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$max_units),]

  ### Now, let's create a variable that's represents the start and end time
  ### corresponding to each row of observation

dat_downsamp.long.itt$time_start <- ave(dat_downsamp.long.itt$fake_mrn, dat_downsamp.long.itt$fake_mrn, FUN=seq_along)
dat_downsamp.long.itt$time_start <- (dat_downsamp.long.itt$time_start-1)*time_unit
dat_downsamp.long.itt$time_end <- dat_downsamp.long.itt$time_start+time_unit

  ### Last, we have to modify the Y and C variables so that they are only
  ### equal to 1 if the event/censoring happened in that time interval

dat_downsamp.long.itt$Y_pos <- ifelse(
  dat_downsamp.long.itt$Y3_itt_trunc == 2 &
    dat_downsamp.long.itt$Y3_itt_t_trunc == dat_downsamp.long.itt$time_start,
  1, 0
)

dat_downsamp.long.itt$Y_neg <- ifelse(
  dat_downsamp.long.itt$Y3_itt_trunc == 1 &
    dat_downsamp.long.itt$Y3_itt_t_trunc == dat_downsamp.long.itt$time_start,
  1, 0
)

dat_downsamp.long.itt$C <- ifelse(
  dat_downsamp.long.itt$Y3_itt_trunc == 0 &
    dat_downsamp.long.itt$Y3_itt_t_trunc == dat_downsamp.long.itt$time_start,
  1, 0
)


dat_downsamp.long.itt$Y_pos <- ifelse(dat_downsamp.long.itt$C==1, NA, dat_downsamp.long.itt$Y_pos)
dat_downsamp.long.itt$Y_neg <- ifelse(dat_downsamp.long.itt$C==1, NA, dat_downsamp.long.itt$Y_neg)


# ITT IPCW ----------------------------------------------------------------


### No IPCW applied to ITT analysis


# ITT Pooled Logistic -----------------------------------------------------


eqc_pooled_itt_fit1 <- glm(Y_neg ~ ns(time_end, knots = c(10,20,30))*treatment +
                                  # demographic
                                  sex_admin + age_years + bmi + race + charlson_cat_fac +
                                  # other
                                  ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                  # NEC
                                  flu_vax,
                                data=dat_downsamp.long.itt,
                                family=binomial())
summary(eqc_pooled_itt_fit1)

eqc_pooled_itt_fit2 <- glm(Y_pos ~ ns(time_end, knots = c(10,20,30))*treatment +
                                  # demographic
                                  sex_admin + age_years + bmi + race + charlson_cat_fac +
                                  # other
                                  ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                                  # NEC
                                  flu_vax,
                                data=dat_downsamp.long.itt,
                                family=binomial())
summary(eqc_pooled_itt_fit2)


# ITT Survival and Risk ---------------------------------------------------


dat_downsamp$gmaxt <- 53

### G formula data setup A=0
eqc_itt_A0.long <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$gmaxt),]
eqc_itt_A0.long$time_start <- ave(eqc_itt_A0.long$fake_mrn, eqc_itt_A0.long$fake_mrn, FUN=seq_along)
eqc_itt_A0.long$time_start <- (eqc_itt_A0.long$time_start-1)*time_unit
eqc_itt_A0.long$time_end <- eqc_itt_A0.long$time_start+time_unit
eqc_itt_A0.long$treatment <- 0

### G formula data setup A=1
eqc_itt_A1.long <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$gmaxt),]
eqc_itt_A1.long$time_start <- ave(eqc_itt_A1.long$fake_mrn, eqc_itt_A1.long$fake_mrn, FUN=seq_along)
eqc_itt_A1.long$time_start <- (eqc_itt_A1.long$time_start-1)*time_unit
eqc_itt_A1.long$time_end <- eqc_itt_A1.long$time_start+time_unit
eqc_itt_A1.long$treatment <- 1

  # IPTW data setup
  # eqc_itt_A0 <- data.frame(treatment=0, time_end=unique(dat_downsamp.long.itt$time_end))
  # eqc_itt_A1 <- data.frame(treatment=1, time_end=unique(dat_downsamp.long.itt$time_end))
  # eqc_itt_A0
  # eqc_itt_A1

### Calculate predicted hazards:
eqc_itt_A0.long$hazard_pos <- predict(eqc_pooled_itt_fit2, newdata=eqc_itt_A0.long, type="response")
eqc_itt_A1.long$hazard_pos <- predict(eqc_pooled_itt_fit2, newdata=eqc_itt_A1.long, type="response")
eqc_itt_A0.long$hazard_neg <- predict(eqc_pooled_itt_fit1, newdata=eqc_itt_A0.long, type="response")
eqc_itt_A1.long$hazard_neg <- predict(eqc_pooled_itt_fit1, newdata=eqc_itt_A1.long, type="response")
### Corrected hazards under no treatment
eqc_itt_A0.long$hazard_pos_c <- eqc_itt_A0.long$hazard_pos * (eqc_itt_A1.long$hazard_neg / eqc_itt_A0.long$hazard_neg)

### Calculate (1 - hazard)
eqc_itt_A0.long$pnoevent_pos <- 1 - eqc_itt_A0.long$hazard_pos
eqc_itt_A1.long$pnoevent_pos <- 1 - eqc_itt_A1.long$hazard_pos
eqc_itt_A0.long$pnoevent_neg <- 1 - eqc_itt_A0.long$hazard_neg
eqc_itt_A1.long$pnoevent_neg <- 1 - eqc_itt_A1.long$hazard_neg
### Corrected (1 - hazard) under no treatment
eqc_itt_A0.long$pnoevent_pos_c <- 1 - eqc_itt_A0.long$hazard_pos_c

### Sort the data by ID, time
eqc_itt_A0.long <- eqc_itt_A0.long[order(eqc_itt_A0.long$fake_mrn, eqc_itt_A0.long$time_end),] 
eqc_itt_A1.long <- eqc_itt_A1.long[order(eqc_itt_A1.long$fake_mrn, eqc_itt_A1.long$time_end),] 

### Calculate the cumulative survival by taking the cumulative product
### of (1 - hazard)

  # IPTW approach

  # eqc_itt_A0$survival_pos <- cumprod(eqc_itt_A0$pnoevent_neg * lag(eqc_itt_A0$pnoevent_pos, n = 1, default = 1))
  # eqc_itt_A1$survival_pos <- cumprod(eqc_itt_A1$pnoevent_neg * lag(eqc_itt_A1$pnoevent_pos, n = 1, default = 1))
  # ### Corrected cumulative survival under no treatment 
  # eqc_itt_A0$survival_pos_c <- cumprod(eqc_itt_A0$pnoevent_neg * lag(eqc_itt_A0$pnoevent_pos_c, n = 1, default = 1))

  # G FORM approach

eqc_itt_A0.long <- eqc_itt_A0.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1),
         pnoevent_pos_c_lag = lag(pnoevent_pos_c, n=1, default=1)) |> 
  ungroup()

eqc_itt_A1.long <- eqc_itt_A1.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1)) |> 
  ungroup()

    # product at each time (P(no event neg) * lag P(no event pos))

eqc_itt_A0.long$surv_prod_lag <- eqc_itt_A0.long$pnoevent_neg * eqc_itt_A0.long$pnoevent_pos_lag
eqc_itt_A1.long$surv_prod_lag <- eqc_itt_A1.long$pnoevent_neg * eqc_itt_A1.long$pnoevent_pos_lag
eqc_itt_A0.long$surv_prod_c_lag <- eqc_itt_A0.long$pnoevent_neg * eqc_itt_A0.long$pnoevent_pos_c_lag

    # cumulative product within individual

eqc_itt_A0.long$survival_pos <- ave(eqc_itt_A0.long$surv_prod_lag, eqc_itt_A0.long$fake_mrn, FUN=cumprod)
eqc_itt_A1.long$survival_pos <- ave(eqc_itt_A1.long$surv_prod_lag, eqc_itt_A1.long$fake_mrn, FUN=cumprod)
eqc_itt_A0.long$survival_pos_c <- ave(eqc_itt_A0.long$surv_prod_c_lag, eqc_itt_A0.long$fake_mrn, FUN=cumprod)

### Calculate risk using CIF estimator
  
  # IPTW approach

  # eqc_itt_A0$risk_pos <- cumsum(eqc_itt_A0$hazard_pos * eqc_itt_A0$survival_pos)
  # eqc_itt_A1$risk_pos <- cumsum(eqc_itt_A1$hazard_pos * eqc_itt_A1$survival_pos)
  # eqc_itt_A0$risk_pos_c <- cumsum(eqc_itt_A0$hazard_pos_c * eqc_itt_A0$survival_pos_c)

  # G FORM approach

    # product at each time (haz pos * surv pos)

eqc_itt_A0.long$risk_prod_pos <- eqc_itt_A0.long$hazard_pos * eqc_itt_A0.long$survival_pos
eqc_itt_A1.long$risk_prod_pos <- eqc_itt_A1.long$hazard_pos * eqc_itt_A1.long$survival_pos
eqc_itt_A0.long$risk_prod_pos_c <- eqc_itt_A0.long$hazard_pos_c * eqc_itt_A0.long$survival_pos_c

    # cumulative sum within individual

eqc_itt_A0.long$risk_pos <- ave(eqc_itt_A0.long$risk_prod_pos, eqc_itt_A0.long$fake_mrn, FUN=cumsum)
eqc_itt_A1.long$risk_pos <- ave(eqc_itt_A1.long$risk_prod_pos, eqc_itt_A1.long$fake_mrn, FUN=cumsum)
eqc_itt_A0.long$risk_pos_c <- ave(eqc_itt_A0.long$risk_prod_pos_c, eqc_itt_A0.long$fake_mrn, FUN=cumsum)

# Calculate the average risk at each time point

eqc_itt_A0.long.res <- aggregate(risk_pos ~ time_end, data=eqc_itt_A0.long, FUN=mean)
eqc_itt_A1.long.res <- aggregate(risk_pos ~ time_end, data=eqc_itt_A1.long, FUN=mean)
eqc_itt_A0.long.res.c <- aggregate(risk_pos_c ~ time_end, data=eqc_itt_A0.long, FUN=mean)

### plot the risk curves

# png("results/eqc_Y3_risks_itt_p.png", width = 2400, height = 1800, res=300)

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, eqc_itt_A0.long.res$time_end, eqc_itt_A1.long.res$time_end)),
     ylim = range(c(0, 0.10)),
     xlab="Weeks",
     ylab="Risk",
     main="Risk Curves",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main=1.4
)
mtext("Equi-confounding (ITT)", side = 3, line = 0.5, font = 3, cex=1.2)
# lines(c(0, eqc_itt_A0.long.res$time_end), c(0, eqc_itt_A0.long.res$risk_pos), col='#006663', lty=2)
grid()
lines(c(0, eqc_itt_A1.long.res$time_end), c(0, eqc_itt_A1.long.res$risk_pos), col='#FF6B1A', lty=1, lwd=4)
lines(c(0, eqc_itt_A0.long.res.c$time_end), c(0, eqc_itt_A0.long.res.c$risk_pos_c), col='#006663', lty=1, lwd=4)
legend("topleft",
       legend = c("No Booster", "Booster"),
       col = c("#006663", "#FF6B1A"),
       lty = 1, lwd = 4, bty = "n", cex=1.2)

# legend("topright",
#        legend = c("Corrected", "Original"),
#        lty = c(1, 2), lwd = 2, bty = "n")

# dev.off() # 2025-12-10

### Plot RR over time

setDT(eqc_itt_A0.long.res)
setDT(eqc_itt_A1.long.res)
setDT(eqc_itt_A0.long.res.c)

eqc_RR_itt <- data.table(
  Week = 1:53,
  RR = sapply(1:53, function(wk) {
    num_pos <- eqc_itt_A1.long.res[time_end == wk, risk_pos]
    denom_pos <- eqc_itt_A0.long.res[time_end == wk, risk_pos]
    return(
      as.numeric(num_pos) / as.numeric(denom_pos)
    )
  }),
  RR_c = sapply(1:53, function(wk) {
    num_pos <- eqc_itt_A1.long.res[time_end == wk, risk_pos]
    denom_pos <- eqc_itt_A0.long.res.c[time_end == wk, risk_pos_c]
    return(
      as.numeric(num_pos) / as.numeric(denom_pos)
    )
  })
)

# write_rds(eqc_RR_itt, file = "results/eqc_RR_itt.rds") # 2025-12-10

# png("results/eqc_Y3_RR_itt_p.png", width = 2400, height = 1800, res=300)

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, eqc_RR_itt$Week)),
     ylim = range(c(0.0, 1.2)),
     xlab = "Weeks",
     ylab = "Risk Ratio (RR)",
     main = "Risk Ratio Over Time",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main=1.4)
grid()
lines(eqc_RR_itt$Week, eqc_RR_itt$RR_c, col='black', lty=1, lwd=4)
# lines(eqc_RR_itt$Week, eqc_RR_itt$RR, col='black', lty=2, lwd=2)
# legend("topright",
#        legend = c("Corrected", "Original"),
#        lty = c(1, 2), lwd = 2, bty = "n")
mtext("Equi-confounding (ITT)", side = 3, line = 0.5, font = 3, cex=1.2)
abline(h = 1, col = "black", lty = 1, lwd = 0.5)

# dev.off() # 2025-12-10


# PP long format expansion ------------------------------------------------


### calc number of rows needed for each individual
time_unit <- 1
dat_downsamp$max_units_pp <- ceiling(dat_downsamp$Y3_pp_t_trunc/time_unit)+1
dat_downsamp.long.pp <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$max_units_pp),]

### start and end time corresponding to each row of observation

dat_downsamp.long.pp$time_start <- ave(dat_downsamp.long.pp$fake_mrn, dat_downsamp.long.pp$fake_mrn, FUN=seq_along)
dat_downsamp.long.pp$time_start <- (dat_downsamp.long.pp$time_start-1)*time_unit
dat_downsamp.long.pp$time_end <- dat_downsamp.long.pp$time_start+time_unit

### modify Y and C variables so that they are only equal to 1 if the 
### event/censoring happened in that time interval

dat_downsamp.long.pp$Y_pos <- ifelse(
  dat_downsamp.long.pp$Y3_pp_trunc == 2 &
    dat_downsamp.long.pp$Y3_pp_t_trunc == dat_downsamp.long.pp$time_start,
  1, 0
)

dat_downsamp.long.pp$Y_neg <- ifelse(
  dat_downsamp.long.pp$Y3_pp_trunc == 1 &
    dat_downsamp.long.pp$Y3_pp_t_trunc == dat_downsamp.long.pp$time_start,
  1, 0
)

dat_downsamp.long.pp$C <- ifelse(
  dat_downsamp.long.pp$Y3_pp_trunc == 0 &
    dat_downsamp.long.pp$Y3_pp_t_trunc == dat_downsamp.long.pp$time_start,
  1, 0
)


dat_downsamp.long.pp$Y_pos <- ifelse(dat_downsamp.long.pp$C==1, NA, dat_downsamp.long.pp$Y_pos)
dat_downsamp.long.pp$Y_neg <- ifelse(dat_downsamp.long.pp$C==1, NA, dat_downsamp.long.pp$Y_neg)


# PP IPCW -----------------------------------------------------------------


dat_downsamp.long.pp <- dat_downsamp.long.pp |> arrange(fake_mrn, time_end)

ipcw_denom <- glm(C ~ ns(time_end, knots = c(10,20,30))*treatment +
                         # demographic
                         sex_admin + age_years + bmi + race + charlson_cat_fac +
                         # other
                         ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                         # NEC
                         flu_vax,
                       data = dat_downsamp.long.pp, family = binomial(link = "logit")
                       )

dat_downsamp.long.pp$pd.cens <- 1-predict(ipcw_denom, type = "response", newdata = dat_downsamp.long.pp)

ipcw_num <- speedglm(C ~ treatment,
                     data = dat_downsamp.long.pp,
                     family = binomial(link = "logit")
                     )

dat_downsamp.long.pp$pn.cens <- 1-predict(ipcw_num, type = "response", newdata = dat_downsamp.long.pp)

dat_downsamp.long.pp <- dat_downsamp.long.pp |>
  group_by(fake_mrn) |>
  mutate(
    pdcuml.cens = cumprod(pd.cens),
    pncuml.cens = cumprod(pn.cens)) |>
  mutate(ipcw = pncuml.cens / pdcuml.cens) |>
  ungroup()

# dat_downsamp.long.pp$sw.cw <-  dat_downsamp.long.pp$sw * dat_downsamp.long.pp$ipcw


# PP Pooled Logistic ------------------------------------------------------


eqc_pooled_pp_fit1 <- glm(Y_neg ~ ns(time_end, knots = c(10,20,30))*treatment +
                             # demographic
                             sex_admin + age_years + bmi + race + charlson_cat_fac +
                             # other
                             ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                             # NEC
                             flu_vax,
                           data=dat_downsamp.long.pp,
                           family=binomial(), weights=ipcw)


eqc_pooled_pp_fit2 <- glm(Y_pos ~ ns(time_end, knots = c(10,20,30))*treatment +
                            # demographic
                            sex_admin + age_years + bmi + race + charlson_cat_fac +
                            # other
                            ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
                            # NEC
                            flu_vax,
                          data=dat_downsamp.long.pp,
                          family=binomial(), weights=ipcw)

# PP survival and risk ----------------------------------------------------

dat_downsamp$gmaxt <- 53

### G formula data setup A=0
eqc_pp_A0.long <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$gmaxt),]
eqc_pp_A0.long$time_start <- ave(eqc_pp_A0.long$fake_mrn, eqc_pp_A0.long$fake_mrn, FUN=seq_along)
eqc_pp_A0.long$time_start <- (eqc_pp_A0.long$time_start-1)*time_unit
eqc_pp_A0.long$time_end <- eqc_pp_A0.long$time_start+time_unit
eqc_pp_A0.long$treatment <- 0

### G formula data setup A=1
eqc_pp_A1.long <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$gmaxt),]
eqc_pp_A1.long$time_start <- ave(eqc_pp_A1.long$fake_mrn, eqc_pp_A1.long$fake_mrn, FUN=seq_along)
eqc_pp_A1.long$time_start <- (eqc_pp_A1.long$time_start-1)*time_unit
eqc_pp_A1.long$time_end <- eqc_pp_A1.long$time_start+time_unit
eqc_pp_A1.long$treatment <- 1

### Calculate predicted hazards:
eqc_pp_A0.long$hazard_pos <- predict(eqc_pooled_pp_fit2, newdata=eqc_pp_A0.long, type="response")
eqc_pp_A1.long$hazard_pos <- predict(eqc_pooled_pp_fit2, newdata=eqc_pp_A1.long, type="response")
eqc_pp_A0.long$hazard_neg <- predict(eqc_pooled_pp_fit1, newdata=eqc_pp_A0.long, type="response")
eqc_pp_A1.long$hazard_neg <- predict(eqc_pooled_pp_fit1, newdata=eqc_pp_A1.long, type="response")
  ### Corrected hazards under no treatment
eqc_pp_A0.long$hazard_pos_c <- eqc_pp_A0.long$hazard_pos * (eqc_pp_A1.long$hazard_neg / eqc_pp_A0.long$hazard_neg)

### Calculate (1 - hazard)
eqc_pp_A0.long$pnoevent_pos <- 1 - eqc_pp_A0.long$hazard_pos
eqc_pp_A1.long$pnoevent_pos <- 1 - eqc_pp_A1.long$hazard_pos
eqc_pp_A0.long$pnoevent_neg <- 1 - eqc_pp_A0.long$hazard_neg
eqc_pp_A1.long$pnoevent_neg <- 1 - eqc_pp_A1.long$hazard_neg
### Corrected (1 - hazard) under no treatment
eqc_pp_A0.long$pnoevent_pos_c <- 1 - eqc_pp_A0.long$hazard_pos_c

### Sort the data by ID, time
eqc_pp_A0.long <- eqc_pp_A0.long[order(eqc_pp_A0.long$fake_mrn, eqc_pp_A0.long$time_end),] 
eqc_pp_A1.long <- eqc_pp_A1.long[order(eqc_pp_A1.long$fake_mrn, eqc_pp_A1.long$time_end),]

### Calculate the cumulative survival by taking the cumulative product
### of (1 - hazard)

  # G FORM approach

eqc_pp_A0.long <- eqc_pp_A0.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1),
         pnoevent_pos_c_lag = lag(pnoevent_pos_c, n=1, default=1)) |> 
  ungroup()

eqc_pp_A1.long <- eqc_pp_A1.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1)) |> 
  ungroup()

  # product at each time (P(no event neg) * lag P(no event pos))

eqc_pp_A0.long$surv_prod_lag <- eqc_pp_A0.long$pnoevent_neg * eqc_pp_A0.long$pnoevent_pos_lag
eqc_pp_A1.long$surv_prod_lag <- eqc_pp_A1.long$pnoevent_neg * eqc_pp_A1.long$pnoevent_pos_lag
eqc_pp_A0.long$surv_prod_c_lag <- eqc_pp_A0.long$pnoevent_neg * eqc_pp_A0.long$pnoevent_pos_c_lag

  # cumulative product within individual

eqc_pp_A0.long$survival_pos <- ave(eqc_pp_A0.long$surv_prod_lag, eqc_pp_A0.long$fake_mrn, FUN=cumprod)
eqc_pp_A1.long$survival_pos <- ave(eqc_pp_A1.long$surv_prod_lag, eqc_pp_A1.long$fake_mrn, FUN=cumprod)
eqc_pp_A0.long$survival_pos_c <- ave(eqc_pp_A0.long$surv_prod_c_lag, eqc_pp_A0.long$fake_mrn, FUN=cumprod)

### Calculate risk using CIF estimator

  # product at each time (haz pos * surv pos)

eqc_pp_A0.long$risk_prod_pos <- eqc_pp_A0.long$hazard_pos * eqc_pp_A0.long$survival_pos
eqc_pp_A1.long$risk_prod_pos <- eqc_pp_A1.long$hazard_pos * eqc_pp_A1.long$survival_pos
eqc_pp_A0.long$risk_prod_pos_c <- eqc_pp_A0.long$hazard_pos_c * eqc_pp_A0.long$survival_pos_c

  # cumulative sum within individual

eqc_pp_A0.long$risk_pos <- ave(eqc_pp_A0.long$risk_prod_pos, eqc_pp_A0.long$fake_mrn, FUN=cumsum)
eqc_pp_A1.long$risk_pos <- ave(eqc_pp_A1.long$risk_prod_pos, eqc_pp_A1.long$fake_mrn, FUN=cumsum)
eqc_pp_A0.long$risk_pos_c <- ave(eqc_pp_A0.long$risk_prod_pos_c, eqc_pp_A0.long$fake_mrn, FUN=cumsum)

  # Calculate the average risk at each time point

eqc_pp_A0.long.res <- aggregate(risk_pos ~ time_end, data=eqc_pp_A0.long, FUN=mean)
eqc_pp_A1.long.res <- aggregate(risk_pos ~ time_end, data=eqc_pp_A1.long, FUN=mean)
eqc_pp_A0.long.res.c <- aggregate(risk_pos_c ~ time_end, data=eqc_pp_A0.long, FUN=mean)


### plot the risk curves

# png("results/eqc_Y3_risks_pp_p.png", width = 2400, height = 1800, res=300)

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, eqc_pp_A0.long.res$time_end, eqc_pp_A1.long.res$time_end)),
     ylim = range(c(0, 0.1)),
     xlab="Weeks",
     ylab="Risk",
     main="Risk Curves",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main=1.4
)
mtext("Equi-confounding (PP)", side = 3, line = 0.5, font = 3, cex=1.2)
grid()
# lines(c(0, eqc_pp_A0.long.res$time_end), c(0, eqc_pp_A0.long.res$risk_pos), col='#006663', lty=2, lwd=4)
lines(c(0, eqc_pp_A1.long.res$time_end), c(0, eqc_pp_A1.long.res$risk_pos), col='#FF6B1A', lty=1, lwd=4)
lines(c(0, eqc_pp_A0.long.res.c$time_end), c(0, eqc_pp_A0.long.res.c$risk_pos_c), col='#006663', lty=1, lwd=4)
legend("topleft",
       legend = c("No Booster", "Booster"),
       col = c("#006663", "#FF6B1A"),
       lty = 1, lwd = 4, bty = "n", cex=1.2)

# legend("topright",
#        legend = c("Corrected", "Original"),
#        lty = c(1, 2), lwd = 4, bty = "n")

# dev.off() # 2025-12-10


### Plot RR over time

setDT(eqc_pp_A0.long.res)
setDT(eqc_pp_A1.long.res)
setDT(eqc_pp_A0.long.res.c)

eqc_RR_pp <- data.table(
  Week = 1:53,
  RR = sapply(1:53, function(wk) {
    num_pos <- eqc_pp_A1.long.res[time_end == wk, risk_pos]
    denom_pos <- eqc_pp_A0.long.res[time_end == wk, risk_pos]
    return(
      as.numeric(num_pos) / as.numeric(denom_pos)
    )
  }),
  RR_c = sapply(1:53, function(wk) {
    num_pos <- eqc_pp_A1.long.res[time_end == wk, risk_pos]
    denom_pos <- eqc_pp_A0.long.res.c[time_end == wk, risk_pos_c]
    return(
      as.numeric(num_pos) / as.numeric(denom_pos)
    )
  })
)

# write_rds(eqc_RR_pp, file = "results/eqc_RR_pp.rds") # 2025-12-10

# png("results/eqc_Y3_RR_pp_p.png", width = 2400, height = 1800, res=300)

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, eqc_RR_pp$Week)),
     ylim = range(c(0, 1.2)),
     xlab = "Weeks",
     ylab = "Risk Ratio (RR)",
     main = "Risk Ratio Over Time",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main=1.4)
grid()
lines(eqc_RR_pp$Week, eqc_RR_pp$RR_c, col='black', lty=1, lwd=4)
# lines(eqc_RR_pp$Week, eqc_RR_pp$RR, col='black', lty=2, lwd=4)
# legend("topright",
#        legend = c("Corrected", "Original"),
#        lty = c(1, 2), lwd = 2, bty = "n")
mtext("Equi-confounding (PP)", side = 3, line = 0.5, font = 3, cex=1.2)
abline(h = 1, col = "black", lty = 1, lwd = 0.5)

# dev.off() # 2025-12-10




