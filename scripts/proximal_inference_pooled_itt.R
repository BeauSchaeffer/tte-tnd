##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Proximal Inference Analysis Pooled

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

time_unit <- 1
dat_downsamp$max_units <- ceiling(dat_downsamp$Y3_itt_t_trunc/time_unit)+1
dat_downsamp.long.itt <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$max_units),]

dat_downsamp.long.itt$time_start <- ave(dat_downsamp.long.itt$fake_mrn, dat_downsamp.long.itt$fake_mrn, FUN=seq_along)
dat_downsamp.long.itt$time_start <- (dat_downsamp.long.itt$time_start-1)*time_unit
dat_downsamp.long.itt$time_end <- dat_downsamp.long.itt$time_start+time_unit

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

  ### time interacting with all variables, note ()
  ### takes a while to fit
prox_pooled_itt_s1 <- glm(Y_neg ~ ns(time_end, knots = c(10,20,30))*(treatment +
                             # demographic
                             sex_admin + age_years + bmi + race + charlson_cat_fac +
                             # other
                             ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                             # NEC
                             flu_vax),
                           data=dat_downsamp.long.itt,
                           family=binomial())

dat_downsamp.long.itt$p_itt <- predict(prox_pooled_itt_s1, newdata = dat_downsamp.long.itt)

prox_pooled_itt_s2 <- glm(Y_pos ~ ns(time_end, knots = c(10,20,30))*treatment +
                            # demographic
                            sex_admin + age_years + bmi + race + charlson_cat_fac +
                            # other
                            ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                            # predictions from stage 1
                            p_itt,
                            # no NEC
                          data=dat_downsamp.long.itt,
                          family=binomial())

prox_pooled_itt_obs <- glm(Y_pos ~ ns(time_end, knots = c(10,20,30))*treatment +
                             # demographic
                             sex_admin + age_years + bmi + race + charlson_cat_fac +
                             # other
                             ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                             # NEC
                             flu_vax,
                           data=dat_downsamp.long.itt,
                           family=binomial())


df_ref_A1 <- data.frame(time_end=seq(1,53,1),
                        treatment=1,
                        sex_admin=factor("F"),
                        age_years=0,
                        bmi=0,
                        race=factor("White"),
                        charlson_cat_fac=factor("0"),
                        ndi=0,
                        prior_inf=0,
                        tests_count=0,
                        service_region=factor("Central valley"),
                        last_vax_infect_weeks=0,
                        p_itt=0)

df_ref_A0 <- data.frame(time_end=seq(1,53,1),
           treatment=0,
           sex_admin=factor("F"),
           age_years=0,
           bmi=0,
           race=factor("White"),
           charlson_cat_fac=factor("0"),
           ndi=0,
           prior_inf=0,
           tests_count=0,
           service_region=factor("Central valley"),
           last_vax_infect_weeks=0,
           p_itt=0)

haz_ref_A1 <- predict(prox_pooled_itt_s2, newdata=df_ref_A1, type = "link")
haz_ref_A0 <- predict(prox_pooled_itt_s2, newdata=df_ref_A0, type = "link")

exp(haz_ref_A1)/exp(haz_ref_A0)

  ### compare no interaction to Cox HR
  # prox_pooled_itt_s1_noinx <- glm(Y_neg ~ ns(time_end, knots = c(10,20,30)) + treatment +
  #                             # demographic
  #                             sex_admin + age_years + bmi + race + charlson_cat_fac +
  #                             # other
  #                             ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
  #                             # NEC
  #                             flu_vax,
  #                           data=dat_downsamp.long.itt,
  #                           family=binomial())
  # 
  # dat_downsamp.long.itt$p_itt_noinx <- predict(prox_pooled_itt_s1_noinx, newdata = dat_downsamp.long.itt)
  # 
  # prox_pooled_itt_s2_noinx <- glm(Y_pos ~ ns(time_end, knots = c(10,20,30)) + treatment +
  #                             # demographic
  #                             sex_admin + age_years + bmi + race + charlson_cat_fac +
  #                             # other
  #                             ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
  #                             # predictions from stage 1
  #                             p_itt_noinx,
  #                           # no NEC
  #                           data=dat_downsamp.long.itt,
  #                           family=binomial())
  # 
  # summary(prox_pooled_itt_s2_noinx)
  # exp(-0.379471) # 0.6842233 good :)
  
  ### compare no interaction to Cox HR
  ### if use non intx as models, should run correctly as proportional hazards

# ITT Survival and Risk ---------------------------------------------------

dat_downsamp$gmaxt <- 53

### G formula data setup A=0
prox_itt_A0.long <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$gmaxt),]
prox_itt_A0.long$time_start <- ave(prox_itt_A0.long$fake_mrn, prox_itt_A0.long$fake_mrn, FUN=seq_along)
prox_itt_A0.long$time_start <- (prox_itt_A0.long$time_start-1)*time_unit
prox_itt_A0.long$time_end <- prox_itt_A0.long$time_start+time_unit
prox_itt_A0.long$treatment_obs <- prox_itt_A0.long$treatment
prox_itt_A0.long$treatment <- 0

### G formula data setup A=1
prox_itt_A1.long <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$gmaxt),]
prox_itt_A1.long$time_start <- ave(prox_itt_A1.long$fake_mrn, prox_itt_A1.long$fake_mrn, FUN=seq_along)
prox_itt_A1.long$time_start <- (prox_itt_A1.long$time_start-1)*time_unit
prox_itt_A1.long$time_end <- prox_itt_A1.long$time_start+time_unit
prox_itt_A1.long$treatment_obs <- prox_itt_A1.long$treatment
prox_itt_A1.long$treatment <- 1

### predictions from stage 1 model
prox_itt_A0.long$p_itt <- predict(prox_pooled_itt_s1, newdata=prox_itt_A0.long, type="link") 
prox_itt_A1.long$p_itt <- predict(prox_pooled_itt_s1, newdata=prox_itt_A1.long, type="link") 

### predicted hazards testing NEGATIVE stage 1 model
prox_itt_A0.long$hazard_neg <- predict(prox_pooled_itt_s1, newdata=prox_itt_A0.long, type="response") 
prox_itt_A1.long$hazard_neg <- predict(prox_pooled_itt_s1, newdata=prox_itt_A1.long, type="response")

# ### predicted hazards testing POSITIVE from stage 2 model
  ### overwrite below
# prox_itt_A0.long$hazard_pos <- predict(prox_pooled_itt_s2, newdata=prox_itt_A0.long, type="response")
# prox_itt_A1.long$hazard_pos <- predict(prox_pooled_itt_s2, newdata=prox_itt_A1.long, type="response")

### predicted hazards testing POSITIVE from observed model
prox_itt_A0.long$hazard_pos_obs <- predict(prox_pooled_itt_obs, newdata=prox_itt_A0.long, type="response")
prox_itt_A1.long$hazard_pos_obs <- predict(prox_pooled_itt_obs, newdata=prox_itt_A1.long, type="response")

### new section

time_df <- data.frame(time_end=seq(1,53,1),
                      logHR=haz_ref_A1-haz_ref_A0)

prox_itt_A0.long <- left_join(prox_itt_A0.long, time_df, by="time_end")
prox_itt_A1.long <- left_join(prox_itt_A1.long, time_df, by="time_end")

### switching function

### removing treatment from treated
prox_itt_A0.long$hazard_pos <- prox_itt_A0.long$hazard_pos_obs * exp(-prox_itt_A0.long$logHR * prox_itt_A0.long$treatment_obs)
### adding treated to untreated
prox_itt_A1.long$hazard_pos <- prox_itt_A1.long$hazard_pos_obs * exp(prox_itt_A1.long$logHR * (1-prox_itt_A1.long$treatment_obs))

### calculate (1 - hazard POSITIVE)
prox_itt_A0.long$pnoevent_pos <- 1 - prox_itt_A0.long$hazard_pos
prox_itt_A1.long$pnoevent_pos <- 1 - prox_itt_A1.long$hazard_pos

### calculate (1 - hazard NEGATIVE)
prox_itt_A0.long$pnoevent_neg <- 1 - prox_itt_A0.long$hazard_neg
prox_itt_A1.long$pnoevent_neg <- 1 - prox_itt_A1.long$hazard_neg

### sort the data by ID, time
prox_itt_A0.long <- prox_itt_A0.long[order(prox_itt_A0.long$fake_mrn, prox_itt_A0.long$time_end),] 
prox_itt_A1.long <- prox_itt_A1.long[order(prox_itt_A1.long$fake_mrn, prox_itt_A1.long$time_end),]

### lag (1 - hazard POSITIVE)
prox_itt_A0.long <- prox_itt_A0.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1)) |> 
  ungroup()

prox_itt_A1.long <- prox_itt_A1.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_pos_lag = lag(pnoevent_pos, n=1, default=1)) |> 
  ungroup()

### lag (1 - hazard NEGATIVE)
prox_itt_A0.long <- prox_itt_A0.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_neg_lag = lag(pnoevent_neg, n=1, default=1)) |> 
  ungroup()

prox_itt_A1.long <- prox_itt_A1.long |> 
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |> 
  mutate(pnoevent_neg_lag = lag(pnoevent_neg, n=1, default=1)) |> 
  ungroup()



### AJ estimator

prox_itt_A0.long$surv_prod <- prox_itt_A0.long$pnoevent_neg * prox_itt_A0.long$pnoevent_pos_lag 
prox_itt_A1.long$surv_prod <- prox_itt_A1.long$pnoevent_neg * prox_itt_A1.long$pnoevent_pos_lag

prox_itt_A0.long$survival <- ave(prox_itt_A0.long$surv_prod, prox_itt_A0.long$fake_mrn, FUN=cumprod)
prox_itt_A1.long$survival <- ave(prox_itt_A1.long$surv_prod, prox_itt_A1.long$fake_mrn, FUN=cumprod)

prox_itt_A0.long$risk_prod <- prox_itt_A0.long$hazard_pos * prox_itt_A0.long$survival
prox_itt_A1.long$risk_prod <- prox_itt_A1.long$hazard_pos * prox_itt_A1.long$survival

prox_itt_A0.long$risk_pos <- ave(prox_itt_A0.long$risk_prod, prox_itt_A0.long$fake_mrn, FUN=cumsum)
prox_itt_A1.long$risk_pos <- ave(prox_itt_A1.long$risk_prod, prox_itt_A1.long$fake_mrn, FUN=cumsum)

prox_itt_A0.long.res <- aggregate(risk_pos ~ time_end, data=prox_itt_A0.long, FUN=mean)
prox_itt_A1.long.res <- aggregate(risk_pos ~ time_end, data=prox_itt_A1.long, FUN=mean)

### AJ estimator



# ### Other estimator
# 
# ### product at each time (P(no event neg) * lag P(no event pos))
# prox_itt_A0.long$surv_prod_lag <- prox_itt_A0.long$pnoevent_neg_lag * prox_itt_A0.long$pnoevent_pos_lag # updated w lag neg
# prox_itt_A1.long$surv_prod_lag <- prox_itt_A1.long$pnoevent_neg_lag * prox_itt_A1.long$pnoevent_pos_lag # updated w lag neg
# 
# ### cumulative product within individual
# # prox_itt_A0.long$survival <- ave(prox_itt_A0.long$pnoevent_pos_lag, prox_itt_A0.long$fake_mrn, FUN=cumprod)
# # prox_itt_A1.long$survival <- ave(prox_itt_A1.long$pnoevent_pos_lag, prox_itt_A1.long$fake_mrn, FUN=cumprod)
# prox_itt_A0.long$survival <- ave(prox_itt_A0.long$surv_prod_lag, prox_itt_A0.long$fake_mrn, FUN=cumprod) # keep
# prox_itt_A1.long$survival <- ave(prox_itt_A1.long$surv_prod_lag, prox_itt_A1.long$fake_mrn, FUN=cumprod) # keep
# 
# ### product at each time (haz pos * surv)
# prox_itt_A0.long$risk_prod_pos <- prox_itt_A0.long$hazard_pos * prox_itt_A0.long$survival # keep
# prox_itt_A1.long$risk_prod_pos <- prox_itt_A1.long$hazard_pos * prox_itt_A1.long$survival # keep 
# 
# ### product at each time (haz neg * surv)
# prox_itt_A0.long$risk_prod_neg <- prox_itt_A0.long$hazard_neg * prox_itt_A0.long$survival # keep
# prox_itt_A1.long$risk_prod_neg <- prox_itt_A1.long$hazard_neg * prox_itt_A1.long$survival # keep 
# 
# prox_itt_A0.long.res.pos <- aggregate(cbind(risk_prod_pos, survival) ~ time_end,
#                                   data = prox_itt_A0.long, FUN = mean) # keep
# prox_itt_A1.long.res.pos <- aggregate(cbind(risk_prod_pos, survival) ~ time_end,
#                                   data = prox_itt_A1.long, FUN = mean) # keep
# 
# prox_itt_A0.long.res.neg <- aggregate(cbind(risk_prod_neg, survival) ~ time_end,
#                                   data = prox_itt_A0.long, FUN = mean) # keep
# prox_itt_A1.long.res.neg <- aggregate(cbind(risk_prod_neg, survival) ~ time_end,
#                                   data = prox_itt_A1.long, FUN = mean) # keep
# 
# # lambda T specific to cause
# prox_itt_A0.long.res.pos$lambda_T <- prox_itt_A0.long.res.pos$risk_prod_pos / prox_itt_A0.long.res.pos$survival # keep
# prox_itt_A1.long.res.pos$lambda_T <- prox_itt_A1.long.res.pos$risk_prod_pos / prox_itt_A1.long.res.pos$survival # keep
# 
# # prox_itt_A0.long.res.neg$lambda_T <- prox_itt_A0.long.res.neg$risk_prod_neg / prox_itt_A0.long.res.neg$survival # keep
# # prox_itt_A1.long.res.neg$lambda_T <- prox_itt_A1.long.res.neg$risk_prod_neg / prox_itt_A1.long.res.neg$survival # keep
# 
# prox_itt_A0.long.res.pos$cuminc <- cumsum(prox_itt_A0.long.res.pos$lambda_T * prox_itt_A0.long.res.pos$survival) # keep
# prox_itt_A1.long.res.pos$cuminc <- cumsum(prox_itt_A1.long.res.pos$lambda_T * prox_itt_A1.long.res.pos$survival) # keep
# 
# # ### cumulative sum within individual
# # prox_itt_A0.long$risk <- ave(prox_itt_A0.long$risk_prod, prox_itt_A0.long$fake_mrn, FUN=cumsum)
# # prox_itt_A1.long$risk <- ave(prox_itt_A1.long$risk_prod, prox_itt_A1.long$fake_mrn, FUN=cumsum)
# # 
# # ### calculate the average risk at each time point
# # prox_itt_A0.long.res <- aggregate(risk ~ time_end, data=prox_itt_A0.long, FUN=mean)
# # prox_itt_A1.long.res <- aggregate(risk ~ time_end, data=prox_itt_A1.long, FUN=mean)


### plot the risk curves

# png("results/prox_Y3_risks_itt_p.png", width = 2400, height = 1800, res=300)

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, prox_itt_A0.long.res$time_end, prox_itt_A1.long.res$time_end)),
     ylim = range(c(0, 0.10)),
     xlab="Weeks",
     ylab="Risk",
     main="Risk Curves",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main=1.4
)
mtext("Proximal Causal Inference (ITT)", side = 3, line = 0.5, font = 3, cex=1.2)
grid()
lines(c(0, prox_itt_A0.long.res$time_end), c(0, prox_itt_A0.long.res$risk_pos), col='#006663', lty=1, lwd=4)
lines(c(0, prox_itt_A1.long.res$time_end), c(0, prox_itt_A1.long.res$risk_pos), col='#FF6B1A', lty=1, lwd=4)
legend("topleft",
       legend = c("No Booster", "Booster"),
       col = c("#006663", "#FF6B1A"),
       lty = 1, lwd = 4, bty = "n", cex=1.2)

# dev.off() # 2025-12-10

# rm(prox_itt_A0.long.res, prox_itt_A1.long.res, prox_itt_A0.long, prox_itt_A1.long)

### Plot RR over time

setDT(prox_itt_A0.long.res)
setDT(prox_itt_A1.long.res)

prox_RR_itt <- data.table(
  Week = 1:53,
  RR = sapply(1:53, function(wk) {
    num_pos <- prox_itt_A1.long.res[time_end == wk, risk_pos]
    denom_pos <- prox_itt_A0.long.res[time_end == wk, risk_pos]
    return(
      as.numeric(num_pos) / as.numeric(denom_pos)
    )
  })
)

# write_rds(prox_RR_itt, file = "results/prox_RR_itt.rds") # 2025-12-10

png("results/prox_Y3_RR_itt_p.png", width = 2400, height = 1800, res=300)

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, prox_RR_itt$Week)),
     ylim = range(c(0.0, 1.2)),
     xlab = "Weeks",
     ylab = "Risk Ratio (RR)",
     main = "Risk Ratio Over Time",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main=1.4)
lines(prox_RR_itt$Week, prox_RR_itt$RR, col='black', lty=1, lwd=4)
mtext("Proximal Causal Inference (ITT)", side = 3, line = 0.5, font = 3, cex=1.2)
grid()
abline(h = 1, col = "black", lty = 1, lwd = 0.5)

# dev.off() # 2025-12-10

### Plot EQC and PCI together

# plot(NULL,
#      xlim = range(c(0, eqc_itt_A0.long.res$time_end, eqc_itt_A1.long.res$time_end)),
#      ylim = range(c(0, 0.15)),
#      xlab="Weeks",
#      ylab="Risk",
#      main="Risk Curves"
# )
# mtext("Combined PCI EQC", side = 3, line = 0.5, font = 3)
# lines(c(0, eqc_itt_A0.long.res$time_end), c(0, eqc_itt_A0.long.res$risk_pos), col='gray', lty=2)
# lines(c(0, eqc_itt_A1.long.res$time_end), c(0, eqc_itt_A1.long.res$risk_pos), col='#FF6B1A', lty=2)
# lines(c(0, eqc_itt_A0.long.res.c$time_end), c(0, eqc_itt_A0.long.res.c$risk_pos_c), col='#006663', lty=2)
# lines(c(0, prox_itt_A0.long.res$time_end), c(0, prox_itt_A0.long.res$risk_pos), col='#006663', lty=1)
# lines(c(0, prox_itt_A1.long.res$time_end), c(0, prox_itt_A1.long.res$risk_pos), col='#FF6B1A', lty=1)
# legend("topleft",
#        legend = c("No Booster", "Booster"),
#        col = c("#006663", "#FF6B1A"),
#        lty = 1, lwd = 2, bty = "n")
# legend("topright",
#        legend = c("PCI", "EQC"),
#        lty = c(1, 2), lwd = 2, bty = "n")



# STOP here ---------------------------------------------------------------
# PP effect needs work ----------------------------------------------------




# PP long format expansion ------------------------------------------------

### calc number of rows needed for each individual
time_unit <- 1
dat_downsamp$max_units_pp <- ceiling(dat_downsamp$Y3_pp_t/time_unit)+1
dat_downsamp.long.pp <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$max_units_pp),]

### start and end time corresponding to each row of observation
dat_downsamp.long.pp$time_start <- ave(dat_downsamp.long.pp$fake_mrn, dat_downsamp.long.pp$fake_mrn, FUN=seq_along)
dat_downsamp.long.pp$time_start <- (dat_downsamp.long.pp$time_start-1)*time_unit
dat_downsamp.long.pp$time_end <- dat_downsamp.long.pp$time_start+time_unit

dat_downsamp.long.pp$Y_pos <- ifelse(
  dat_downsamp.long.pp$Y3 == 2 &
    dat_downsamp.long.pp$Y3_pp_t == dat_downsamp.long.pp$time_start,
  1, 0
)

dat_downsamp.long.pp$Y_neg <- ifelse(
  dat_downsamp.long.pp$Y3 == 1 &
    dat_downsamp.long.pp$Y3_pp_t == dat_downsamp.long.pp$time_start,
  1, 0
)

dat_downsamp.long.pp$C <- ifelse(
  dat_downsamp.long.pp$Y3 == 0 &
    dat_downsamp.long.pp$Y3_pp_t == dat_downsamp.long.pp$time_start,
  1, 0
)


dat_downsamp.long.pp$Y_pos <- ifelse(dat_downsamp.long.pp$C==1, NA, dat_downsamp.long.pp$Y_pos)
dat_downsamp.long.pp$Y_neg <- ifelse(dat_downsamp.long.pp$C==1, NA, dat_downsamp.long.pp$Y_neg)


# PP IPCW -----------------------------------------------------------------
# PP Pooled Logistic ------------------------------------------------------

### ipcw stage 1

ipcw1_denom <- glm(C ~ ns(time_end, knots = c(1,15,30,60))*treatment +
                     # demographic
                     sex_admin + age_years + bmi + race + charlson_cat_fac +
                     # other
                     ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                     # NEC
                     flu_vax,
                   data=dat_downsamp.long.pp,
                   family=binomial())

dat_downsamp.long.pp$pd.cens1 <- 1-predict(ipcw1_denom, type = "response", newdata = dat_downsamp.long.pp)

ipcw1_num <- speedglm(C ~ treatment,
                     data = dat_downsamp.long.pp,
                     family = binomial(link = "logit")
)

dat_downsamp.long.pp$pn.cens1 <- 1-predict(ipcw1_num, type = "response", newdata = dat_downsamp.long.pp)

dat_downsamp.long.pp <- dat_downsamp.long.pp |>
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |>
  mutate(
    pdcuml.cens1 = cumprod(pd.cens1),
    pncuml.cens1 = cumprod(pn.cens1)) |>
  mutate(ipcw1 = pncuml.cens1 / pdcuml.cens1) |>
  ungroup()


### pooled logistic stage 1

prox_pooled_pp_s1 <- glm(Y_neg ~ ns(time_end, knots = c(1,15,30,60))*treatment +
                            # demographic
                            sex_admin + age_years + bmi + race + charlson_cat_fac +
                            # other
                            ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                            # NEC
                            flu_vax,
                          data=dat_downsamp.long.pp,
                          family=binomial(), weights = ipcw1)


### preds

dat_downsamp.long.pp$p_pp <- predict(prox_pooled_pp_s1, newdata = dat_downsamp.long.pp)


### ipcw stage 2

ipcw2_denom <- glm(C ~ ns(time_end, knots = c(1,15,30,60))*treatment +
                     # demographic
                     sex_admin + age_years + bmi + race + charlson_cat_fac +
                     # other
                     ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                     # predictions from stage 1
                     p_pp,
                     # no NEC
                   data=dat_downsamp.long.pp,
                   family=binomial())

dat_downsamp.long.pp$pd.cens2 <- 1-predict(ipcw2_denom, type = "response", newdata = dat_downsamp.long.pp)

ipcw2_num <- speedglm(C ~ treatment,
                      data = dat_downsamp.long.pp,
                      family = binomial(link = "logit")
)

dat_downsamp.long.pp$pn.cens2 <- 1-predict(ipcw2_num, type = "response", newdata = dat_downsamp.long.pp)

dat_downsamp.long.pp <- dat_downsamp.long.pp |>
  arrange(fake_mrn, time_end) |> 
  group_by(fake_mrn) |>
  mutate(
    pdcuml.cens2 = cumprod(pd.cens2),
    pncuml.cens2 = cumprod(pn.cens2)) |>
  mutate(ipcw2 = pncuml.cens2 / pdcuml.cens2) |>
  ungroup()

### pooled logistic stage 2

prox_pooled_pp_s2 <- glm(Y_pos ~ ns(time_end, knots = c(1,15,30,60))*treatment +
                            # demographic
                            sex_admin + age_years + bmi + race + charlson_cat_fac +
                            # other
                            ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                            # predictions from stage 1
                            p_pp,
                          # no NEC
                          data=dat_downsamp.long.pp,
                          family=binomial(), weights = ipcw2)


# PP survival and risk ----------------------------------------------------

dat_downsamp$gmaxt <- 79

### G formula data setup A=0
prox_pp_A0.long <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$gmaxt),]
prox_pp_A0.long$time_start <- ave(prox_pp_A0.long$fake_mrn, prox_pp_A0.long$fake_mrn, FUN=seq_along)
prox_pp_A0.long$time_start <- (prox_pp_A0.long$time_start-1)*time_unit
prox_pp_A0.long$time_end <- prox_pp_A0.long$time_start+time_unit
prox_pp_A0.long$treatment <- 0

### G formula data setup A=1
prox_pp_A1.long <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$gmaxt),]
prox_pp_A1.long$time_start <- ave(prox_pp_A1.long$fake_mrn, prox_pp_A1.long$fake_mrn, FUN=seq_along)
prox_pp_A1.long$time_start <- (prox_pp_A1.long$time_start-1)*time_unit
prox_pp_A1.long$time_end <- prox_pp_A1.long$time_start+time_unit
prox_pp_A1.long$treatment <- 1













# cut below ---------------------------------------------------------------
# cut below ---------------------------------------------------------------
# cut below ---------------------------------------------------------------
# 
# 
# # PP survival and risk ----------------------------------------------------
# 
# dat_downsamp$gmaxt <- 79
# 
# ### G formula data setup A=0
# prox_pp_A0.long <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$gmaxt),]
# prox_pp_A0.long$time_start <- ave(prox_pp_A0.long$fake_mrn, prox_pp_A0.long$fake_mrn, FUN=seq_along)
# prox_pp_A0.long$time_start <- (prox_pp_A0.long$time_start-1)*time_unit
# prox_pp_A0.long$time_end <- prox_pp_A0.long$time_start+time_unit
# prox_pp_A0.long$treatment <- 0
# 
# ### G formula data setup A=1
# prox_pp_A1.long <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$gmaxt),]
# prox_pp_A1.long$time_start <- ave(prox_pp_A1.long$fake_mrn, prox_pp_A1.long$fake_mrn, FUN=seq_along)
# prox_pp_A1.long$time_start <- (prox_pp_A1.long$time_start-1)*time_unit
# prox_pp_A1.long$time_end <- prox_pp_A1.long$time_start+time_unit
# prox_pp_A1.long$treatment <- 1
# 
# ### predictions from stage 1 model ("adjusting" for U)
# prox_pp_A0.long$p_pp <- predict(prox_pooled_pp_s1, newdata=prox_pp_A0.long, type="link") # check type
# prox_pp_A1.long$p_pp <- predict(prox_pooled_pp_s1, newdata=prox_pp_A1.long, type="link") # check type
# 
# ### predicted hazards from stage 2 model
# prox_pp_A0.long$hazard <- predict(prox_pooled_pp_s2, newdata=prox_pp_A0.long, type="response") # check type
# prox_pp_A1.long$hazard <- predict(prox_pooled_pp_s2, newdata=prox_pp_A1.long, type="response") # check type
# 
# ### calculate (1 - hazard)
# prox_pp_A0.long$pnoevent <- 1 - prox_pp_A0.long$hazard
# prox_pp_A1.long$pnoevent <- 1 - prox_pp_A1.long$hazard
# 
# ### sort the data by ID, time
# prox_pp_A0.long <- prox_pp_A0.long[order(prox_pp_A0.long$fake_mrn, prox_pp_A0.long$time_end),]
# prox_pp_A1.long <- prox_pp_A1.long[order(prox_pp_A1.long$fake_mrn, prox_pp_A1.long$time_end),]
# 
# ### lag (1 - hazard)
# prox_pp_A0.long <- prox_pp_A0.long |>
#   arrange(fake_mrn, time_end) |>
#   group_by(fake_mrn) |>
#   mutate(pnoevent_lag = lag(pnoevent, n=1, default=1)) |>
#   ungroup()
# 
# prox_pp_A1.long <- prox_pp_A1.long |>
#   arrange(fake_mrn, time_end) |>
#   group_by(fake_mrn) |>
#   mutate(pnoevent_lag = lag(pnoevent, n=1, default=1)) |>
#   ungroup()
# 
# ### cumulative product within individual
# prox_pp_A0.long$survival <- ave(prox_pp_A0.long$pnoevent_lag, prox_pp_A0.long$fake_mrn, FUN=cumprod)
# prox_pp_A1.long$survival <- ave(prox_pp_A1.long$pnoevent_lag, prox_pp_A1.long$fake_mrn, FUN=cumprod)
# 
# ### product at each time (haz * surv)
# prox_pp_A0.long$risk_prod <- prox_pp_A0.long$hazard * prox_pp_A0.long$survival
# prox_pp_A1.long$risk_prod <- prox_pp_A1.long$hazard * prox_pp_A1.long$survival
# 
# ### cumulative sum within individual
# prox_pp_A0.long$risk <- ave(prox_pp_A0.long$risk_prod, prox_pp_A0.long$fake_mrn, FUN=cumsum)
# prox_pp_A1.long$risk <- ave(prox_pp_A1.long$risk_prod, prox_pp_A1.long$fake_mrn, FUN=cumsum)
# 
# ### calculate the average risk at each time point
# prox_pp_A0.long.res <- aggregate(risk ~ time_end, data=prox_pp_A0.long, FUN=mean)
# prox_pp_A1.long.res <- aggregate(risk ~ time_end, data=prox_pp_A1.long, FUN=mean)
# 
# ### plot the risk curves
# 
# # png("results/prox_Y3_risks_pp_p.png", width = 2400, height = 1800, res=300)
# 
# plot(NULL,
#      xlim = range(c(0, prox_pp_A0.long.res$time_end, prox_pp_A1.long.res$time_end)),
#      ylim = range(c(0, 0.25)),
#      xlab="Weeks",
#      ylab="Risk",
#      main="Risk Curves"
# )
# mtext("Per-Protocol (PCI)", side = 3, line = 0.5, font = 3)
# lines(c(0, prox_pp_A0.long.res$time_end), c(0, prox_pp_A0.long.res$risk), col='#006663', lty=1)
# lines(c(0, prox_pp_A1.long.res$time_end), c(0, prox_pp_A1.long.res$risk), col='#FF6B1A', lty=1)
# legend("topleft",
#        legend = c("No Booster", "Booster"),
#        col = c("#006663", "#FF6B1A"),
#        lty = 1, lwd = 2, bty = "n")
# 
# # dev.off() # 2025-10-03
# 
# ### Plot RR over time
# 
# setDT(prox_pp_A0.long.res)
# setDT(prox_pp_A1.long.res)
# 
# prox_RR_pp <- data.table(
#   Week = 1:79,
#   RR = sapply(1:79, function(wk) {
#     num_pos <- prox_pp_A1.long.res[time_end == wk, risk]
#     denom_pos <- prox_pp_A0.long.res[time_end == wk, risk]
#     return(
#       as.numeric(num_pos) / as.numeric(denom_pos)
#     )
#   })
# )
# 
# # png("results/prox_Y3_RR_pp_p.png", width = 2400, height = 1800, res=300)
# 
# plot(NULL,
#      xlim = range(c(0, prox_RR_pp$Week)),
#      ylim = range(c(0.2, 1.5)),
#      xlab = "Weeks",
#      ylab = "Risk Ratio (RR)",
#      main = "Risk Ratio Over Time")
# lines(prox_RR_pp$Week, prox_RR_pp$RR, col='black', lty=1, lwd=2)
# # legend("topright",
# #        legend = c("Corrected", "Original"),
# #        lty = c(1, 2), lwd = 2, bty = "n")
# mtext("Per-Protocol (PCI)", side = 3, line = 0.5, font = 3)
# grid()
# abline(h = 1, col = "black", lty = 1, lwd = 0.5)
# 
# # dev.off() # 2025-10-03












