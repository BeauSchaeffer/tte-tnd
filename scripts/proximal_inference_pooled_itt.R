##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Proximal Inference Analysis Pooled


# Packages ----------------------------------------------------------------


library(tidyverse)
library(data.table)
library(speedglm)
library(splines)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_Y3.rds")
dat <- data_Y3
setDT(dat)

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results/"


# Downsample --------------------------------------------------------------


# subclass_ids <- data_Y3 |> dplyr::select(subclass) |> unique()
# set.seed(345)
# subclass_ids_subset <- dplyr::slice_sample(subclass_ids, n=10000)
# dat_downsamp <- data_Y3 |> dplyr::filter(subclass %in% subclass_ids_subset$subclass) |> droplevels()
# rm(subclass_ids, subclass_ids_subset)


# ITT long format expansion -----------------------------------------------


  ### calc number of rows needed for each individual
time_unit <- 1

### ensure at least 1 row for each individual
dat$max_units <- ceiling(dat$Y3_itt_t_trunc/time_unit)+1
dat.long.itt <- dat[rep(1:nrow(dat), dat$max_units),]

### variable that represents the start and end time corresponding to each row of observation
dat.long.itt$time_start <- ave(dat.long.itt$fake_mrn, dat.long.itt$fake_mrn, FUN=seq_along)
dat.long.itt$time_start <- (dat.long.itt$time_start-1)*time_unit
dat.long.itt$time_end <- dat.long.itt$time_start+time_unit

### modify the Y and C variables so that they are only equal to 1 if the 
  ### event/censoring happened in that time interval
dat.long.itt$Y_pos <- ifelse(
  dat.long.itt$Y3_itt_trunc == 2 &
    dat.long.itt$Y3_itt_t_trunc == dat.long.itt$time_start,
  1, 0
)

dat.long.itt$Y_neg <- ifelse(
  dat.long.itt$Y3_itt_trunc == 1 &
    dat.long.itt$Y3_itt_t_trunc == dat.long.itt$time_start,
  1, 0
)

dat.long.itt$C <- ifelse(
  dat.long.itt$Y3_itt_trunc == 0 &
    dat.long.itt$Y3_itt_t_trunc == dat.long.itt$time_start,
  1, 0
)


dat.long.itt$Y_pos <- ifelse(dat.long.itt$C==1, NA, dat.long.itt$Y_pos)
dat.long.itt$Y_neg <- ifelse(dat.long.itt$C==1, NA, dat.long.itt$Y_neg)


# ITT Pooled Logistic -----------------------------------------------------

  ### time interacting with all variables, note ns()*()
  ### mem pressure ~54-58
  ### crashed, will need to request more than 64gb

prox_pooled_itt_s1 <- glm(Y_neg ~ ns(time_end, knots = c(10,20,30))*(treatment +
                             # demographic
                             sex_admin + age_years + bmi + race + charlson_cat_fac +
                             # other
                             ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                             # NEC
                             flu_vax),
                           data=dat.long.itt,
                           family=binomial())

dat.long.itt$p_itt <- predict(prox_pooled_itt_s1, newdata = dat.long.itt)

prox_pooled_itt_s2 <- glm(Y_pos ~ ns(time_end, knots = c(10,20,30))*treatment +
                            # demographic
                            sex_admin + age_years + bmi + race + charlson_cat_fac +
                            # other
                            ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                            # predictions from stage 1
                            p_itt,
                            # no NEC
                          data=dat.long.itt,
                          family=binomial())

prox_pooled_itt_obs <- glm(Y_pos ~ ns(time_end, knots = c(10,20,30))*treatment +
                             # demographic
                             sex_admin + age_years + bmi + race + charlson_cat_fac +
                             # other
                             ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                             # NEC
                             flu_vax,
                           data=dat.long.itt,
                           family=binomial())



# ITT Survival and Risk ---------------------------------------------------

dat$gmaxt <- 53

### G formula data setup A=0
prox_itt_A0.long <- dat[rep(1:nrow(dat), dat$gmaxt),]
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












