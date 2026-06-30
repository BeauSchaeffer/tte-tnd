##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Standard Analysis Pooled
##----- Per-protocol, no censoring weights


# Packages ----------------------------------------------------------------


library(tidyverse)
library(data.table)
library(speedglm)
library(splines)
library(geepack)


# Data --------------------------------------------------------------------


data_Y2 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y2_weekmatch.rds")
dat <- data_Y2

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results_weekmatch_pp/"

# Downsample --------------------------------------------------------------


# subclass_ids <- data_Y2 |> dplyr::select(subclass) |> unique()
# set.seed(345)
# subclass_ids_subset <- slice_sample(subclass_ids, n=150000) # works with 150000
# dat_downsamp <- data_Y2 |> filter(subclass %in% subclass_ids_subset$subclass) |> droplevels()
# rm(subclass_ids, subclass_ids_subset)


# PP long format expansion -----------------------------------------------


### calc number of rows needed for each individual
time_unit <- 1

### ensure at least 1 row for each individual
dat$max_units <- ceiling(dat$Y2_pp_t_trunc/time_unit)+1
dat.long.pp <- dat[rep(1:nrow(dat), dat$max_units),]

### variable that represents the start and end time corresponding to each row of observation
dat.long.pp$time_start <- ave(dat.long.pp$fake_mrn, dat.long.pp$fake_mrn, FUN=seq_along)
dat.long.pp$time_start <- (dat.long.pp$time_start-1)*time_unit
dat.long.pp$time_end <- dat.long.pp$time_start+time_unit

### modify the Y and C variables so that they are only equal to 1 if the 
### event/censoring happened in that time interval
dat.long.pp$Y <- ifelse(
  dat.long.pp$Y2_pp_trunc == 1 &
    dat.long.pp$Y2_pp_t_trunc == dat.long.pp$time_start,
  1, 0
)

dat.long.pp$C <- ifelse(
  dat.long.pp$Y2_pp_trunc == 0 &
    dat.long.pp$Y2_pp_t_trunc == dat.long.pp$time_start,
  1, 0
)

dat.long.pp$Y <- ifelse(dat.long.pp$C==1, NA, dat.long.pp$Y)


# PP Pooled Logistic -----------------------------------------------------


std_pooled_ppp <- speedglm(Y ~ ns(time_end, knots = c(10,20,30))*treatment +
                             # demographic
                             sex_admin + age_years + bmi + race + charlson_cat_fac +
                             # other
                             ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
                             # NEC
                             flu_vax,
                           data=dat.long.pp,
                           family=binomial())

summary(std_pooled_ppp)
# saveRDS(std_pooled_ppp, paste0(res_path,"std_pooled_ppp_model.rds")) # 2026-06-30

### sanity check
### fit same model without interaction terms
### coefficients should be equivalent/similar to Cox

# std_pooled_itt_noint <- speedglm(Y ~ ns(time_end, knots = c(10,20,30)) + treatment +
#                              # demographic
#                              sex_admin + age_years + bmi + race + charlson_cat_fac +
#                              # other
#                              ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks + 
#                              # NEC
#                              flu_vax,
#                            data=dat.long.itt,
#                            family=binomial())
# 
# std_pooled_itt_noint_tidy <- geepack::tidy(std_pooled_itt_noint, conf.int = TRUE, exponentiate = TRUE)
# std_pooled_itt_noint_tidy |> print(n=100)


# ITT Survival and Risk ---------------------------------------------------


dat$gmaxt <- 53

### G formula data setup A=0
std_pp_A0.long <- dat[rep(1:nrow(dat), dat$gmaxt),]
std_pp_A0.long$time_start <- ave(std_pp_A0.long$fake_mrn, std_pp_A0.long$fake_mrn, FUN=seq_along)
std_pp_A0.long$time_start <- (std_pp_A0.long$time_start-1)*time_unit
std_pp_A0.long$time_end <- std_pp_A0.long$time_start+time_unit
std_pp_A0.long$treatment <- 0

### G formula data setup A=1
std_pp_A1.long <- dat[rep(1:nrow(dat), dat$gmaxt),]
std_pp_A1.long$time_start <- ave(std_pp_A1.long$fake_mrn, std_pp_A1.long$fake_mrn, FUN=seq_along)
std_pp_A1.long$time_start <- (std_pp_A1.long$time_start-1)*time_unit
std_pp_A1.long$time_end <- std_pp_A1.long$time_start+time_unit
std_pp_A1.long$treatment <- 1

### Calculate predicted hazards:
std_pp_A0.long$hazard <- predict(std_pooled_pp, newdata=std_pp_A0.long, type="response")
std_pp_A1.long$hazard <- predict(std_pooled_pp, newdata=std_pp_A1.long, type="response")

### Calculate (1 - hazard)
std_pp_A0.long$pnoevent <- 1 - std_pp_A0.long$hazard
std_pp_A1.long$pnoevent <- 1 - std_pp_A1.long$hazard

### Sort the data by time

#* sort by ID, time

std_pp_A0.long <- std_pp_A0.long[order(std_pp_A0.long$time_end),]
std_pp_A1.long <- std_pp_A1.long[order(std_pp_A1.long$time_end),]

# ### Calculate the cumulative survival by taking the cumulative product
# ### of (1 - hazard)
# std_itt_A0$survival <- cumprod(std_itt_A0$pnoevent)
# std_itt_A1$survival <- cumprod(std_itt_A1$pnoevent)

std_pp_A0.long$survival <- ave(std_pp_A0.long$pnoevent, std_pp_A0.long$fake_mrn, FUN=cumprod)
std_pp_A1.long$survival <- ave(std_pp_A1.long$pnoevent, std_pp_A1.long$fake_mrn, FUN=cumprod)

### Calculate risk = 1 - survival
std_pp_A0.long$risk <- 1 - std_pp_A0.long$survival
std_pp_A1.long$risk <- 1 - std_pp_A1.long$survival

# Calculate the average risk at each time point. Here, we're going to use
# the aggregate function to do so:
std_pp_A0.long <- aggregate(risk ~ time_end, data=std_pp_A0.long, FUN=mean)
std_pp_A1.long <- aggregate(risk ~ time_end, data=std_pp_A1.long, FUN=mean)

# save point estimate risk curves in bootstrap-compatible format
std.pp.risk.pointest <- tibble(
  sim = 0L,  # 0 = main analysis (bootstraps are 1..B)
  time_end = std_pp_A0.long$time_end,
  risk0 = std_pp_A0.long$risk,
  risk1 = std_pp_A1.long$risk
)

# saveRDS(std.itt.risk.pointest, paste0(res_path, "std.itt.risk.pointest.rds")) # 2026-04-03

### Plot for the risk curves:

# png(paste0(res_path,"std_Y2_risks_itt_p.png"), width = 2400, height = 1800, res=300)

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, std_itt_A0.long$time_end, std_itt_A1.long$time_end)),
     ylim = range(c(0, 0.1)),
     xlab="Weeks",
     ylab="Risk",
     main="Risk Curves",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main=1.4
)
mtext("Standard TTE (ITT)", side = 3, line = 0.5, font = 3, cex=1.2)
grid()
lines(c(0, std_itt_A0.long$time_end), c(0, std_itt_A0.long$risk), col='#006663', lty=1, lwd=4)
lines(c(0, std_itt_A1.long$time_end), c(0, std_itt_A1.long$risk), col='#FF6B1A', lty=1, lwd=4)
legend("topleft",
       legend = c("No Booster", "Booster"),
       col = c('#006663', '#FF6B1A'),
       lty = 1, lwd = 4, cex=1.2,
       bty = "n")

# dev.off()

### Plot RR over time

setDT(std_itt_A0.long)
setDT(std_itt_A1.long)

std_RR_itt <- data.table(
  Week = 1:53,
  RR = sapply(1:53, function(wk) {
    num <- std_itt_A1.long[time_end == wk, risk]
    denom <- std_itt_A0.long[time_end == wk, risk]
    return(as.numeric(num) / as.numeric(denom))
  })
)

# write_rds(std_RR_itt, file = "results/std_RR_itt.rds") # 2025-12-10


# png("results/std_Y2_RR_itt_p.png", width = 2400, height = 1800, res=300)

par(mar = c(5.1, 5.5, 4.1, 2.1))
plot(NULL,
     xlim = range(c(0, std_RR_itt$Week)),
     # ylim = range(std_RR_itt$RR-0.25, std_RR_itt$RR+0.25),
     ylim = range(c(0.0, 1.2)),
     xlab = "Weeks",
     ylab = "Risk Ratio (RR)",
     main = "Risk Ratio Over Time",
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex.main=1.4)
grid()
lines(std_RR_itt$Week, std_RR_itt$RR, col='black', lty=1, lwd=4)
mtext("Standard TTE (ITT)", side = 3, line = 0.5, font = 3, cex=1.2)
abline(h = 1, col = "black", lty = 1, lwd = 0.5)

# dev.off() # 2025-12-10


# # PP long format expansion ------------------------------------------------
# 
#   ### calc number of rows needed for each individual
# time_unit <- 1
# dat_downsamp$max_units_pp <- ceiling(dat_downsamp$Y2_pp_t_trunc/time_unit)+1
#   # analysis_data_Y_TWO$max_units <- pmax(1, ceiling(analysis_data_Y_TWO$Y2_itt_t / time_unit))
#   ### above line adjusted from causal lab code to ensure that
#   ### at least 1 row created for each individual
# dat_downsamp.long.pp <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$max_units_pp),]
# 
#   ### Now, let's create a variable that's represents the start and end time 
#   ### corresponding to each row of observation
# 
# dat_downsamp.long.pp$time_start <- ave(dat_downsamp.long.pp$fake_mrn, dat_downsamp.long.pp$fake_mrn, FUN=seq_along)
# dat_downsamp.long.pp$time_start <- (dat_downsamp.long.pp$time_start-1)*time_unit
# dat_downsamp.long.pp$time_end <- dat_downsamp.long.pp$time_start+time_unit
# 
#   ### Last, we have to modify the Y and C variables so that they are only 
#   ### equal to 1 if the event/censoring happened in that time interval
# 
# dat_downsamp.long.pp$Y <- ifelse(
#   dat_downsamp.long.pp$Y2_pp_trunc == 1 &
#     dat_downsamp.long.pp$Y2_pp_t_trunc == dat_downsamp.long.pp$time_start,
#   1, 0
# )
# 
# dat_downsamp.long.pp$C <- ifelse(
#   dat_downsamp.long.pp$Y2_pp_trunc == 0 &
#     dat_downsamp.long.pp$Y2_pp_t_trunc == dat_downsamp.long.pp$time_start,
#   1, 0
# )
# 
# dat_downsamp.long.pp$Y <- ifelse(dat_downsamp.long.pp$C==1, NA, dat_downsamp.long.pp$Y)
# 
# # PP IPCW -----------------------------------------------------------------
# 
# dat_downsamp.long.pp <- dat_downsamp.long.pp |> dplyr::arrange(fake_mrn, time_end)
# 
# ipcw_denom <- speedglm(
#   C ~ ns(time_end, knots = c(10,20,30)) + treatment +
#     # demographic
#     sex_admin + age_years + bmi + race + charlson_cat_fac +
#     # other
#     ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
#     # NEC
#     flu_vax,
#   data = dat_downsamp.long.pp, family = binomial(link = "logit")
# )
# 
# dat_downsamp.long.pp$pd.cens <- 1-predict(ipcw_denom, type = "response", newdata = dat_downsamp.long.pp)
# 
# ipcw_num <- speedglm(
#   C ~ treatment,
#   data = dat_downsamp.long.pp,
#   family = binomial(link = "logit")
# )
# 
# dat_downsamp.long.pp$pn.cens <- 1-predict(ipcw_num, type = "response", newdata = dat_downsamp.long.pp)
# 
# dat_downsamp.long.pp <- dat_downsamp.long.pp |>
#   group_by(fake_mrn) |>
#   mutate(
#     pdcuml.cens = cumprod(pd.cens),
#     pncuml.cens = cumprod(pn.cens)) |>
#   mutate(ipcw = pncuml.cens / pdcuml.cens) |>
#   ungroup()
# 
# # dat_downsamp.long.pp$sw.cw <-  dat_downsamp.long.pp$sw * dat_downsamp.long.pp$ipcw ### from IPTW approach
# 
# 
# # PP Pooled Logistic ------------------------------------------------------
# 
# 
# std_pooled_pp <- speedglm(
#   Y ~ ns(time_end, knots = c(10,20,30))*treatment +
#   # demographic
#   sex_admin + age_years + bmi + race + charlson_cat_fac +
#   # other
#   ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
#   # NEC
#   flu_vax,
#   data=dat_downsamp.long.pp,
#   weights = ipcw,
#   family=binomial())
# 
# summary(std_pooled_pp)
# 
# std_pooled_pp_noint <- speedglm(Y ~ ns(time_end, knots = c(10,20,30)) + treatment +
#                                    # demographic
#                                    sex_admin + age_years + bmi + race + charlson_cat_fac +
#                                    # other
#                                    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
#                                    # NEC
#                                    flu_vax,
#                                  data=dat_downsamp.long.pp,
#                                  family=binomial())
# 
# std_pooled_pp_noint_tidy <- geepack::tidy(std_pooled_pp_noint, conf.int = TRUE, exponentiate = TRUE)
# std_pooled_pp_noint_tidy |> print(n=100)
# # readr::write_rds(std_pooled_pp_noint_tidy, file = "results/std_pooled_pp_noint_tidy.rds") # 2025-12-10
# 
# 
# # PP Survival and Risk ----------------------------------------------------
# 
# 
# ### G formula data setup A=0
# std_pp_A0.long <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$gmaxt),]
# std_pp_A0.long$time_start <- ave(std_pp_A0.long$fake_mrn, std_pp_A0.long$fake_mrn, FUN=seq_along)
# std_pp_A0.long$time_start <- (std_pp_A0.long$time_start-1)*time_unit
# std_pp_A0.long$time_end <- std_pp_A0.long$time_start+time_unit
# std_pp_A0.long$treatment <- 0
# 
# ### G formula data setup A=1
# std_pp_A1.long <- dat_downsamp[rep(1:nrow(dat_downsamp), dat_downsamp$gmaxt),]
# std_pp_A1.long$time_start <- ave(std_pp_A1.long$fake_mrn, std_pp_A1.long$fake_mrn, FUN=seq_along)
# std_pp_A1.long$time_start <- (std_pp_A1.long$time_start-1)*time_unit
# std_pp_A1.long$time_end <- std_pp_A1.long$time_start+time_unit
# std_pp_A1.long$treatment <- 1
# 
#   # std_pp_A0 <- data.frame(treatment=0, time_end=unique(dat_downsamp.long.pp$time_end))
#   # std_pp_A1 <- data.frame(treatment=1, time_end=unique(dat_downsamp.long.pp$time_end))
#   # std_pp_A0
#   # std_pp_A1
# 
#   ### Calculate predicted hazards:
# std_pp_A0.long$hazard <- predict(std_pooled_pp, newdata=std_pp_A0.long, type="response")
# std_pp_A1.long$hazard <- predict(std_pooled_pp, newdata=std_pp_A1.long, type="response")
# 
#   ### Calculate (1 - hazard)
# std_pp_A0.long$pnoevent <- 1 - std_pp_A0.long$hazard
# std_pp_A1.long$pnoevent <- 1 - std_pp_A1.long$hazard
# 
#   ### Sort the data by time
# std_pp_A0.long <- std_pp_A0.long[order(std_pp_A0.long$time_end),]
# std_pp_A1.long <- std_pp_A1.long[order(std_pp_A1.long$time_end),]
# 
#   #   ### Calculate the cumulative survival by taking the cumulative product
#   #   ### of (1 - hazard)
#   # std_pp_A0$survival <- cumprod(std_pp_A0$pnoevent)
#   # std_pp_A1$survival <- cumprod(std_pp_A1$pnoevent)
# 
# std_pp_A0.long$survival <- ave(std_pp_A0.long$pnoevent, std_pp_A0.long$fake_mrn, FUN=cumprod)
# std_pp_A1.long$survival <- ave(std_pp_A1.long$pnoevent, std_pp_A1.long$fake_mrn, FUN=cumprod)
# 
#   ### Calculate risk = 1 - survival
# std_pp_A0.long$risk <- 1 - std_pp_A0.long$survival
# std_pp_A1.long$risk <- 1 - std_pp_A1.long$survival
# 
# # Calculate the average risk at each time point. Here, we're going to use
# # the aggregate function to do so:
# std_pp_A0.long <- aggregate(risk ~ time_end, data=std_pp_A0.long, FUN=mean)
# std_pp_A1.long <- aggregate(risk ~ time_end, data=std_pp_A1.long, FUN=mean)
# 
#   ### Plot for the risk curves:
# 
# # png("results/std_Y2_risks_pp_p.png", width = 2400, height = 1800, res=300)
# 
# par(mar = c(5.1, 5.5, 4.1, 2.1))
# plot(NULL,
#      xlim = range(c(0, std_pp_A0.long$time_end, std_pp_A1.long$time_end)),
#      ylim = range(c(0, 0.1)),
#      xlab="Weeks",
#      ylab="Risk",
#      main="Risk Curves",
#      cex.axis = 1.5,
#      cex.lab = 1.5,
#      cex.main=1.4
# )
# mtext("Standard TTE (PP)", side = 3, line = 0.5, font = 3, cex=1.2)
# grid()
# lines(c(0, std_pp_A0.long$time_end), c(0, std_pp_A0.long$risk), col='#006663', lwd=4)
# lines(c(0, std_pp_A1.long$time_end), c(0, std_pp_A1.long$risk), col='#FF6B1A', lwd=4)
# legend("topleft",
#        legend = c("No Booster", "Booster"),
#        col = c('#006663', '#FF6B1A'),
#        lty = c(1, 1), lwd=4,
#        bty = "n", cex=1.2)
# 
# # dev.off() # 2025-12-10
# 
# setDT(std_pp_A0.long)
# setDT(std_pp_A1.long)
# 
# std_RR_pp <- data.table(
#   Week = 1:53,
#   RR = sapply(1:53, function(wk) {
#     num <- std_pp_A1.long[time_end == wk, risk]
#     denom <- std_pp_A0.long[time_end == wk, risk]
#     return(as.numeric(num) / as.numeric(denom))
#   })
# )
# 
# # write_rds(std_RR_pp, file = "results/std_RR_pp.rds") # 2025-12-10
# 
# # png("results/std_Y2_RR_pp_p.png", width = 2400, height = 1800, res=300)
# 
# par(mar = c(5.1, 5.5, 4.1, 2.1))
# plot(NULL,
#      xlim = range(c(0, std_RR_pp$Week)),
#      # ylim = range(std_RR_pp$RR-0.25, std_RR_pp$RR+0.25),
#      ylim = range(c(0, 1.2)),
#      xlab = "Weeks",
#      ylab = "Risk Ratio (RR)",
#      main = "Risk Ratio Over Time",
#      cex.axis = 1.5,
#      cex.lab = 1.5,
#      cex.main=1.4)
# grid()
# lines(std_RR_pp$Week, std_RR_pp$RR, col='black', lty=1, lwd=4)
# mtext("Standard TTE (PP)", side = 3, line = 0.5, font = 3, cex=1.2)
# abline(h = 1, col = "black", lty = 1, lwd = 0.5)
# 
# # dev.off() # 2025-12-10








