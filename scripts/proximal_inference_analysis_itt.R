##----- Beau Schaeffer
##----- Kaiser Causal TTE-TND
##----- Proximal Inference ITT Analysis

## in progress

# Packages ----------------------------------------------------------------


library(tidyverse)
library(tidycmprsk)
library(survival)
library(data.table)


# Data --------------------------------------------------------------------


data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_Y3.rds")

data_Y3 <- data_Y3 |> 
  mutate(Y3_itt_factor = case_when(
    Y3_itt_trunc==0 ~ "Censor",
    Y3_itt_trunc==1 ~ "Test Negative",
    Y3_itt_trunc==2 ~ "Test Positive"
  )) |> 
  mutate(Y3_itt_factor = factor(Y3_itt_factor, levels = c("Censor", "Test Negative", "Test Positive")),
         subclass=as.character(subclass))

res_path <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/results/"



# Cox model ---------------------------------------------------------------

# Intention to treat

# Stage 1

prox_itt_s1 <- coxph(
  Surv(Y3_itt_t_trunc, Y3_itt_factor == "Test Negative") ~ treatment +
    # demog
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # NEC
    flu_vax +
    cluster(subclass),
  data = data_Y3
)

# Predictions

data_Y3$p <- predict(prox_itt_s1, newdata = data_Y3)

# Stage 2

prox_itt_s2 <- coxph(
  Surv(Y3_itt_t_trunc, Y3_itt_factor == "Test Positive") ~ treatment +
    # demog
    sex_admin + age_years + bmi + race + charlson_cat_fac +
    # other
    ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
    # predictions from S1
    p +
    # NO NEC
    cluster(subclass),
  data = data_Y3
)

pci.itt.cox.pointest <- tidy(prox_itt_s2, conf.int = TRUE, exponentiate = TRUE)
saveRDS(pci.itt.cox.pointest, paste0(res_path,"pci.itt.cox.pointest.rds")) # 2026-01-12 


# Bootstrap CIs for ITT ---------------------------------------------------

num.boot <- 100

set.seed(1155)
seed <- floor(runif(num.boot)*10^8)

setkey(data_Y3, subclass)
subclasses <- data_Y3[, unique(subclass)]
n_sub <- length(subclasses)

boot.results <- lapply(1:num.boot, function(i){
  
  t0 <- Sys.time()
  
  set.seed(seed[i])
  
  message("Starting PCI ITT Cox bootstrap ", i, " (seed=", seed[i], ")")
  
  # select matched pairs
  samp_sub <- sample(subclasses, size = n_sub, replace = TRUE)
  
  # build boot dataset efficiently via one join:
  # map draw index j -> sampled subclass, then join to replicate all rows per subclass
  map <- data.table(j = seq_along(samp_sub), subclass = samp_sub)
  dat.boot <- data_Y3[map, on = "subclass", allow.cartesian = TRUE]
  # new cluster id per draw
  dat.boot[, subclass_boot := paste0(subclass, ".", j)] 
  
  # run stage 1 model
  prox_itt_s1 <- coxph(
    Surv(Y3_itt_t_trunc, Y3_itt_factor == "Test Negative") ~ treatment +
      # demog
      sex_admin + age_years + bmi + race + charlson_cat_fac +
      # other
      ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
      # NEC
      flu_vax +
      cluster(subclass),
    data = dat.boot
  )
  
  # generate predictions
  
  dat.boot$p <- predict(prox_itt_s1, newdata = dat.boot)
  
  # run stage 2 model
  
  prox_itt_s2 <- coxph(
    Surv(Y3_itt_t_trunc, Y3_itt_factor == "Test Positive") ~ treatment +
      # demog
      sex_admin + age_years + bmi + race + charlson_cat_fac +
      # other
      ndi + prior_inf + tests_count + service_region + last_vax_infect_weeks +
      # predictions from S1
      p +
      # NO NEC
      cluster(subclass),
    data = dat.boot
  )
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message("Finished bootstrap ", i, " in ", round(elapsed, 2), " minutes")
  
  return(c(
    sim=i,
    treatHR=unname(exp(prox_itt_s2$coefficients["treatment"]))
    ))
  
})

boot.long <- bind_rows(lapply(boot.results, tibble::as_tibble_row))
saveRDS(boot.long, paste0(res_path, "pci.itt.cox.boot.long.rds"))





