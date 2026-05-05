##----- Beau Schaeffer
##----- Kaiser TTE-TND
##----- Study period epi curve

library(tidyverse)

# Load data ---------------------------------------------------------------

data_Y3 <- read_rds("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/kaiser/data_weekmatch/data_Y3_weekmatch.rds")

# Prepare data ------------------------------------------------------------

start_date <- as.Date("2022-09-01")

epi_data <- data_Y3 |>
  filter(Y3_itt_trunc == 2) |>
  mutate(
    event_week = index_time + Y3_itt_t_trunc,
    arm = if_else(treatment == 1, "Booster", "No Booster")
  )

weekly_cases <- epi_data |>
  count(event_week, arm) |>
  pivot_wider(names_from = arm, values_from = n, values_fill = 0) |>
  arrange(event_week)

# x-axis labels -----------------------------------------------------------

x_breaks <- seq(0, max(weekly_cases$event_week), by = 4)
x_labels <- format(start_date + x_breaks * 7, "%m/%Y")

# Plot --------------------------------------------------------------------

col0 <- "#006663"
col1 <- "#FF6B1A"

png("figures_draft_wm/epi_curve.png", width = 2400, height = 1800, res = 300)

par(mar = c(5.1, 5.5, 4.1, 2.1))

barplot(
  rbind(weekly_cases$Booster, weekly_cases$`No Booster`),
  col     = c(col1, col0),
  border  = NA,
  space   = 0.15,
  beside  = FALSE,
  axes    = FALSE,
  xlab    = "",
  ylab    = "",
  ylim    = c(0,1000),
  cex.main = 1.4
)

axis(2, cex.axis = 1.5, las = 1)
mtext("COVID-19 Cases (Positive NAAT)", side = 2, line = 4, cex = 1.5)

axis(1,
     at     = (x_breaks + 0.5) * 1.15,
     labels = x_labels,
     cex.axis = 1.1,
     tick   = TRUE)

grid(nx = NA, ny = NULL)

legend("topright",
       legend = c("No Booster", "Booster"),
       fill   = c(col0, col1),
       border = NA,
       bty    = "n",
       cex    = 1.2)

dev.off()


# Greyscale only ----------------------------------------------------------

weekly_cases_total <- epi_data |>
  count(event_week) |>
  arrange(event_week)

# png("figures_draft_wm/epi_curve_grey.png", width = 2400, height = 1800, res = 300)

par(mar = c(5.1, 5.5, 4.1, 2.1))

barplot(
  weekly_cases_total$n,
  col      = "grey40",
  border   = NA,
  space    = 0.15,
  beside   = FALSE,
  axes     = FALSE,
  xlab     = "",
  ylab     = "",
  ylim     = c(0, 1000),
  cex.main = 1.4
)

axis(2, cex.axis = 1.5, las = 1)
mtext("COVID-19 Cases (Positive NAAT)", side = 2, line = 4, cex = 1.5)
axis(1,
     at     = (x_breaks + 0.5) * 1.15,
     labels = x_labels,
     cex.axis = 1.1,
     tick   = TRUE)
grid(nx = NA, ny = NULL)

# dev.off()


