
## Install and load required packages 
# Function to install packages if not already installed
install_if_missing <- function(packages) {
  for(pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

# List of  required packages 
required_packages <- c(
  "terra",
  "dplyr", 
  "lubridate",
  "ggplot2",
  "scales",
  "pbapply",
  "eurostat",
  "tidyr",
  "sf"
)

# Install and load all packages
install_if_missing(required_packages)

# Load all packages
library(terra)
library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)
library(pbapply)
library(eurostat)
library(tidyr)
library(sf)

cat("All packages successfully installed and loaded!\n")
cat("Packages: terra, dplyr, lubridate, ggplot2, scales, pbapply, eurostat, tidyr, sf\n")

# ----------------------------
# 1. Load NUTS3 polygons
# ----------------------------
nuts0 <- get_eurostat_geospatial(
  output_class = "sf",
  resolution = "60",
  nuts_level = "0",
  year = 2021,
  crs = "4326"
)
nuts0 <- nuts0[nuts0$NUTS_ID != "TR", ]
ext <- ext(-10, 35, 20, 73)  # terra::ext not extent
nuts0 <- st_crop(nuts0, ext)

# ----------------------------
# 2. Thermal-response equations 
# ----------------------------
briere <- function(c, T0, Tm, T) {
  out <- c * T * (T - T0) * sqrt(pmax(0, Tm - T))
  out[T < T0 | T > Tm] <- 0
  out
}

quad <- function(c, T0, Tm, T) {
  out <- c * (T - T0) * (Tm - T)
  out[T < T0 | T > Tm] <- 0
  out
}

# ----------------------------
# 3. Trait functions
# ----------------------------
EFD   <- function(T, c, T0, Tm) briere(c, T0, Tm, T)
pEA   <- function(T, c, T0, Tm) quad(c, T0, Tm, T)
MDR   <- function(T, c, T0, Tm) briere(c, T0, Tm, T)
lf <- function(T, c_lf, T0_lf, Tm_lf) quad(c_lf, T0_lf, Tm_lf, T)
mu <- function(T, c_lf, T0_lf, Tm_lf) {
  x <- lf(T, c_lf, T0_lf, Tm_lf)
  x[x <= 0] <- NA
  1 / x
}

EV    <- function(T) 0.5070 * exp(-((T - 30.85)/12.82)^2)
# ----------------------------
# 4. Parameters
# ----------------------------
### calibrated parameters using priors from Mordecai et al., 2017; 2019
params <- c(
  c_efd = 3.304e-02, T0_efd = 3.412,   Tm_efd = 36.76,
  c_pea = 4.832e-03, T0_pea = 6.309,   Tm_pea = 41.53,
  c_mdr = 8.345e-05, T0_mdr = 4.311,   Tm_mdr = 41.77,
  c_lf  = 2.305e-01, T0_lf  = 10.41,   Tm_lf  = 33.88
)

# ----------------------------
# 5. Load temperature
# ----------------------------
## EG.here we used 2024
this_year = 2024
prev_year = this_year-1
tmax_full = rast("~/Aedes_simulation/t_max.tif") # load maximum temperature
tmax = tmax_full[[which(format(time(tmax_full), "%Y") == this_year)]]
tmax = crop(tmax,nuts0, mask = TRUE)
tmax_1 = tmax_full[[which(format(time(tmax_full), "%Y") == prev_year)]]
tmax_1 = crop(tmax_1,nuts0, mask = TRUE)


tmean = rast("~/Aedes_simulation/tmean.tif")# also load mean temperature
tmean = tmean[[which(format(time(tmean), "%Y") == this_year)]]
tmean = crop(tmean, nuts0, mask = TRUE)
# ----------------------------
# 7. Compute reg_m (ecophysiological limit)
# ----------------------------
cat("Computing ecophysiological regulation...\n")
effect_fun <- function(T, alpha=10, lower_threshold=10.3, upper_threshold=35){
  ifelse(T <= lower_threshold | T >= upper_threshold, 0, 
         1 - exp(-alpha * (T - lower_threshold)^2))
}
##calculate ecophysiological threshold for current year
effect_r <- app(tmax[[1:90]], effect_fun)
mean_effect <- crop(mean(effect_r), nuts0, mask=TRUE)
effect_r_1 <- app(tmax_1[[1:90]], effect_fun)
mean_effect_1 <- crop(mean(effect_r_1), nuts0, mask=TRUE)
ecophy_mean = mean(mean_effect,mean_effect_1)
pal <- c("#2c7bb6","#9FCAD6", "#ffffbf", "#fdae61", "#d7191c")
plot(ecophy_mean,col = pal, breaks = c(0,0.3,0.4,0.6,0.8,1))##plot ecophysiological threshold
plot(st_geometry(nuts0$geometry),add=TRUE)

# Transform to EPSG:3035
ecophy_mean_3035 <- project(mean_effect, "EPSG:3035")
nuts0_3035 <- st_transform(nuts0$geometry, "EPSG:3035")
#classify to boolean
classified_r <- app(ecophy_mean, fun=function(x) ifelse(x>=0.30,1,0))
reg_m <- classified_r
reg_m = resample(reg_m,tmean)
plot(reg_m)


# 8. Mosquito Control effect parameters
#Load SIT control values from data extracted from Balastos et al., 
# ----------------------------

daily_factors <- read.csv("~/Aedes_simulation/SIT_final.csv")
ED <- daily_factors$percentage_density_change
Hr <- daily_factors$percentage.hatched

# Compute Chemical control effect
n_month <- 1:30
cma <- 0.9580 * exp(-0.30 * n_month)
cma_rep <- rep(cma, 4)  
Cc <- rep(0, nlyr(tmean))
Cc[121:240] <- cma_rep[1:120]

cat("Data loaded. Layers:\n")
cat("  Temperature:", nlyr(tmean), "\n")
cat("  ED length:", length(ED), "\n")
cat("  Hr length:", length(Hr), "\n")
cat("  Cc length:", length(Cc), "\n")

# ----------------------------
# Raster implementation 
# ----------------------------
# ----------------------------
# Single time-layer mosquito density (for raster context)
# ----------------------------
mosquito_density_raster <- function(T_vals, params, ED_val, Hr_val, Cc_val, reg_m_vec){
  # Calculate base traits
  efd <- EFD(T_vals, params["c_efd"], params["T0_efd"], params["Tm_efd"])
  pea <- pEA(T_vals, params["c_pea"], params["T0_pea"], params["Tm_pea"])
  mdr <- MDR(T_vals, params["c_mdr"], params["T0_mdr"], params["Tm_mdr"])
  mu_v <- mu(T_vals, params["c_lf"], params["T0_lf"], params["Tm_lf"])
  
  # Apply SIT controls (exact same order as your working code)
  efd_regulated <- efd * ED_val
  efd_regulated <- efd_regulated * (EV(T_vals) * Hr_val)
  
  # Calculate density with chemical control
  out <- (efd_regulated * pea * mdr * reg_m_vec) / ((mu_v + Cc_val)^2)
  
  out[is.na(out)] <- 0
  return(out)
}

# ----------------------------
# Compute across all time layers
# ----------------------------
compute_mosquito_density_fast <- function(T_stack, params, ED, Hr, Cc, reg_m){
  cat("Computing mosquito density...\n")
  
  n_layers <- nlyr(T_stack)
  reg_m_vec <- values(reg_m)[, 1]  # Spatial regulation vector
  
  # Ensure control vectors match number of layers
  if(length(ED) < n_layers){
    ED <- rep(ED, length.out = n_layers)
  } else if(length(ED) > n_layers){
    ED <- ED[1:n_layers]
  }
  
  if(length(Hr) < n_layers){
    Hr <- rep(Hr, length.out = n_layers)
  } else if(length(Hr) > n_layers){
    Hr <- Hr[1:n_layers]
  }
  
  if(length(Cc) < n_layers){
    Cc <- rep(Cc, length.out = n_layers)
  } else if(length(Cc) > n_layers){
    Cc <- Cc[1:n_layers]
  }
  
  out_list <- pblapply(1:n_layers, function(i){
    T_vals <- values(T_stack[[i]])
    
    # Extract scalar control values for this time step
    ED_val <- ED[i]
    Hr_val <- Hr[i]
    Cc_val <- Cc[i]
    
    # Compute density
    dens <- mosquito_density_raster(T_vals, params, ED_val, Hr_val, Cc_val, reg_m_vec)
    # Return as raster
    out_r <- T_stack[[1]]
    values(out_r) <- dens
    out_r
  })
  
  rast(out_list)
}

# ----------------------------
# 12. Compute results
# ----------------------------
#create empty raster for represent simulation without ecophysiologcal threshold
reg_m_none <- tmean[[1]]
values(reg_m_none) <- 1

laea_crs <- "EPSG:3035"
nuts1_laea <- st_transform(nuts0, crs = laea_crs)

# Compute mosquito density with ecophysiology
mosq_stack <- compute_mosquito_density_fast(tmean, params, ED=1, Hr=1, Cc=0, reg_m)
mosq_stack<- project(mosq_stack, laea_crs)
mosq_stack = crop(mosq_stack,nuts1_laea,mask = TRUE)
pal <- c("#2c7bb6","#74add1","#c0e1eb", "#ffffbf", "#fdae61", "#d7191c")
plot(sum(mosq_stack>1000), col = pal, breaks = c(0,30,60,90,120,150,Inf))
plot(nuts1_laea$geometry, add = TRUE)
assign(paste0("mosq_", this_year), mosq_stack)

##compute mosquito density without ecophysiology
mosq_stack_raw <- compute_mosquito_density_fast(tmean, params, ED=1, Hr=1, Cc=0, reg_m_none)
mosq_stack_raw<- project(mosq_stack_raw, laea_crs)
mosq_stack_raw = crop(mosq_stack_raw,nuts1_laea,mask = TRUE)
plot(sum(mosq_stack_raw>1000), col = pal, breaks = c(0,30,60,90,120,150,Inf))
plot(nuts1_laea$geometry, add = TRUE)
assign(paste0("mosq_raw_", this_year), mosq_stack_raw)

#load mosquito occurance
occ <- read.csv("~/Aedes_simulation/occ_aedes.csv")
occ <- occ %>% filter(year == this_year)
# Convert to SpatVector and project
occ_vect <- vect(occ, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
occ_vect <- terra::project(occ_vect, "EPSG:3035")

occ_sf <- st_as_sf(occ_vect)
# Keep only points inside NUTS polygons
occ_sf <- st_intersection(occ_sf, nuts1_laea)

# Plot NUTS1 background (white fill, dark border)
plot(nuts1_laea$geometry,
     col = "white",
     border = "grey20",
     lwd = 1)

# Points with a tiny black outline (pch=21)
points(occ_sf,
       pch = 21,
       bg = "#FF8C00",   # dark orange fill
       col = "black",    # thin boundary
       lwd = 0.3,        # tiny boundary
       cex = 0.5)


####### REPEAT THE ABOVE STEPS YEARLY ###############

####. EVALUATION. ###
## Spatial evaluation:: Evaluate one simulated raster against occurrences at a time
#load simulated mosquito occurances 
raster_eval <- sum(get(paste0("mosq_", this_year)) > 200)
raster_binary <- raster_eval >= 100
# Filter 2024 data
occ_eval <- occ %>% filter(year == this_year)
# Convert to SpatVector and project
occ_vect <- vect(occ_eval, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
occ_vect <- terra::project(occ_vect, "EPSG:3035")
# 1. Extract predicted values at presences
occ_buffer <- buffer(occ_vect, width = 5000)  # buffer 5 km
occ_raster <- rasterize(occ_buffer, raster_eval, field=1)
occ_raster[is.na(occ_raster)] <- 0
# absence = 0
occ_values <- terra::extract(raster_eval, occ_vect)[, 2]
# 1. Apply threshold to raster: values < 100 -> 0
raster_reclass <- terra::ifel(raster_eval <= 100, 0, raster_eval)
pred_raster <- raster_reclass
pred_raster[pred_raster > 0] <- 1
pred_raster[pred_raster == 0] <- 0  # optional, ensures explicit zeros
# Make sure pred_raster has same extent/resolution as occ_raster
pred_raster <- terra::resample(pred_raster, occ_raster, method="near")

# Flatten rasters
obs <- values(occ_raster)
pred <- values(pred_raster)

# Remove NAs
valid_idx <- !is.na(obs) & !is.na(pred)
obs <- obs[valid_idx]
pred <- pred[valid_idx]

# Confusion matrix
TP <- sum(pred == 1 & obs == 1)
FN <- sum(pred == 0 & obs == 1)
FP <- sum(pred == 1 & obs == 0)
TN <- sum(pred == 0 & obs == 0)

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
TSS <- sensitivity + specificity - 1

# --------------------------
# 5. Over- and underprediction fractions
# --------------------------
overprediction_fraction <- FP / (FP + TN + TP + FN)       # fraction of raster predicted but no occurrence
underprediction_fraction <- FN / (TP + FN)               # fraction of occurrences missed

# --------------------------
# 6. ROC / AUC
# --------------------------
# Extract predicted values at occurrences
occ_values <- terra::extract(raster_reclass, occ_vect)[,2]

# Generate background points (pseudo-absences)
set.seed(42)
n_background <- length(occ_values) * 2
bg_points <- spatSample(raster_reclass, size=n_background, method="random", na.rm=TRUE, as.points=TRUE)
bg_values <- terra::extract(raster_reclass, bg_points)[,2]

# Combine for ROC
values_auc <- c(occ_values, bg_values)
labels_auc <- c(rep(1, length(occ_values)), rep(0, length(bg_values)))
library(pROC)
# Compute ROC and AUC
roc_obj <- roc(labels_auc, values_auc, quiet=TRUE)
auc_value <- as.numeric(roc_obj$auc)

# Optional: ROC plot
roc_df <- data.frame(
  specificity = rev(roc_obj$specificities),
  sensitivity = rev(roc_obj$sensitivities)
)

roc_plot <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color="blue", size=1.2) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="gray50") +
  labs(
    x="1 - Specificity",
    y="Sensitivity",
    title=paste("ROC Curve - AUC =", round(auc_value, 3))
  ) +
  theme_minimal()

# --------------------------
# 7. Output all metrics
# --------------------------
evaluation_stats <- list(
  TP = TP,
  FN = FN,
  FP = FP,
  TN = TN,
  Sensitivity = sensitivity,
  Specificity = specificity,
  TSS = TSS,
  Overprediction_fraction = overprediction_fraction,
  Underprediction_fraction = underprediction_fraction,
  AUC = auc_value
  
)
print(evaluation_stats)
#print(roc_plot)

this_year =2024
##Temporal evaluation
year_selected <- this_year

# --------------------------
# 1. Aggregate daily raster to monthly mean
# --------------------------
dates <- seq(as.Date(paste0(this_year, "-01-01")), by = "day", length.out = nlyr(mosq_stack))
assign(paste0("mosq_", this_year), { r <- get(paste0("mosq_", this_year)); time(r) <- dates; r })#change between mosq_raw_/mosq_

year_selected <- this_year
monthly_mean <- tapp(get(paste0("mosq_", this_year)), index = "months", fun = mean, na.rm = TRUE)#change between mosq_raw_/mosq_

# Compute global mean per month
global_monthly_mean <- sapply(1:nlyr(monthly_mean), function(i) {
  mean(values(monthly_mean[[i]]), na.rm = TRUE)
})

# --------------------------
# 2. Summarize occurrences per month
# --------------------------


occ_month_summary <- occ %>%
  filter(year == year_selected) %>%
  mutate(month = month(as.Date(eventDate))) %>%
  filter(!is.na(month)) %>%
  count(month) %>%
  mutate(
    month_name = month.abb[month],
    global_raster_mean = global_monthly_mean[month]
  )


# --------------------------
# 1. Ensure all months 1-12 are present
# --------------------------
all_months <- data.frame(month = 1:12)

occ_month_summary_full <- all_months %>%
  left_join(occ_month_summary, by = "month") %>%
  mutate(
    n = ifelse(is.na(n), 0, n),  # missing months -> 0
    global_raster_mean = ifelse(is.na(global_raster_mean), 0, global_raster_mean),
    month_name = factor(month.abb[month], levels = month.abb)  # order Jan → Dec
  )

# --------------------------
# 2. Rescale observed and predicted
# --------------------------
occ_month_summary_full <- occ_month_summary_full %>%
  mutate(
    n_norm = n / max(n),
    pred_norm = global_raster_mean / max(global_raster_mean)
  )

# --------------------------
# 3. Compute R²
# --------------------------
r_squared <- cor(occ_month_summary_full$n_norm, occ_month_summary_full$pred_norm)^2

# --------------------------
# 4. Plot seasonal comparison
# --------------------------
seasonal_plot <- ggplot(occ_month_summary_full, aes(x = month_name)) +
  geom_col(aes(y = n_norm), fill = "steelblue", alpha = 0.8, width = 0.4, color = "black") +
  geom_line(aes(y = pred_norm, group = 1), color = "red", size = 1.2) +
  geom_point(aes(y = pred_norm), color = "red", shape = 15, size = 3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  
  labs(
    x = "Month",
    y = "Rescaled count (0–1)",
    title = NULL          # removes the title completely
  ) +
  
  # R² top-left
  annotate("text", 
           x = 1, y = 0.98, 
           label = paste0("R² = ", round(r_squared, 3)),
           hjust = 0, vjust = 1, 
           size = 6, fontface = "bold", color = "black") +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 25, face = "bold", margin = margin(t = 15)),
    axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 15)),
    axis.text.x  = element_text(size = 19),
    axis.text.y  = element_text(size = 19),
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.margin  = margin(15, 15, 15, 15)   # gives the big labels a bit of breathing room
  )

print(seasonal_plot)
####REPEAT FOR ALL YEARS ######

### CONTROL IMPLEMENTATION_ WE USED 2024 AS EXAMPLE
# we simulate mosquito density without control, then with chemical and Sterile insect control techniques


mosq_2024_cc <- compute_mosquito_density_fast(tmean, params, ED=1, Hr=1, Cc=Cc, reg_m)###Change Cc from 0
mosq_2024_csit <- compute_mosquito_density_fast(tmean, params, ED=ED, Hr=Hr, Cc=0, reg_m)###Change ED and Hr from 1

mosq_2024_cc <- project(mosq_2024_cc, "EPSG:3035")
mosq_2024_csit <- project(mosq_2024_csit, "EPSG:3035")

##Visualize/Plot
mosq_2024_cc = crop(mosq_2024_cc,nuts1_laea,mask = TRUE)
plot(sum(mosq_2024_cc>1000), col = pal, breaks = c(0,30,60,90,120,150,Inf))
plot(nuts1_laea$geometry, add = TRUE)

mosq_2024_csit = crop(mosq_2024_csit,nuts1_laea,mask = TRUE)
plot(sum(mosq_2024_csit>1000), col = pal, breaks = c(0,30,60,90,120,150,Inf))
plot(nuts1_laea$geometry, add = TRUE)

