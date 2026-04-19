# Multimodel inferevce, following Grueber et al. (2011)
# https://doi.org/10.1111/j.1420-9101.2010.02210.x

library(MuMIn)
options(na.action = "na.fail")

# Fit full model
m_full <- glmer(SR ~ Month + MowFreq +  
                  scale(n_mow_events_befre_sampling) +
                  scale(Litter_Cover) + 
                  scale(road_density_km_per_ha) + 
                  log(patch_size_m2) +
                  #  Bare_Ground_Cover + 
                  scale(slope_degr) + scale(dist_tree) + 
                  scale(sky_view_factor) +
                  scale(Biotop_richness_specific) + 
                  log1p(green_cover_pct) + 
                  (1|PlotNo),
                family = poisson, data = Data_1m2)

# All subsets
all_models <- dredge(m_full)

# View top models
head(all_models)

# Model averaging
avg_model <- model.avg(all_models, subset = delta < 2)
summary(avg_model)



# using Grueber et al. (2011) approach for model averaging  ---------

# Generate all subsets
# all_models <- dredge(m_full)

# Top models (delta AIC < 2)
top_models <- get.models(all_models, subset = delta < 2)

# Model averaging (following Grueber et al.)
avg_model <- model.avg(all_models, subset = delta < 2)

drop1(avg_model, test = "Chisq")

# Full average (shrinkage estimates) - recommended by Grueber et al.
summary(avg_model)
confint(avg_model, full = TRUE)

# Importance of each variable (sum of Akaike weights)
sw(avg_model)



# Create new data for predictions
newdata <- data.frame(
  road_density_km_per_ha = seq(0.069, 0.151, by = 0.001),
  Month = levels(Data_1m2$Month)[1],
  MowFreq = levels(Data_1m2$MowFreq)[1],
  Biotop_richness_specific = mean(Data_1m2$Biotop_richness_specific, na.rm = TRUE),
  slope_degr = mean(Data_1m2$slope_degr, na.rm = TRUE),
  sky_view_factor = mean(Data_1m2$sky_view_factor, na.rm = TRUE),
  dist_tree = mean(Data_1m2$dist_tree, na.rm = TRUE),
  green_cover_pct = mean(Data_1m2$green_cover_pct, na.rm = TRUE),
  patch_size_m2 = mean(Data_1m2$patch_size_m2, na.rm = TRUE),
  PlotNo = Data_1m2$PlotNo[1]  # use first PlotNo as placeholder
)

# Predict with re.form = NA to ignore random effects
pred <- predict(avg_model, newdata = newdata, se.fit = TRUE, type = "response", re.form = NA)

# Create data frame for plotting
pred_df <- data.frame(
  road_density_km_per_ha = newdata$road_density_km_per_ha,
  fit = pred$fit,
  se = pred$se.fit,
  lower = pred$fit - 1.96 * pred$se.fit,
  upper = pred$fit + 1.96 * pred$se.fit
)

# Plot
ggplot(pred_df, aes(x = road_density_km_per_ha, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
  geom_line(linewidth = 1) +
  geom_point(data = Data_1m2,
             aes(road_density_km_per_ha, SR, color = Month),
             pch = 19, size = 1.5, alpha = 0.6) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = "Road density (km/ha)", y = "Species richness")





# Get best model
best_model <- get.models(all_models, subset = 1)[[1]]

# Use effects package
library(effects)
eff <- Effect("road_density_km_per_ha", best_model,
              xlevels = list(road_density_km_per_ha = seq(0.069, 0.151, by = 0.001)))

eff_df <- as.data.frame(eff)

ggplot(eff_df, aes(x = road_density_km_per_ha, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
  geom_line(linewidth = 1) +
  geom_point(data = Data_1m2,
             aes(road_density_km_per_ha, SR, color = Month),
             pch = 19, size = 1.5, alpha = 0.6) +
  theme_bw() +
  scale_color_manual(values = Month_col) +
  labs(x = "Road density (km/ha)", y = "Species richness")
