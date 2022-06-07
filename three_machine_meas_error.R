# Machine comparison
# Linear, with measurement error in x and y, and including random effects

library(data.table)
library(ggplot2)
library(lme4)
library(emmeans)
library(multcomp)
library(brms)

options(mc.cores = 4, brms.backend = 'cmdstanr', brms.file_refit = 'on_change')

theme_set(theme_bw() +
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_blank()))

dat <- fread('soybean_machine_comparison/project/Three Machine Anova Individual Measurements.csv', na.strings = '.')[, 1:8]
dat <- dat[!(Trial %in% 'U7' & Plot %in% 212 & Year %in% 2019)]
dat_wide <- dcast(dat, Trial + Year + Plot + Rep ~ Machine, value.var = 'TW') 

# Get mean and SD for each of the sets of three samples

dat_by_plot <- dat[, .(TW_mean = mean(TW, na.rm = TRUE), TW_sd = sd(TW, na.rm = TRUE)), by = .(Trial, Machine, Year, Plot)]
dat_by_plot_wide <- dcast(dat_by_plot, Trial + Year + Plot ~ Machine, value.var = c('TW_mean', 'TW_sd'))

# Fit measurement error model with proper random effects structure
lmer(TW ~ Machine + (1|Year) + (1|Trial) + (1|Plot:Trial), data = dat)

gac_vol_formula <- bf(TW_mean_GAC ~ me(TW_mean_Vol, TW_sd_Vol) + (1 | Year) + (1 | Trial) + (1|Plot:Trial:Year))

fit_measerror_gac_vol <- brm(gac_vol_formula, 
                             data = dat_by_plot_wide, 
                             family = gaussian(link = 'log'),
                             prior = c(prior(normal(0, 2), class = b)),
                             chains = 4, iter = 3000, warmup = 2000, seed = 524,
                             file = 'soybean_machine_comparison/project/measerror_gac_vol')

fit_measerror_all <- brm(TW ~ Machine + (1 | Trial) + (1 | Year) + (1 | Plot:Trial:Year))

# Hierarchical model formula
# Just random intercept for the single sample (plot nested within trial nested within year)
gac_vol_formula <- bf(GACtrue ~ Voltrue)
gac_formula <- bf(GACtrue ~ GAC + (1|Plot:Trial:Year))
vol_formula <- bf(Voltrue ~ Vol + (1|Plot:Trial:Year))

fit_measerror_gac_vol <- brm(gac_vol_formula + gac_formula + vol_formula,
                             data = dat_wide, 
                             prior = c(prior(normal(0, 2), class = b)),
                             chains = 4, iter = 3000, warmup = 2000, seed = 524,
                             file = 'soybean_machine_comparison/project/measerror_gac_vol')


# Meas error model, linear, no R.E. ---------------------------------------

gac_vol_ME_formula_no_RE <- bf(
  TW_mean_GAC | resp_se(TW_sd_GAC, sigma = TRUE) ~ me(TW_mean_Vol, TW_sd_Vol)
)

gac_vol_ME_fit_no_RE <- brm(gac_vol_ME_formula_no_RE, 
                            data = dat_by_plot_wide, 
                            chains = 4, iter = 3000, warmup = 2000, seed = 524,
                            file = 'soybean_machine_comparison/project/gac_vol_ME_fit_no_RE')
# Fails due to zero measurement error for some of the values.


# Overall model -----------------------------------------------------------


overall_model_byplot <- lmer(TW ~ Machine + (1|Plot:Trial:Year), data = dat)

# brm model where sigma varies by machine

fixedsigma_model <- bf(TW ~ Machine + (1|Plot:Trial:Year))
variablesigma_model <-  bf(TW ~ Machine + (1|Plot:Trial:Year), sigma ~ Machine)

fit_fixedsigma <- brm(fixedsigma_model,
                      data = dat, 
                      prior = c(prior(normal(0, 2), class = b),
                                prior(normal(56.4, 2), class = Intercept)),
                      chains = 4, iter = 5000, warmup = 4000, seed = 524,
                      file = 'soybean_machine_comparison/project/fit_fixedsigma')

fit_variablesigma <- brm(variablesigma_model,
                         data = dat, 
                         prior = c(prior(normal(0, 2), class = b),
                                   prior(exponential(1), class = sd),
                                   prior(normal(56.4, 2), class = Intercept)),
                         chains = 4, iter = 5000, warmup = 4000, seed = 524,
                         file = 'soybean_machine_comparison/project/fit_variablesigma')

fit_fixedsigma <- add_criterion(fit_fixedsigma, 'loo')
fit_variablesigma <- add_criterion(fit_variablesigma, 'loo')

loo_compare(fit_fixedsigma, fit_variablesigma)


# Try at structural equation model ----------------------------------------

dat_wide[, TWtrue := as.numeric(NA)]

bf_gac <- bf(GAC ~ mi(TWtrue))
bf_vol <- bf(Vol ~ mi(TWtrue))
bf_per <- bf(Perten ~ mi(TWtrue))
bf_tw <- bf(TWtrue | mi() ~ (1|Plot:Trial:Year))

# Set priors such that the intercepts are around the overall mean and the slopes are around 1

overall_means <- dat[, .(TW = mean(TW,na.rm=T)), by=Machine]
stanvars <- stanvar(overall_means$TW[1], 'mean_gac') + stanvar(overall_means$TW[2], 'mean_per') + stanvar(overall_means$TW[3], 'mean_vol')

fit_latent <- brm(bf_gac + bf_vol + bf_per + bf_tw + set_rescor(TRUE), data = dat_wide,
                  prior = c(
                    prior(normal(mean_gac, 2), class = Intercept, resp = GAC),
                    prior(normal(mean_per, 2), class = Intercept, resp = Perten),
                    prior(normal(mean_vol, 2), class = Intercept, resp = Vol),
                    prior(normal(0, 2), class = Intercept, resp = TWtrue),
                    prior(normal(1, 1), class = b, resp = GAC),
                    prior(normal(1, 1), class = b, resp = Perten),
                    prior(normal(1, 1), class = b, resp = Vol)
                  ),
                  chains = 4, iter = 5000, warmup = 4000, seed = 525, stanvars = stanvars,
                  control = list(adapt_delta = 0.9, max_treedepth = 20),
                  file = 'soybean_machine_comparison/project/fit_latent')
