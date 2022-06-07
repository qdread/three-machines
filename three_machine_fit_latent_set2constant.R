# Machine comparison
# Linear, with measurement error in x and y, and including random effects
# Third version that sets both the variance of the latent variable and the slope coefficient on the volumetric machine to be 1
# See responses at https://discourse.mc-stan.org/t/issues-with-convergence-for-latent-variable-model-in-brms/27695/5

### run remotely with the following command
### sbatch --job-name=set2const --ntasks=4 --mem=4gb --export=scriptname="soybean_machine_comparison/three_machine_fit_latent_set2constant.R" runjob.sh


library(data.table)
library(brms)

# Set file path
fp <- ifelse(Sys.info()['sysname'] == 'Windows', 'soybean_machine_comparison/project', '/project/qdr/ars_misc_data')

options(mc.cores = 4, brms.backend = 'cmdstanr', brms.file_refit = 'always')

dat <- fread(file.path(fp, 'Three Machine Anova Individual Measurements.csv'), na.strings = '.')[, 1:8]
dat <- dat[!(Trial %in% 'U7' & Plot %in% 212 & Year %in% 2019)]
dat_wide <- dcast(dat, Trial + Year + Plot + Rep ~ Machine, value.var = 'TW') 

# Latent variable structural equation model ----------------------------------------

dat_wide[, TWtrue := as.numeric(NA)]

bf_gac <- bf(GAC ~ mi(TWtrue))
bf_vol <- bf(Vol ~ mi(TWtrue))
bf_per <- bf(Perten ~  mi(TWtrue))
bf_tw <- bf(TWtrue | mi() ~ 0 + (1|Plot:Trial:Year))

# Set priors such that the intercepts are around the overall mean and the slopes are around 1

overall_means <- dat[, .(TW = mean(TW, na.rm = TRUE)), by = Machine]
stanvars <- stanvar(overall_means$TW[1], 'mean_gac') + stanvar(overall_means$TW[2], 'mean_per') + stanvar(overall_means$TW[3], 'mean_vol')

# Priors assume that the intercept for each machine will be around its mean, and the slope relating any machine to the true value will be roughly 1.
latent_model_priors <- c(
  prior(normal(mean_gac, 1), class = Intercept, resp = GAC),
  prior(normal(mean_per, 1), class = Intercept, resp = Perten),
  prior(normal(mean_vol, 1), class = Intercept, resp = Vol),
  prior(normal(1, 0.25), class = b, resp = GAC),
  prior(normal(1, 0.25), class = b, resp = Perten),
  prior(constant(1), class = b, resp = Vol),
  prior(gamma(5, 5), class = sigma, resp = GAC),
  prior(gamma(5, 5), class = sigma, resp = Perten),
  prior(gamma(5, 5), class = sigma, resp = Vol),
  prior(constant(1), class = sigma, resp = TWtrue)
)

fit_latent_rescor_2constant <- brm(bf_gac + bf_vol + bf_per + bf_tw + set_rescor(TRUE), data = dat_wide,
                  prior = latent_model_priors,
                  chains = 4, iter = 10000, warmup = 7500, seed = 65432, stanvars = stanvars,
                  control = list(adapt_delta = 0.95, max_treedepth = 20),
                  file = file.path(fp, 'three_machine_fit_latent_rescor_2constant'), file_refit = 'never')

fit_latent_norescor_2constant <- brm(bf_gac + bf_vol + bf_per + bf_tw + set_rescor(FALSE), data = dat_wide,
                  prior = latent_model_priors,
                  chains = 4, iter = 10000, warmup = 7500, seed = 23456, stanvars = stanvars,
                  control = list(adapt_delta = 0.95, max_treedepth = 20),
                  file = file.path(fp, 'three_machine_fit_latent_norescor_2constant'), file_refit = 'never')


# Play with output ----------------------------------------------------------

library(tidybayes)
library(dplyr)
library(tidyr)

# norescor model is used

pp_check(fit_latent_norescor_2constant, resp = 'GAC')
pp_check(fit_latent_norescor_2constant, resp = 'Perten')
pp_check(fit_latent_norescor_2constant, resp = 'Vol')

# Get all draws from the model object for the missing values and random effects for plot+trial+year combination
TWtrue_draws <- fit_latent_norescor_2constant %>%
  spread_draws(Ymi_TWtrue[term])

r_draws <- fit_latent_norescor_2constant %>%
  spread_draws(`r_Plot:Trial:Year__TWtrue`[plot_trial_year, term]) %>%
  separate(plot_trial_year, into = c('Plot', 'Trial', 'Year'))

z_draws <- fit_latent_norescor_2constant %>%
  spread_draws(z_1[i,j])

qprobs <- c(.5, .17, .83, .05, .95, .025, .975, .005, .995)

TWtrue_quants <- TWtrue_draws %>%
  summarize(p = qprobs, q = quantile(Ymi_TWtrue, probs = qprobs)) %>%
  pivot_wider(names_from = p, values_from = q, names_prefix = 'q')

# Relationship between observed reading from the model and the "true value"