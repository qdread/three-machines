# Machine comparison
# Linear, with measurement error in x and y, and including random effects
# Now the loading on Volumetric is set to 1 to consider it the standard and make the model identifiable
# See responses at https://discourse.mc-stan.org/t/issues-with-convergence-for-latent-variable-model-in-brms/27695/5

### run remotely with the following command
### sbatch --job-name=set1const --ntasks=4 --mem=4gb --export=scriptname="soybean_machine_comparison/three_machine_fit_latent_set1constant.R" runjob.sh


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
  prior(gamma(2, 5), class = sigma, resp = TWtrue)
)

fit_latent_rescor_volconstant <- brm(bf_gac + bf_vol + bf_per + bf_tw + set_rescor(TRUE), data = dat_wide,
                  prior = latent_model_priors,
                  chains = 4, iter = 10000, warmup = 7500, seed = 91919, stanvars = stanvars,
                  control = list(adapt_delta = 0.95, max_treedepth = 20),
                  file = file.path(fp, 'three_machine_fit_latent_rescor_volconstant'))

fit_latent_norescor_volconstant <- brm(bf_gac + bf_vol + bf_per + bf_tw + set_rescor(FALSE), data = dat_wide,
                  prior = latent_model_priors,
                  chains = 4, iter = 10000, warmup = 7500, seed = 19191, stanvars = stanvars,
                  control = list(adapt_delta = 0.95, max_treedepth = 20),
                  file = file.path(fp, 'three_machine_fit_latent_norescor_volconstant'))
