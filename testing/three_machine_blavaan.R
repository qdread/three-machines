library(lavaan)
library(data.table)
library(semPlot)
library(ggplot2)

# Set file path
fp <- ifelse(Sys.info()['sysname'] == 'Windows', 'project', '/project/qdr/ars_misc_data')

dat <- fread(file.path(fp, 'Three Machine Anova Individual Measurements.csv'), na.strings = '.')[, 1:8]
dat <- dat[!(Trial %in% 'U7' & Plot %in% 212 & Year %in% 2019)]
dat_wide <- dcast(dat, Trial + Year + Plot + Rep ~ Machine, value.var = 'TW') 

dat_wide[, Plot_Trial_Year := interaction(Plot, Trial, Year)]

cfa_model <- 'TWtrue =~ Vol + GAC + Perten' # Because Vol comes first, it will set the variance of it to 1

# # Try with all default priors. This does not work because blavaan does not support multilevel models
# fit_default_priors <- bsem(model = cfa_model, data = dat_wide, cluster = 'Plot_Trial_Year',
#                            target = 'stan', n.chains = 1, burnin = 1000, adapt = 1000, sample = 2000, seed = 919,
#                            mcmcfile = fp, save.lvs = TRUE)

fit_lavaan <- cfa(model = cfa_model, data = dat_wide, cluster = 'Plot_Trial_Year')
semPaths(fit_lavaan)

# Predict latent variables for each combination of 

# Plot the results, mean +/- 95% confidence interval (1.96*SE)

fit_summ <- summary(fit_lavaan)
param_df <- as.data.table(fit_summ$PE)
param_df[, param := fcase(op == "=~", 'slope',
                          op == "~1", 'intercept',
                          op == '~~', 'variance')]

ggplot(param_df[param == "slope"], aes(x = rhs, y = est, ymin = est - qnorm(.975)*se, ymax = est + qnorm(.975)*se)) +
  geom_pointrange() + coord_flip()

ggplot(param_df[param == "intercept" & !lhs %in% "TWtrue"], aes(x = lhs, y = est, ymin = est - qnorm(.975)*se, ymax = est + qnorm(.975)*se)) +
  geom_pointrange() + coord_flip()

ggplot(param_df[param == "variance" & !lhs %in% "TWtrue"], aes(x = lhs, y = est, ymin = est - qnorm(.975)*se, ymax = est + qnorm(.975)*se)) +
  geom_pointrange() + coord_flip()


# Try very wide data ------------------------------------------------------

dat_verywide <- dcast(dat, Trial + Year + Plot ~ Machine + Rep, value.var = "TW")
dat_verywide[, Plot_Trial_Year := interaction(Plot, Trial, Year)]
cfa_model_wide <- 
  'TWtrue =~ Vol + GAC + Perten
   Vol =~ Vol_1 + Vol_2 + Vol_3
   GAC =~ GAC_1 + GAC_2 + GAC_3
   Perten =~ Perten_1 + Perten_2 + Perten_3' 
fit_verywide <- cfa(model = cfa_model_wide, data = dat_verywide)

fit_verywide_summ <- as.data.table(summary(fit_verywide)$PE)

# Here, slope between TWtrue and Vol is set to 1.

lavPredict(fit_verywide, type = 'lv')

dat_verywide_rename <- copy(dat_verywide)
dat_verywide_rename[, Plot_Trial_Year := NULL]
setnames(dat_verywide_rename, c('Trial','Year','Plot','G1','G2','G3','P1','P2','P3','V1','V2','V3'))

cfa_model <- 
  'W_true =~ V + G + P
   V =~ V1 + V2 + V3
   G =~ G1 + G2 + G3
   P =~ P1 + P2 + P3' 

cfa_fit <- cfa(model = cfa_model, data = dat_verywide_rename)
